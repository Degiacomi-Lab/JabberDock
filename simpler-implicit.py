import scipy.ndimage.morphology as SNM
import numpy as np
import time
import skfmm
import torch
import torch.nn as nn
import visdom

vis = visdom.Visdom(port=8345)

def dirac(x):
    return (1/np.pi)/(1+x**2)

def heaviside(x,epsilon=1):
    return 0.5*(1+(2./np.pi)*torch.atan(x/epsilon))

def dirac_eps(x, epsilon=1):
    return (epsilon/np.pi)/(epsilon**2+x**2)

def hacknrm(x):
    return ((x-x.min())/(x.max()-x.min()))

def quat2mat(quat):
    # see https://raw.githubusercontent.com/ClementPinard/SfmLearner-Pytorch/master/inverse_warp.py
    norm_quat = torch.cat([quat[:,:1].detach()*0 + 1, quat], dim=1)
    norm_quat = norm_quat/norm_quat.norm(p=2, dim=1, keepdim=True)
    w, x, y, z = norm_quat[:,0], norm_quat[:,1], norm_quat[:,2], norm_quat[:,3]

    B = quat.size(0)

    w2, x2, y2, z2 = w.pow(2), x.pow(2), y.pow(2), z.pow(2)
    wx, wy, wz = w*x, w*y, w*z
    xy, xz, yz = x*y, x*z, y*z

    rotMat = torch.stack([w2 + x2 - y2 - z2, 2*xy - 2*wz, 2*wy + 2*xz,
                          2*wz + 2*xy, w2 - x2 + y2 - z2, 2*yz - 2*wx,
                          2*xz - 2*wy, 2*wx + 2*yz, w2 - x2 - y2 + z2], dim=1).reshape(B, 3, 3)
    return rotMat

def visualise(a,b, win='v', normalise=True, sf=1.0, pre_f=lambda x: x, f=lambda x: x):
    # a,b = signed distance functions
    # sf = scale_factor (e.g. 8 = 8x bigger image)
    # f  = post_processing function (e.g. heaviside)
    slice_a = f(torch.nn.functional.interpolate(pre_f(a), scale_factor=sf, mode='trilinear', align_corners=True)[0,0,sf*a.size(2)//2,:,:])
    slice_b = f(torch.nn.functional.interpolate(pre_f(b), scale_factor=sf, mode='trilinear', align_corners=True)[0,0,sf*b.size(2)//2,:,:])
    img = torch.cat((slice_a.unsqueeze(0),slice_b.unsqueeze(0),slice_b.unsqueeze(0)), dim=0)
    vis.image(hacknrm(img) if normalise else img, win=win)

# helper function for sobel orientation generation
def circular(iterable):
    count = len(iterable)
    combi = list(iterable) + list(iterable)
    for i in range(count):
        yield combi[i:i + count]

def sobel_filts(dims, normalised=True):
    # Sobel-Feldman Operator (n dimensional) consists of two operations:
    # Triangle-filter smoothing (1,2,1) in directions perpendicular to derivative,
    # Central difference in derivative direction (-1,0,1)
    # elementwise multiplication of these ops generates filter
    tri_filt = torch.tensor([1.0, 2.0, 1.0])
    c_diff = torch.tensor([-1.0, 0.0, 1.0])
    if normalised:
        tri_filt /= 4
        c_diff /= 2

    # Output shape is going to be a 3x3x3x... filter
    out_size = [3] * dims

    tri_filts = list()

    for i in range(dims - 1):
        view_shape = [1] * dims
        view_shape[i] = 3
        tri_filts.append(tri_filt.view(view_shape).expand(out_size))

    # Assume derivative is occuring along last dimension for initial tensor
    view_shape = [1] * dims
    view_shape[-1] = 3
    c_diff = c_diff.view(view_shape).expand(out_size)

    # Multiply filters together to generate sobel filter
    out_filt = torch.ones(out_size)
    for filt in tri_filts:
        out_filt *= filt
    out_filt *= c_diff
    # torch tensor (num_filters, filter shape)
    # 2D: [2x3x3]
    # 3D: [3x3x3x3]...
    filtlist = torch.empty([dims, *[3] * dims])

    for i, perm in enumerate(circular(range(dims))):
        filtlist[i] = out_filt.permute(*perm)
    return filtlist

# Define these as classes so we don't have to re-initialise the sobel filter tensor on every iteration.
class Sobel3d(nn.Module):
    def __init__(self, magnitude=False):
        super(Sobel3d, self).__init__()
        self.magnitude = magnitude
        filters = sobel_filts(3).reshape([3, 1, 3, 3, 3])
        self.conv = torch.nn.Conv3d(
            in_channels=1,
            out_channels=3,
            kernel_size=3,
            bias=False,
            groups=1)
        print(self.conv.weight.shape)
        for i in range(self.conv.weight.shape[0]):
            self.conv.weight.data[i] = filters[i % 3]

    def forward(self, x):
        inshape = x.shape
        outshape = (*inshape[:-3], 3, *(x-2 for x in inshape[-3:]))
        out = self.conv(x.reshape((-1, 1, *inshape[-3:]))).reshape(outshape)
        if self.magnitude:
            out = out.norm(dim=-4)
        return out

def gaussian_kernel(sigma, kernel_size, dim):
    if (kernel_size % 2 == 0):
        kernel_size += 1

    x = torch.arange(kernel_size).float().to(device)
    kernel_size = len(x)
    mean = (kernel_size - 1)/2.
    variance = sigma**2.

    x = (1./(2.*math.pi*variance)) * torch.exp(-(x-mean)**2 / (2*variance))
    x /= torch.sum(x)

    if dim==2:
        gx = x.unsqueeze(0).repeat(kernel_size, 1)
        gy = x.unsqueeze(1).repeat(1, kernel_size)
        g = gx*gy
        return g
    else:
        gx = x.unsqueeze(0).unsqueeze(0).repeat(kernel_size, kernel_size, 1)
        gy = x.unsqueeze(1).unsqueeze(0).repeat(kernel_size, 1, kernel_size)
        gz = x.unsqueeze(1).unsqueeze(2).repeat(1, kernel_size, kernel_size)
        g = gx*gy*gz
        return g

def gaussian_blur(x, sigma, kernel_size=17, dim=2):
    g = gaussian_kernel(sigma, kernel_size, dim)
    c = torch.nn.functional.conv3d if dim == 3 else torch.nn.functional.conv2d
    x = c(x, g.unsqueeze(0).unsqueeze(0), bias=None, stride=1, padding=kernel_size//2)
    return x

#####################################################################################

device = torch.device('cuda') if torch.cuda.is_available() else torch.device('cpu')
sobel = Sobel3d().to(device)

receptor = np.load("receptor.npy")
ligand = np.load("good_ligand.npy")

# 0) 
# precompute sdf of two proteins (negative on inside)
phi_1 = -torch.from_numpy(skfmm.distance(receptor-0.43)).float().to(device).unsqueeze(0).unsqueeze(0)
phi_2 = -torch.from_numpy(skfmm.distance(ligand-0.43)).float().to(device).unsqueeze(0).unsqueeze(0)

# precompute gradients of phi_1, phi_2
grad_phi_1 = sobel(phi_1.squeeze(0))
grad_phi_2 = -1*sobel(phi_2.squeeze(0))

# padding just spatial dimension (1,1), (1,1), (1,1)
grad_phi_1 = torch.nn.functional.pad(grad_phi_1, pad=(1,1,1,1,1,1), mode='replicate').data
grad_phi_2 = torch.nn.functional.pad(grad_phi_2, pad=(1,1,1,1,1,1), mode='replicate').data

# precompute integral
surface_area_2 = dirac(phi_2).sum()

# 1) initialise rotate/translate parameterisation
quat  = geoopt.ManifoldParameter(torch.randn(1,4).to(device), manifold=geoopt.manifolds.Sphere().to(device)) # q1..4
quat.proj_()
trans = geoopt.ManifoldParameter(torch.randn(1,3,1).to(device)) # xyz 
optim = geoopt.optim.RiemannianSGD([quat,trans], lr=0.1) # parameters to optimise and the learning rate


vis.line(X=np.array([0]), Y=np.array([[np.nan, np.nan, np.nan]]),
    win='loss', opts=dict(title='loss'))

def build_matrix():
    return torch.cat((quat2mat(quat), trans), dim=2)

#parameters to play with
contact_c = 1.
intersect_c = 1.


zero_intersect = torch.zeros(torch.max(phi_2).size()).to(device)
for step in range(1000):

    optim.zero_grad()

    grid = torch.nn.functional.affine_grid(build_matrix(), phi_2.size())
    s_phi_1 = torch.nn.functional.grid_sample(phi_1, grid, mode='bilinear', padding_mode='border')
    
    # Contact
    contact = -heaviside(torch.min(s_phi_1, phi_2))#dirac(torch.min(s_phi_1, phi_2))

    # Intersection
    # Remember that the property data values now just contain the distance (essentially) to the surface.
    # Hence why we use max and min here. The heaviside penalises distances far from the contact / 
    # intersection site (which is what we want)
    k = torch.max(zero_intersect.data, -torch.max(s_phi_1, phi_2))
    intersection = heaviside(k) #-heaviside(torch.max(s_phi_1, phi_2), 1) 
    distance = torch.sum(trans**2, dim=1)
    loss = contact_c * contact.mean() + intersect_c * intersection.mean() + distance

    loss.backward()
    optim.step()

    if step % 100 == 0:
        # vis.image((intersection   [0,0,21]), win='inter', opts={'title':'intersection'})
        # vis.image((hacknrm(contact)        [0,0,21]), win='conthetatact', opts={'title':'contact'})
        visualise(s_phi_1, phi_2, sf=8, f=dirac)
        vis.line(X=np.array([step]), Y=np.array([[intersection.mean().item(), contact.mean().item(), distance.mean().item()]]),
            win='loss', opts=dict(legend=['intersection', 'contact', 'distance']), update='append')


