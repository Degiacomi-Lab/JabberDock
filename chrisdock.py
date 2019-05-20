import scipy.ndimage.morphology as SNM
import numpy as np
import time
import matplotlib.pyplot as plt
import skfmm
import torch
import torch.nn as nn
import visdom

vis = visdom.Visdom(port=8613)

def dirac(x):
    return (1/np.pi)/(1+x**2)

def heaviside(x,epsilon):
    return 0.5*(1+(2./np.pi)*torch.atan(x/epsilon))

def dirac_eps(x, epsilon):
    return (epsilon/np.pi)/(epsilon**2+x**2)

def hacknrm(x):
    return ((x-x.min())/(x.max()-x.min()))

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


def grad_mag_3d(x):
    base = np.array([[1, 0, -1], [2, 0, -2], [1, 0, -1]])
    base = torch.from_numpy(base).float().unsqueeze(0)
    base = torch.cat((base, torch.zeros(1, 3, 3), -base), dim=0)

    gx_weights = base.permute(0, 1, 2).unsqueeze(0).unsqueeze(0).to(device=device)
    gy_weights = base.permute(2, 0, 1).unsqueeze(0).unsqueeze(0).to(device=device)
    gz_weights = base.permute(1, 2, 0).unsqueeze(0).unsqueeze(0).to(device=device)

    gx = torch.nn.functional.conv3d(x, gx_weights, bias=None, stride=1, padding=0, dilation=1) / 8.0
    gy = torch.nn.functional.conv3d(x, gy_weights, bias=None, stride=1, padding=0, dilation=1) / 8.0
    gz = torch.nn.functional.conv3d(x, gz_weights, bias=None, stride=1, padding=0, dilation=1) / 8.0

    gm = torch.sqrt(gx**2 + gy**2 + gz**2)
    return gm

def grad_mag_2d(x):
    base = np.array([[1, 0, -1], [2, 0, -2], [1, 0, -1]])
    base = torch.from_numpy(base).float()

    gx_weights = base.permute(0, 1).unsqueeze(0).unsqueeze(0).repeat(x.size(1), x.size(1), 1, 1).to(device=device)
    gy_weights = base.permute(1, 0).unsqueeze(0).unsqueeze(0).repeat(x.size(1), x.size(1), 1, 1).to(device=device)

    # PADDING ERROR!?
    gx = torch.nn.functional.conv2d(x, gx_weights, bias=None, stride=1) / 8.0
    gy = torch.nn.functional.conv2d(x, gy_weights, bias=None, stride=1) / 8.0

    gm = torch.sqrt(gx**2 + gy**2)
    return gm

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


 # todo: initialise randomly, particle-swarm blah blah

# theta[0,:,3] = 0.5

# grid = torch.nn.functional.affine_grid(theta, phi_2.size())

# tmp1 = torch.nn.functional.grid_sample(phi_1, grid, mode='bilinear', padding_mode='border')
# vis.image(hacknrm(tmp1[0,0,23]), win='blah2')


# learning algorithm

# 0) 
# precompute sdf of two proteins
phi_1 = torch.from_numpy(skfmm.distance(receptor-0.43)).float().to(device).unsqueeze(0).unsqueeze(0)
phi_2 = torch.from_numpy(skfmm.distance(ligand-0.43)).float().to(device).unsqueeze(0).unsqueeze(0)

# precompute gradients of phi_1, phi_2
grad_phi_1 = sobel(phi_1.squeeze(0))
grad_phi_2 = -1*sobel(phi_2.squeeze(0))

# padding just spatial dimension (1,1), (1,1), (1,1)
grad_phi_1 = torch.nn.functional.pad(grad_phi_1, pad=(1,1,1,1,1,1), mode='replicate').data
grad_phi_2 = torch.nn.functional.pad(grad_phi_2, pad=(1,1,1,1,1,1), mode='replicate').data

# precompute integral
surface_area_2 = dirac(phi_2).sum()

# 1) initialise theta
theta = torch.eye(3,4).unsqueeze(0).to(device)
theta_eye = torch.eye(3,4).unsqueeze(0).to(device)
theta_eye[0,:,-1] = 1
theta = nn.Parameter(theta)

optim = torch.optim.Adam([theta], lr=0.001)

vis.line(X=np.array([0]), Y=np.array([[np.nan, np.nan, np.nan]]),
    win='loss', opts=dict(title='loss', legend=['1', '2', '3']))

for step in range(500):

    optim.zero_grad()

    # sample domain
    grid = torch.nn.functional.affine_grid(theta*theta_eye, phi_2.size())
    s_phi_1 = torch.nn.functional.grid_sample(phi_1, grid, mode='bilinear', padding_mode='border')
    s_grad_phi_1 = torch.nn.functional.grid_sample(grad_phi_1, grid, mode='bilinear', padding_mode='border')

    mask = dirac(phi_2)
    align = ((s_grad_phi_1*grad_phi_2).sum(dim=1).unsqueeze(1))
    dist = (s_phi_1 - phi_2)**2
    contact = dirac_eps(s_phi_1, 1.6)/surface_area_2

    loss = (mask * align * torch.exp(-0.5*dist) * contact).mean()
    # loss = align.mean()
    loss.backward()
    optim.step()

    vis.image(hacknrm(mask   [0,0,21]), win='m', opts={'title':'mask'})
    vis.image(hacknrm(align  [0,0,21]), win='a', opts={'title':'align'})
    vis.image(hacknrm(dist   [0,0,21]), win='d', opts={'title':'dist'})
    vis.image(hacknrm(contact[0,0,21]), win='c', opts={'title':'contact'})
    vis.line(X=np.array([step]), Y=np.array([[align.mean().item(), dist.mean().item(), contact.mean().item()]]),
        win='loss', opts=dict(legend=['align', 'dist', 'contact']), update='append')



# rec_data = torch.from_numpy(receptor).to(device).float()
# d_rec3 = sobel.forward(rec_data)
# rec_result =  d_rec3.data.cpu().numpy()

# lig_data = torch.from_numpy(ligand).to(device).float()
# d_lig3 = sobel.forward(lig_data)
# lig_result =  d_rec3.data.cpu().numpy()




# newSDF_A = torch.nn.grid_sample(sdfA, thetaA)



# '''
# fig = plt.figure()
# ax = fig.add_subplot(1,3,1)
# plt.imshow(d_rec1[45])
# plt.colorbar()
# ax = fig.add_subplot(1,3,2)
# plt.imshow(d_rec2[45])
# plt.colorbar()
# ax = fig.add_subplot(1,3,3)
# plt.imshow(result[:,45])
# plt.colorbar()
# '''

# fig = plt.figure()
# plt.imshow(rec_result[:,:,45])
# plt.colorbar()

# plt.show()

