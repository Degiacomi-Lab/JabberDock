import scipy.ndimage.morphology as SNM
import numpy as np
import time
import skfmm
import torch
import torch.nn as nn
import visdom
import geoopt
import biobox as bb
from copy import deepcopy

vis = visdom.Visdom(port=8345)

#### Define necessary functions #####

def dirac(x):
    return (1/np.pi)/(1+x**2)

def heaviside(x,epsilon=1):
    return 0.5*(1+(2./np.pi)*torch.atan(x/epsilon))

def dirac_eps(x, epsilon=1):
    return (epsilon/np.pi)/(epsilon**2+x**2)

def hacknrm(x):
    return ((x-x.min())/(x.max()-x.min()))

def quat2mat(quat):
    a, b, c, d = quat[0][0], quat[0][1], quat[0][2], quat[0][3]

    a2, b2, c2, d2 = a.pow(2), b.pow(2), c.pow(2), d.pow(2)

    rotMat = torch.stack([a2 + b2 - c2 - d2, 2*b*c - 2*a*d, 2*b*d + 2*a*c,
                         2*b*c + 2*a*d, a2 - b2 + c2 - d2, 2*c*d - 2*a*b,
                         2*b*d - 2*a*c, 2*c*d + 2*a*b, a2 - b2 - c2 + d2], dim=0).reshape(3,3)

    return rotMat.unsqueeze(0)

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

def build_matrix(quat, trans):
    return torch.cat((quat2mat(quat), trans), dim=2)

def fit_model(phi_1, phi_2, contact_c, intersect_c, distance_c, Q, T, steps=1000):
    """
    Try and find a solution for the phi_1 docking into phi_2. Initialise with a random
    translation and quaternion. Constants are supplied externally (this is used to 
    optmise these parameters for specific solutions - e.g. guiding toward getting our complex)

    :param phi_1: signed distance function for protein 1
    :param phi_2: signed distance function for protein 2
    :param contact_c: contact score constant
    :param intersect_c: intersection score constant
    :param distance_c: distance score constant
    :param Q: Quaternion seed (in torch.tensor form of (1,4) shape)
    :param T: Translation seed (in torch.tensor form of (1, 3, 1) shape)
    :param steps: number of steps to run during optimisation
    :returns: quaternion and translation used to find the optimum orientation to satify the loss function
    """

    # initialise rotate/translate parameterisation
    quat  = geoopt.ManifoldParameter(Q.to(device), manifold=geoopt.manifolds.Sphere().to(device)) # q1..4
    quat.proj_()
    trans = geoopt.ManifoldParameter(T.to(device)) # xyz 
    optim = geoopt.optim.RiemannianAdam([quat,trans], lr=0.1) # parameters to optimise and the learning rate

    zero_intersect = torch.zeros_like(phi_2)

    vis.line(X=np.array([0]), Y=np.array([[np.nan, np.nan, np.nan]]),
        win='loss', opts=dict(title='loss'))

    for step in range(steps):
    
        optim.zero_grad()
    
        grid = torch.nn.functional.affine_grid(build_matrix(quat, trans), phi_2.size())
        s_phi_1 = torch.nn.functional.grid_sample(phi_1, grid, mode='bilinear', padding_mode='border')
        
        # Contact
        # try without the minus to so it better tends to zero? - good contact is still minimum
        contact = heaviside(torch.min(s_phi_1, phi_2))#dirac(torch.min(s_phi_1, phi_2))
    
        # Intersection
        # Remember that the property data values now just contain the distance (essentially) to the surface.
        # Hence why we use max and min here. The heaviside penalises distances far from the contact / 
        # intersection site (which is what we want)
        k = torch.max(zero_intersect, -torch.max(s_phi_1, phi_2))
        
        # this lets us access the nonzero values
        if len(k.nonzero()) == 0:
            intersection = heaviside(k)
        else:
            print('here')
            k_tmp = torch.tensor([]).to(device)
            for i in tuple(k.nonzero()):
                k_tmp = torch.cat((k_tmp, k[tuple(i)].unsqueeze(0)))
            intersection = heaviside(k_tmp)        

        # distance loss is just sum of the square of the translation vector (from orig)
        distance = torch.sum(trans**2, dim=1)

        # clamp distance_c so we don't explode by going negative (but also we don't want 0!)     distance_c.clamp(0.000001
        # remove distance_c parmeter? contact_c so it stops going -ve all the time?
        # absolute on all parameters so functional form matters more. -> need to constrain?
        loss = torch.abs(contact_c) * contact.mean() + torch.abs(intersect_c) * intersection.mean() + torch.abs(distance_c) * distance

        loss.backward()
        optim.step()

        I = intersection.mean().item() * intersect_c
        C = contact.mean().item() * contact_c
        
        if step % 100 == 0:quat2mat
            visualise(s_phi_1, phi_2, sf=8, f=dirac)quat2mat
            vis.line(X=np.array([step]), Y=np.array([[I.item(), C.item(), distance.mean().item()]])quat2mat,
                win='loss', opts=dict(legend=['intersection', 'contact', 'distance']), update='appequat2matnd')

    return quat, trans, s_phi_1

def rotate_pdb_(pdb, points, R):
    """
    Rotate pdb coordinates by rotation matrix R
    :param pdb: biobob pdb of this host coordinates
    :param points: tensor points corresponding to the pdb coordinates
    :param R: rotation matrix (tensor object)
    """

    # First center molecule before applying rotation
    if 'center' not in pdb.properties:
        pdb.get_center()

    pdb_center = torch.tensor([pdb.properties['center']], dtype=torch.float32).to(device)
    points = points - pdb_center

    # Rotate molecule using R
    points = torch.matmul(points, R)

    # translate back to its position
    return points + pdb_center

def torch_rmsd(points0, points1):
    """
    Calculate the RMSD between two sets of points - this is essentially our loss function for the second layer of optimisation
    :param points0: Set of coordinates 0
    :param points1: Set of coordinates 1
    :returns: rmsd between the two sets of points
    """

    # Get coordinates and center them
    L = len(points0)
    COM1 = torch.sum(points0, dim=0) / float(L)
    points0 = points0 - COM1
    points0_sum = torch.sum(torch.sum(points0 * points0, dim=0), dim=0)

    COM2 = torch.sum(points1, dim=0) / float(L)
    points1 = points1 - COM2

    E0 = points0_sum + torch.sum(torch.sum(points1 * points1, dim=0), dim=0)

    V, S, W = torch.svd(torch.matmul(torch.t(points1), points0))

    reflect = float(str(float(np.linalg.det(V.cpu().detach().numpy()) * np.linalg.det(torch.t(W).cpu().detach().numpy()))))
   
    if reflect == -1.0:
        S[-1] = -S[-1]
    
    rmsdval = E0 - (2.0 * sum(S))

    return torch.sqrt(rmsdval / L)


    
#####################################################################################

# visdom
vis = visdom.Visdom(port=8345)

# Setup CUDA
device = torch.device('cuda') if torch.cuda.is_available() else torch.device('cpu')
sobel = Sobel3d().to(device)

####### Load in proteins ######
# Densities to move around
receptor = np.load("test_models/receptor_points.npy")
ligand = np.load("test_models/ligand_points.npy")

# PDBs to check RMSD
Ru = bb.Molecule()
Lu = bb.Molecule()
Ru.import_pdb('test_models/receptor.pdb')
Lu.import_pdb('test_models/ligand.pdb')

# bound case to compare with
Rb = bb.Molecule()
Lb = bb.Molecule()
Rb.import_pdb('test_models/receptor_bound.pdb')
Lb.import_pdb('test_models/ligand_bound.pdb')

# Assign atomtypes
Ru.assign_atomtype()
Lu.assign_atomtype()
Rb.assign_atomtype()
Lb.assign_atomtype()

###### Begin optimisation #####

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

#parameters to optimise with RMSDR
contact_c = torch.nn.Parameter(torch.tensor(1.).to(device))   # these might need to be global parameters 
intersect_c = torch.nn.Parameter(torch.tensor(1000.).to(device))
distance_c = torch.nn.Parameter(torch.tensor(0.1).to(device))

opts = torch.optim.Adam([contact_c, intersect_c, distance_c], lr=1.0, weight_decay = 0.1)

#### Begin looping through optimisation ####
Rb_res, Ru_res = Rb.match_residue(Ru)
Lb_res, Lu_res = Lb.match_residue(Lu)

# Get correct residues for rmsd
atoms_Rb, Rb_idxs = Rb.atomselect("*", Rb_res, "CA", get_index=True)
atoms_Ru, Ru_idxs = Ru.atomselect("*", Ru_res, "CA", get_index=True)
atoms_Lb, Lb_idxs = Lb.atomselect("*", Lb_res, "CA", get_index=True)
atoms_Lu, Lu_idxs = Lu.atomselect("*", Lu_res, "CA", get_index=True)

B = Lb.get_subset(Lb_idxs) + Rb.get_subset(Rb_idxs) # bound state molecule
Tb = torch.Tensor(B.points).to(device) # setup points tensor for RMSD
Tb.requires_grad = True

# Prep unbound structure
Lu_prep = torch.tensor(Lu.get_subset(Lu_idxs).points, dtype=torch.float32).to(device)
Lu_prep.requires_grad = True
Ru_small = Ru.get_subset(Ru_idxs)

# Number of steps for optimiser
steps = 20

# prep visdom figure
vis.line(X=np.array([0]), Y=np.array([[np.nan]]),
    win='rmsd', opts=dict(title='rmsd'))

# Prepare seed parameters for quaternion and translation
# We want to converge our coefficients - which isn't possible with random seeds at every step.
Q = torch.randn(1,4)
T = torch.randn(1,3,1)

for i in range(steps):

    opts.zero_grad()

    # Run gradient descent 
    quat, trans, s_phi_1 = fit_model(phi_1, phi_2, contact_c, intersect_c, distance_c, Q, T)
    
    # Get points and roto-translate them
    R_mat = quat2mat(quat)
    Ru_tensor = torch.tensor(Ru_small.points, dtype=torch.float32).to(device)
    Ru_tensor.requires_grad = True
    
    Ru_tensor = Ru_tensor + trans.view(-1)
    # rotate pdb in place 
    Ru_rot = rotate_pdb_(Ru_small, Ru_tensor, R_mat)

    Tu = torch.cat((Lu_prep, Ru_rot[0]), dim=0)

    rmsd = torch_rmsd(Tb, Tu) # rmsd is effectivly our loss function
    print('#############')
    print('loss is: ', rmsd.item())
    print('intersect is: ', intersect_c.item())
    print('contact is: ', contact_c.item())
    print('distance is: ', distance_c.item())

    rmsd.backward()
    opts.step()
    
    vis.line(X=np.array([i]), Y=np.array([rmsd.mean().item()]),
             win='rmsd', opts=dict(legend=['rmsd']), update='append')

    ### Activate if you want a record of the roto-translated PDBs
    if i % 5 == 0:
        protein = Tu.cpu().detach().numpy()
        S = bb.Structure(p=protein)
        S.write_pdb('test_protein_%i.pdb'%(i))
        