import torch
import torch.nn as nn


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

class Sobel2d(nn.Module):
    def __init__(self, magnitude=True):
        super(Sobel2d, self).__init__()
        self.magnitude = magnitude
        filters = sobel_filts(2).reshape([2, 1, 3, 3])
        self.conv = torch.nn.Conv2d(
            in_channels=1,
            out_channels=2,
            kernel_size=3,
            bias=False,
            groups=1)
        for i in range(self.conv.weight.shape[0]):
            self.conv.weight.data[i] = filters[i % 2]

    def forward(self, x):
        inshape = x.shape
        outshape = (*inshape[:-2],2, *(x-2 for x in inshape[-2:]))
        out = self.conv(x.reshape((-1, 1, *inshape[-2:]))).reshape(outshape)
        if self.magnitude:
            out = out.norm(dim=-3)
        return out


class Sobel3d(nn.Module):
    def __init__(self, magnitude=True):
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