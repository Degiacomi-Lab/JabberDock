import numpy as np
import numba as nb
import math
from numba import cuda

bits = 64

@cuda.jit("void(float{}[:, :], float{}[:, :], float{}[:, :])".format(bits, bits, bits))
def distance_matrix(mat1, mat2, out):
    '''
    Optional distance matrix fucntion that takes advantage of a computers GPU

    :param mat1: array of points 1
    :param mat2: array of points 2
    :param out: array of distances output, of size mat1[0] x mat2[0]
    '''
    m = mat1.shape[0]
    ndim = mat1.shape[1]
    n = mat2.shape[0]
    i, j = cuda.grid(2)
   
    if i < m and j < n:
        d = 0
        for k in range(ndim):
            tmp = mat1[i, k] - mat2[j, k]
            d += tmp * tmp
        out[i, j] = math.sqrt(d)

def pairwise_distance(M1, M2):
    "Host Code"
    rows = M1.shape[0]
    cols = M2.shape[0]
    block_dim = (16, 16)
    grid_dim = (int(rows/block_dim[0] + 1), int(cols/block_dim[1] + 1))

    stream = cuda.stream()
    M1_mem = cuda.to_device(np.asarray(M1, dtype=np.float64), stream=stream)
    M2_mem = cuda.to_device(np.asarray(M2, dtype=np.float64), stream=stream)
    out2 = cuda.device_array((rows, cols))
    distance_matrix[grid_dim, block_dim](M1_mem, M2_mem, out2)
    out = out2.copy_to_host(stream=stream)

    return out

