

class Rotation(object):

    dicMatrices = {}  # store matrix if necessary, for each dimension
    def __init__(self):
        self.dicMatrices = {} # otherwise there might be shared bases which is probably not what we want
    def __call__(self, x, inverse=False): # function when calling an object
        """Rotates the input array `x` with a fixed rotation matrix
           (``self.dicMatrices['str(len(x))']``)
        """
        N = x.shape[0]  # can be an array or matrix, TODO: accept also a list of arrays?
        if str(N) not in self.dicMatrices: # create new N-basis for once and all
            B = np.random.randn(N, N)
            for i in xrange(N):
                for j in xrange(0, i):
                    B[i] -= np.dot(B[i], B[j]) * B[j]
                B[i] /= sum(B[i]**2)**0.5
            self.dicMatrices[str(N)] = B
        if inverse:
            return np.dot(self.dicMatrices[str(N)].T, x)  # compute rotation
        else:
            return np.dot(self.dicMatrices[str(N)], x)  # compute rotation
# Use rotate(x) to rotate x
rotate = Rotation()
