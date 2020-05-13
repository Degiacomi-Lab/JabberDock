import collections, numpy as np
from numpy import inf, array, dot, exp, log, sqrt, sum

class GenoPheno(object):

    def __init__(self, dim, scaling=None, typical_x=None, bounds=None, fixed_values=None, tf=None):

        self.N = dim
        self.bounds = bounds
        self.fixed_values = fixed_values
        if tf is not None:
            self.tf_pheno = tf[0]
            self.tf_geno = tf[1]  # TODO: should not necessarily be needed
            # r = np.random.randn(dim)
            # assert all(tf[0](tf[1](r)) - r < 1e-7)
            # r = np.random.randn(dim)
            # assert all(tf[0](tf[1](r)) - r > -1e-7)
#            self.tf_pheno = lambda x: [xi**2 for xi in x]
#            self.tf_geno = lambda x: [abs(xi)**0.5 for xi in x]

#            print("WARNING in class GenoPheno: user defined transformations have not been tested thoroughly")
        else:
            self.tf_geno = None
            self.tf_pheno = None

        if fixed_values:
            if type(fixed_values) is not dict:
                raise _Error("fixed_values must be a dictionary {index:value,...}")
            if max(fixed_values.keys()) >= dim:
                raise _Error("max(fixed_values.keys()) = " + str(max(fixed_values.keys())) +
                    " >= dim=N=" + str(dim) + " is not a feasible index")
            # convenience commenting functionality: drop negative keys
            for k in list(fixed_values.keys()):
                if k < 0:
                    fixed_values.pop(k)
        if bounds:
            if len(bounds) != 2:
                raise _Error('len(bounds) must be 2 for lower and upper bounds')
            for i in (0,1):
                if bounds[i] is not None:
                    bounds[i] = array(dim * [bounds[i]] if np.isscalar(bounds[i]) else
                                        [b for b in bounds[i]])

        def vec_is_default(vec, default_val=0):
            """return True if `vec` has the value `default_val`,
            None or [None] are also recognized as default"""
            try:
                if len(vec) == 1:
                    vec = vec[0]  # [None] becomes None and is always default
                else:
                    return False
            except TypeError:
                pass  # vec is a scalar

            if vec is None or vec == array(None) or vec == default_val:
                return True
            return False

        self.scales = array(scaling)
        if vec_is_default(self.scales, 1):
            self.scales = 1  # CAVE: 1 is not array(1)
        elif self.scales.shape is not () and len(self.scales) != self.N:
            raise _Error('len(scales) == ' + str(len(self.scales)) +
                         ' does not match dimension N == ' + str(self.N))

        self.typical_x = array(typical_x)
        if vec_is_default(self.typical_x, 0):
            self.typical_x = 0
        elif self.typical_x.shape is not () and len(self.typical_x) != self.N:
            raise _Error('len(typical_x) == ' + str(len(self.typical_x)) +
                         ' does not match dimension N == ' + str(self.N))

        if (self.scales is 1 and
                self.typical_x is 0 and
                self.bounds in (None, [None, None]) and
                self.fixed_values is None and
                self.tf_pheno is None):
            self.isidentity = True
        else:
            self.isidentity = False

    def into_bounds(self, y, bounds=None, copy_never=False, copy_always=False):
        """Argument `y` is a phenotypic vector,
        return `y` put into boundaries, as a copy iff ``y != into_bounds(y)``.

        Note: this code is duplicated in `Solution.repair` and might
        disappear in future.

        """
        bounds = bounds if bounds is not None else self.bounds
        if bounds in (None, [None, None]):
            return y if not copy_always else array(y, copy=True)
        if bounds[0] is not None:
            if len(bounds[0]) not in (1, len(y)):
                raise ValueError('len(bounds[0]) = ' + str(len(bounds[0])) +
                                 ' and len of initial solution (' + str(len(y)) + ') disagree')
            if copy_never:  # is rather slower
                for i in xrange(len(y)):
                    y[i] = max(bounds[0][i], y[i])
            else:
                y = np.max([bounds[0], y], axis=0)
        if bounds[1] is not None:
            if len(bounds[1]) not in (1, len(y)):
                raise ValueError('len(bounds[1]) = ' + str(len(bounds[1])) +
                                    ' and initial solution (' + str(len(y)) + ') disagree')
            if copy_never:
                for i in xrange(len(y)):
                    y[i] = min(bounds[1][i], y[i])
            else:
                y = np.min([bounds[1], y], axis=0)
        return y



    def pheno(self, x, bounds=None, copy=True, copy_always=False):
        """maps the genotypic input argument into the phenotypic space,
        boundaries are only applied if argument ``bounds is not None``, see
        help for class `GenoPheno`

        """
        if copy_always and not copy:
            raise ValueError('arguments copy_always=' + str(copy_always) +
                             ' and copy=' + str(copy) + ' have inconsistent values')
        if self.isidentity and bounds in (None, [None, None], (None, None)):
            return x if not copy_always else array(x, copy=copy_always)

        if self.fixed_values is None:
            y = array(x, copy=copy)  # make a copy, in case
        else:  # expand with fixed values
            y = list(x)  # is a copy
            for i in sorted(self.fixed_values.keys()):
                y.insert(i, self.fixed_values[i])
            y = array(y, copy=False)

        if self.scales is not 1:  # just for efficiency
            y *= self.scales

        if self.typical_x is not 0:
            y += self.typical_x

        if self.tf_pheno is not None:
            y = array(self.tf_pheno(y), copy=False)

        if bounds is not None:
            y = self.into_bounds(y, bounds)

        if self.fixed_values is not None:
            for i, k in list(self.fixed_values.items()):
                y[i] = k

        return y

    def geno(self, y, bounds=None, copy=True, copy_always=False, archive=None):
        """maps the phenotypic input argument into the genotypic space.
        If `bounds` are given, first `y` is projected into the feasible
        domain. In this case ``copy==False`` leads to a copy.

        by default a copy is made only to prevent to modify ``y``

        method geno is only needed if external solutions are injected
        (geno(initial_solution) is depreciated and will disappear)

        TODO: arg copy=True should become copy_never=False

        """
        if archive is not None and bounds is not None:
            try:
                return archive[y]['geno']
            except:
                pass

        x = array(y, copy=(copy and not self.isidentity) or copy_always)

        # bounds = self.bounds if bounds is None else bounds
        if bounds is not None:  # map phenotyp into bounds first
            x = self.into_bounds(x, bounds)

        if self.isidentity:
            return x

        # user-defined transformation
        if self.tf_geno is not None:
            x = array(self.tf_geno(x), copy=False)
        else:
            _Error('t1 of options transformation was not defined but is needed as being the inverse of t0')

        # affine-linear transformation: shift and scaling
        if self.typical_x is not 0:
            x -= self.typical_x
        if self.scales is not 1:  # just for efficiency
            x /= self.scales

        # kick out fixed_values
        if self.fixed_values is not None:
            # keeping the transformed values does not help much
            # therefore it is omitted
            if 1 < 3:
                keys = sorted(self.fixed_values.keys())
                x = array([x[i] for i in range(len(x)) if i not in keys], copy=False)
            else:  # TODO: is this more efficient?
                x = list(x)
                for key in sorted(list(self.fixed_values.keys()), reverse=True):
                    x.remove(key)
                x = array(x, copy=False)
        return x
#____________________________________________________________
#____________________________________________________________
# check out built-in package abc: class ABCMeta, abstractmethod, abstractproperty...
# see http://docs.python.org/whatsnew/2.6.html PEP 3119 abstract base classes
#
