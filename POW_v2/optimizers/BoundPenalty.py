import collections, numpy as np
from numpy import inf, array, dot, exp, log, sqrt, sum

class BoundPenalty(object):

    def __init__(self, bounds=None):
        """Argument bounds can be `None` or ``bounds[0]`` and ``bounds[1]``
        are lower  and upper domain boundaries, each is either `None` or
        a scalar or a list or array of appropriate size.
        """
        ##
        # bounds attribute reminds the domain boundary values
        self.bounds = bounds

        self.gamma = 1  # a very crude assumption
        self.weights_initialized = False  # gamma becomes a vector after initialization
        self.hist = []  # delta-f history

    def has_bounds(self):
        """return True, if any variable is bounded"""
        bounds = self.bounds
        if bounds in (None, [None, None]):
            return False
        for i in xrange(bounds[0]):
            if bounds[0][i] is not None and bounds[0][i] > -np.inf:
                return True
        for i in xrange(bounds[1]):
            if bounds[1][i] is not None and bounds[1][i] < np.inf:
                return True
        return False

    def repair(self, x, bounds=None, copy=False, copy_always=False):
        """sets out-of-bounds components of ``x`` on the bounds.

        Arguments
        ---------
            `bounds`
                can be `None`, in which case the "default" bounds are used,
                or ``[lb, ub]``, where `lb` and `ub`
                represent lower and upper domain bounds respectively that
                can be `None` or a scalar or a list or array of length ``len(self)``

        code is more or less copy-paste from Solution.repair, but never tested

        """
        # TODO (old data): CPU(N,lam,iter=20,200,100): 3.3s of 8s for two bounds, 1.8s of 6.5s for one bound
        # TODO: test whether np.max([bounds[0], x], axis=0) etc is speed relevant

        if bounds is None:
            bounds = self.bounds
        if copy_always:
            x_out = array(x, copy=True)
        if bounds not in (None, [None, None], (None, None)):  # solely for effiency
            x_out = array(x, copy=True) if copy and not copy_always else x
            if bounds[0] is not None:
                if np.isscalar(bounds[0]):
                    for i in xrange(len(x)):
                        x_out[i] = max([bounds[0], x[i]])
                else:
                    for i in xrange(len(x)):
                        if bounds[0][i] is not None:
                            x_out[i] = max([bounds[0][i], x[i]])
            if bounds[1] is not None:
                if np.isscalar(bounds[1]):
                    for i in xrange(len(x)):
                        x_out[i] = min([bounds[1], x[i]])
                else:
                    for i in xrange(len(x)):
                        if bounds[1][i] is not None:
                            x_out[i] = min([bounds[1][i], x[i]])
        return x_out  # convenience return

    #____________________________________________________________
    #
    def __call__(self, x, archive, gp):
        """returns the boundary violation penalty for `x` ,where `x` is a
        single solution or a list or array of solutions.
        If `bounds` is not `None`, the values in `bounds` are used, see `__init__`"""
        if x in (None, (), []):
            return x
        if gp.bounds in (None, [None, None], (None, None)):
            return 0.0 if np.isscalar(x[0]) else [0.0] * len(x) # no penalty

        x_is_single_vector = np.isscalar(x[0])
        x = [x] if x_is_single_vector else x

        pen = []
        for xi in x:
            # CAVE: this does not work with already repaired values!!
            # CPU(N,lam,iter=20,200,100)?: 3s of 10s, array(xi): 1s (check again)
            # remark: one deep copy can be prevented by xold = xi first
            xpheno = gp.pheno(archive[xi]['geno'])
            xinbounds = gp.into_bounds(xpheno)
            fac = 1  # exp(0.1 * (log(self.scal) - np.mean(self.scal)))
            pen.append(sum(self.gamma * ((xinbounds - xpheno) / fac)**2) / len(xi))

        return pen[0] if x_is_single_vector else pen

    #____________________________________________________________
    #
    def feasible_ratio(self, solutions):
        """counts for each coordinate the number of feasible values in
        ``solutions`` and returns an array of length ``len(solutions[0])``
        with the ratios.

        `solutions` is a list or array of repaired `Solution` instances

        """
        count = np.zeros(len(solutions[0]))
        for x in solutions:
            count += x.unrepaired == x
        return count / float(len(solutions))

    #____________________________________________________________
    #
    def update(self, function_values, es, bounds=None):
        """updates the weights for computing a boundary penalty.

        Arguments
        ---------
        `function_values`
            all function values of recent population of solutions
        `es`
            `CMAEvolutionStrategy` object instance, in particular the
            method `into_bounds` of the attribute `gp` of type `GenoPheno`
            is used.
        `bounds`
            not (yet) in use other than for ``bounds == [None, None]`` nothing
            is updated.

        Reference: Hansen et al 2009, A Method for Handling Uncertainty...
        IEEE TEC, with addendum at http://www.lri.fr/~hansen/TEC2009online.pdf

        """
        if bounds is None:
            bounds = self.bounds
        if bounds is None or (bounds[0] is None and bounds[1] is None):  # no bounds ==> no penalty
            return self  # len(function_values) * [0.0]  # case without voilations

        N = es.N
        ### prepare
        # compute varis = sigma**2 * C_ii
        varis = es.sigma**2 * array(N * [es.C] if np.isscalar(es.C) else (  # scalar case
                                es.C if np.isscalar(es.C[0]) else  # diagonal matrix case
                                [es.C[i][i] for i in xrange(N)]))  # full matrix case

        # dmean = (es.mean - es.gp.into_bounds(es.mean)) / varis**0.5
        dmean = (es.mean - es.gp.geno(es.gp.into_bounds(es.gp.pheno(es.mean)))) / varis**0.5

        ### Store/update a history of delta fitness value
        fvals = sorted(function_values)
        l = 1 + len(fvals)
        val = fvals[3*l // 4] - fvals[l // 4] # exact interquartile range apart interpolation
        val = val / np.mean(varis)  # new: val is normalized with sigma of the same iteration
        # insert val in history
        if np.isfinite(val) and val > 0:
            self.hist.insert(0, val)
        elif val == inf and len(self.hist) > 1:
            self.hist.insert(0, max(self.hist))
        else:
            pass  # ignore 0 or nan values
        if len(self.hist) > 20 + (3*N) / es.popsize:
            self.hist.pop()

        ### prepare
        dfit = np.median(self.hist)  # median interquartile range
        damp = min(1, es.sp.mueff/10./N)

        ### set/update weights
        # Throw initialization error
        if len(self.hist) == 0:
            raise _Error('wrongful initialization, no feasible solution sampled. ' +
                'Reasons can be mistakenly set bounds (lower bound not smaller than upper bound) or a too large initial sigma0 or... ' +
                'See description of argument func in help(cma.fmin) or an example handling infeasible solutions in help(cma.CMAEvolutionStrategy). ')
        # initialize weights
        if (dmean.any() and (not self.weights_initialized or es.countiter == 2)):  # TODO
            self.gamma = array(N * [2*dfit])
            self.weights_initialized = True
        # update weights gamma
        if self.weights_initialized:
            edist = array(abs(dmean) - 3 * max(1, N**0.5/es.sp.mueff))
            if 1 < 3:  # this is better, around a factor of two
                # increase single weights possibly with a faster rate than they can decrease
                #     value unit of edst is std dev, 3==random walk of 9 steps
                self.gamma *= exp((edist>0) * np.tanh(edist/3) / 2.)**damp
                # decrease all weights up to the same level to avoid single extremely small weights
                #    use a constant factor for pseudo-keeping invariance
                self.gamma[self.gamma > 5 * dfit] *= exp(-1./3)**damp
                #     self.gamma[idx] *= exp(5*dfit/self.gamma[idx] - 1)**(damp/3)
            elif 1 < 3 and (edist>0).any():  # previous method
                # CAVE: min was max in TEC 2009
                self.gamma[edist>0] *= 1.1**min(1, es.sp.mueff/10./N)
                # max fails on cigtab(N=12,bounds=[0.1,None]):
                # self.gamma[edist>0] *= 1.1**max(1, es.sp.mueff/10./N) # this was a bug!?
                # self.gamma *= exp((edist>0) * np.tanh(edist))**min(1, es.sp.mueff/10./N)
            else:  # alternative version, but not better
                solutions = es.pop  # this has not been checked
                r = self.feasible_ratio(solutions)  # has to be the averaged over N iterations
                self.gamma *= exp(np.max([N*[0], 0.3 - r], axis=0))**min(1, es.sp.mueff/10/N)
        es.more_to_write += list(self.gamma) if self.weights_initialized else N * [1.0]
        ### return penalty
        # es.more_to_write = self.gamma if not np.isscalar(self.gamma) else N*[1]
        return self  # bound penalty values

#____________________________________________________________
#____________________________________________________________
#


#____________________________________________________________
#____________________________________________________________
