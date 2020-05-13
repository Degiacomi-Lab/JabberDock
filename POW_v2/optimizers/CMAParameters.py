from numpy import inf, array, dot, exp, log, sqrt, sum
import collections, numpy as np

class BlancClass(object):
    """blanc container class for having a collection of attributes"""

class CMAParameters(object):

    def __init__(self, N, opts, ccovfac=1, verbose=True):
        """Compute strategy parameters, mainly depending on
        dimension and population size, by calling `set`

        """
        self.N = N
        if ccovfac == 1:
            ccovfac = opts['CMA_on']  # that's a hack
        self.set(opts, ccovfac=ccovfac, verbose=verbose)

    def set(self, opts, popsize=None, ccovfac=1, verbose=True):
        """Compute strategy parameters as a function
        of dimension and population size """

        alpha_cc = 1.0  # cc-correction for mueff, was zero before

        def cone(df, mu, N, alphacov=2.0):
            """rank one update learning rate, ``df`` is disregarded and obsolete, reduce alphacov on noisy problems, say to 0.5"""
            return alphacov / ((N + 1.3)**2 + mu)

        def cmu(df, mu, alphamu=0.0, alphacov=2.0):
            """rank mu learning rate, disregarding the constrant cmu <= 1 - cone"""
            c = alphacov * (alphamu + mu - 2 + 1/mu) / ((N + 2)**2 + alphacov * mu / 2)
            # c = alphacov * (alphamu + mu - 2 + 1/mu) / (2 * (N + 2)**1.5 + alphacov * mu / 2)
            # print 'cmu =', c
            return c

        def conedf(df, mu, N):
            """used for computing separable learning rate"""
            return 1. / (df + 2.*sqrt(df) + float(mu)/N)

        def cmudf(df, mu, alphamu):
            """used for computing separable learning rate"""
            return (alphamu + mu - 2. + 1./mu) / (df + 4.*sqrt(df) + mu/2.)

        # ------------------ DEFINING POPSIZE
        sp = self
        N = sp.N
        if popsize:
            opts.evalall({'N':N, 'popsize':popsize})
        else:
            popsize = opts.evalall({'N':N})['popsize']  # the default popsize is computed in Options()
        sp.popsize = popsize
        if opts['CMA_mirrors'] < 0.5:
            sp.lam_mirr = int(0.5 + opts['CMA_mirrors'] * popsize)
        elif opts['CMA_mirrors'] > 1:
            sp.lam_mirr = int(0.5 + opts['CMA_mirrors'])
        else:
            sp.lam_mirr = int(0.5 + 0.16 * min((popsize, 2 * N + 2)) + 0.29)  # 0.158650... * popsize is optimal
            # lam = arange(2,22)
            # mirr = 0.16 + 0.29/lam
            # print(lam); print([int(0.5 + l) for l in mirr*lam])
            # [ 2  3  4  5  6  7  8  9 10 11 12 13 14 15 16 17 18 19 20 21]
            # [1, 1, 1, 1, 1, 1, 2, 2, 2, 2, 2, 2, 3, 3, 3, 3, 3, 3, 3, 4]

        # ------------------ DEFINING MU
        sp.mu_f = sp.popsize / 2.0  # float value of mu
        if opts['CMA_mu'] is not None:
            sp.mu_f = opts['CMA_mu']

        sp.mu = int(sp.mu_f + 0.499999) # round down for x.5
        # in principle we have mu_opt = popsize/2 + lam_mirr/2,
        # which means in particular weights should only be negative for q > 0.5+mirr_frac/2
        if sp.mu > sp.popsize - 2 * sp.lam_mirr + 1:
            print("WARNING: pairwise selection is not implemented, therefore " +
                  " mu = %d > %d = %d - 2*%d + 1 = popsize - 2*mirr + 1 can produce a bias" % (
                    sp.mu, sp.popsize - 2 * sp.lam_mirr + 1, sp.popsize, sp.lam_mirr))
        if sp.lam_mirr > sp.popsize // 2:
            raise _Error("fraction of mirrors in the population as read from option CMA_mirrors cannot be larger 0.5, " +
                         "theoretically optimal is 0.159")

        # ------------------ DEFINING WEIGHTS
#        sp.weights = log(max([sp.mu, sp.popsize / 2.0]) + 0.5) - log(1 + np.arange(sp.mu))
        sp.weights = np.ones(sp.mu)
#        print sp.weights


#        if 11 < 3:  # equal recombination weights
#            sp.mu = sp.popsize // 4
#            sp.weights = np.ones(sp.mu)
#            print(sp.weights[:10])
        sp.weights /= sum(sp.weights)
        sp.mueff = 1 / sum(sp.weights**2)
        sp.cs = (sp.mueff + 2) / (N + sp.mueff + 3)
        # TODO: clean up (here the cumulation constant is shorter if sigma_vec is used)
        sp.dampsvec = opts['CMA_dampsvec_fac'] * (N + 2) if opts['CMA_dampsvec_fac'] else np.Inf
        sp.dampsvec_fading = opts['CMA_dampsvec_fade']
        if np.isfinite(sp.dampsvec):
            sp.cs = ((sp.mueff + 2) / (N + sp.mueff + 3))**0.5
        # sp.cs = (sp.mueff + 2) / (N + 1.5*sp.mueff + 1)
        sp.cc = (4 + alpha_cc * sp.mueff / N) / (N + 4 + alpha_cc * 2 * sp.mueff / N)
        sp.cc_sep = (1 + 1/N + alpha_cc * sp.mueff / N) / (N**0.5 + 1/N + alpha_cc * 2 * sp.mueff / N) # \not\gg\cc
        sp.rankmualpha = opts['CMA_rankmualpha']
        # sp.rankmualpha = _evalOption(opts['CMA_rankmualpha'], 0.3)
        sp.c1 = ccovfac * min(1, sp.popsize/6) * cone((N**2 + N) / 2, sp.mueff, N) # 2. / ((N+1.3)**2 + sp.mucov)
        sp.c1_sep = ccovfac * conedf(N, sp.mueff, N)
        if 11 < 3:
            sp.c1 = 0.
            print('c1 is zero')
        if opts['CMA_rankmu'] != 0:  # also empty
            sp.cmu = min(1 - sp.c1, ccovfac * cmu((N**2+N)/2, sp.mueff, sp.rankmualpha))
            sp.cmu_sep = min(1 - sp.c1_sep, ccovfac * cmudf(N, sp.mueff, sp.rankmualpha))
        else:
            sp.cmu = sp.cmu_sep = 0

        sp.neg = BlancClass()
        if opts['CMA_active']:
            # in principle we have mu_opt = popsize/2 + lam_mirr/2,
            # which means in particular weights should only be negative for q > 0.5+mirr_frac/2
            sp.neg.mu_f = popsize - (popsize + sp.lam_mirr) / 2  if popsize > 2 else 1
            sp.neg.weights = log(sp.mu_f + 0.5) - log(1 + np.arange(sp.popsize - int(sp.neg.mu_f), sp.popsize))
            sp.neg.mu = len(sp.neg.weights)  # maybe never useful?
            sp.neg.weights /= sum(sp.neg.weights)
            sp.neg.mueff = 1 / sum(sp.neg.weights**2)
            sp.neg.cmuexp = opts['CMA_activefac'] * 0.25 * sp.neg.mueff / ((N+2)**1.5 + 2 * sp.neg.mueff)
            assert sp.neg.mu >= sp.lam_mirr  # not really necessary
            # sp.neg.minresidualvariance = 0.66  # not it use, keep at least 0.66 in all directions, small popsize is most critical
        else:
            sp.neg.cmuexp = 0

        sp.CMA_on = sp.c1 + sp.cmu > 0
        # print(sp.c1_sep / sp.cc_sep)

        if not opts['CMA_on'] and opts['CMA_on'] not in (None,[],(),''):
            sp.CMA_on = False
            # sp.c1 = sp.cmu = sp.c1_sep = sp.cmu_sep = 0

        sp.damps = opts['CMA_dampfac'] * (0.5 +
                                          0.5 * min([1, (sp.lam_mirr/(0.159*sp.popsize) - 1)**2])**1 +
                                          2 * max([0, ((sp.mueff-1) / (N+1))**0.5 - 1]) + sp.cs
                                          )
        if 11 < 3:
            # this is worse than damps = 1 + sp.cs for the (1,10000)-ES on 40D parabolic ridge
            sp.damps = 0.3 + 2 * max([sp.mueff/sp.popsize, ((sp.mueff-1)/(N+1))**0.5 - 1]) + sp.cs
        if 11 < 3:
            # this does not work for lambda = 4*N^2 on the parabolic ridge
            sp.damps = opts['CMA_dampfac'] * (2 - 0*sp.lam_mirr/sp.popsize) * sp.mueff/sp.popsize + 0.3 + sp.cs  # nicer future setting
            print('damps =', sp.damps)
        if 11 < 3:
            sp.damps = 10 * sp.damps  # 1e99 # (1 + 2*max(0,sqrt((sp.mueff-1)/(N+1))-1)) + sp.cs;
            # sp.damps = 20 # 1. + 20 * sp.cs**-1  # 1e99 # (1 + 2*max(0,sqrt((sp.mueff-1)/(N+1))-1)) + sp.cs;
            print('damps is %f' % (sp.damps))

        sp.cmean = float(opts['CMA_cmean'])
        # sp.kappa = 1  # 4-D, lam=16, rank1, kappa < 4 does not influence convergence rate
                        # in larger dim it does, 15-D with defaults, kappa=8 factor 2
        if sp.cmean != 1:
            print('  cmean = %f' % (sp.cmean))

        if verbose:
            if not sp.CMA_on:
                print('covariance matrix adaptation turned off')
            if opts['CMA_mu'] != None:
                print('mu = %f' % (sp.mu_f))

        # return self  # the constructor returns itself

    def disp(self):
        pprint(self.__dict__)
