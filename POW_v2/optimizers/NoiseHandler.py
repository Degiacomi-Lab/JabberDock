import collections, numpy as np

class NoiseHandler(object):

    def __init__(self, N, maxevals=10, aggregate=np.median, reevals=None, epsilon=1e-7, parallel=False):

        self.lam_reeval = reevals  # 2 + popsize/20, see method indices(), originally 2 + popsize/10
        self.epsilon = epsilon
        self.parallel = parallel
        self.theta = 0.5  # originally 0.2
        self.cum = 0.3  # originally 1, 0.3 allows one disagreement of current point with resulting noiseS
        self.alphasigma = 1 + 2 / (N+10)
        self.alphaevals = 1 + 2 / (N+10)  # originally 1.5
        self.alphaevalsdown = self.alphaevals**-0.25  # originally 1/1.5
        self.evaluations = 1  # to aggregate for a single f-evaluation
        self.minevals = 1
        self.maxevals = int(np.max(maxevals))
        if hasattr(maxevals, '__contains__'):  # i.e. can deal with ``in``
            if len(maxevals) > 1:
                self.minevals = min(maxevals)
                self.evaluations = self.minevals
            if len(maxevals) > 2:
                self.evaluations = np.median(maxevals)
        self.f_aggregate = aggregate
        self.evaluations_just_done = 0  # actually conducted evals, only for documentation
        self.noiseS = 0

    def __call__(self, X, fit, func, ask=None, args=()):

        self.evaluations_just_done = 0
        if not self.maxevals or self.lam_reeval == 0:
            return 1.0
        res = self.reeval(X, fit, func, ask, args)
        if not len(res):
            return 1.0
        self.update_measure()
        return self.treat()

    def get_evaluations(self):
        """return ``self.evaluations``, the number of evalutions to get a single fitness measurement"""
        return self.evaluations

    def treat(self):
        """adapt self.evaluations depending on the current measurement value
        and return ``sigma_fac in (1.0, self.alphasigma)``

        """
        if self.noiseS > 0:
            self.evaluations = min((self.evaluations * self.alphaevals, self.maxevals))
            return self.alphasigma
        else:
            self.evaluations = max((self.evaluations * self.alphaevalsdown, self.minevals))
            return 1.0

    def reeval(self, X, fit, func, ask, args=()):
        """store two fitness lists, `fit` and ``fitre`` reevaluating some
        solutions in `X`.
        ``self.evaluations`` evaluations are done for each reevaluated
        fitness value.
        See `__call__()`, where `reeval()` is called.

        """
        self.fit = list(fit)
        self.fitre = list(fit)
        self.idx = self.indices(fit)
        if not len(self.idx):
            return self.idx
        evals = int(self.evaluations) if self.f_aggregate else 1
        fagg = np.median if self.f_aggregate is None else self.f_aggregate
        for i in self.idx:
            if self.epsilon:
                if self.parallel:
                    self.fitre[i] = fagg(func(ask(evals, X[i], self.epsilon), *args))
                else:
                    self.fitre[i] = fagg([func(ask(1, X[i], self.epsilon)[0], *args)
                                            for _k in xrange(evals)])
            else:
                self.fitre[i] = fagg([func(X[i], *args) for _k in xrange(evals)])
        self.evaluations_just_done = evals * len(self.idx)
        return self.fit, self.fitre, self.idx

    def update_measure(self):
        """updated noise level measure using two fitness lists ``self.fit`` and
        ``self.fitre``, return ``self.noiseS, all_individual_measures``.

        Assumes that `self.idx` contains the indices where the fitness
        lists differ

        """
        lam = len(self.fit)
        idx = np.argsort(self.fit + self.fitre)
        ranks = np.argsort(idx).reshape((2, lam))
        rankDelta = ranks[0] - ranks[1] - np.sign(ranks[0] - ranks[1])

        # compute rank change limits using both ranks[0] and ranks[1]
        r = np.arange(1, 2 * lam)  # 2 * lam - 2 elements
        limits = [0.5 * (Mh.prctile(np.abs(r - (ranks[0,i] + 1 - (ranks[0,i] > ranks[1,i]))),
                                      self.theta*50) +
                         Mh.prctile(np.abs(r - (ranks[1,i] + 1 - (ranks[1,i] > ranks[0,i]))),
                                      self.theta*50))
                    for i in self.idx]
        # compute measurement
        #                               max: 1 rankchange in 2*lambda is always fine
        s = np.abs(rankDelta[self.idx]) - Mh.amax(limits, 1)  # lives roughly in 0..2*lambda
        self.noiseS += self.cum * (np.mean(s) - self.noiseS)
        return self.noiseS, s

    def indices(self, fit):
        """return the set of indices to be reevaluted for noise measurement,
        taking the ``lam_reeval`` best from the first ``2 * lam_reeval + 2``
        values.

        Given the first values are the earliest, this is a useful policy also
        with a time changing objective.

        """
        lam = self.lam_reeval if self.lam_reeval else 2 + len(fit) / 20
        reev = int(lam) + ((lam % 1) > np.random.rand())
        return np.argsort(array(fit, copy=False)[:2 * (reev + 1)])[:reev]
