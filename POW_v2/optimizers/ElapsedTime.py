import time

class ElapsedTime(object):
    """32-bit C overflows after int(2**32/1e6) == 4294s about 72 min"""
    def __init__(self):
        self.tic0 = time.clock()
        self.tic = self.tic0
        self.lasttoc = time.clock()
        self.lastdiff = time.clock() - self.lasttoc
        self.time_to_add = 0
        self.messages = 0

    def __call__(self):
        toc = time.clock()
        if toc - self.tic >= self.lasttoc - self.tic:
            self.lastdiff = toc - self.lasttoc
            self.lasttoc = toc
        else:  # overflow, reset self.tic
            if self.messages < 3:
                self.messages += 1
                print('  in cma.ElapsedTime: time measure overflow, last difference estimated from',
                        self.tic0, self.tic, self.lasttoc, toc, toc - self.lasttoc, self.lastdiff)

            self.time_to_add += self.lastdiff + self.lasttoc - self.tic
            self.tic = toc  # reset
            self.lasttoc = toc
        self.elapsedtime = toc - self.tic + self.time_to_add
        return self.elapsedtime

#____________________________________________________________
#____________________________________________________________
#
