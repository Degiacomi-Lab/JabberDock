import collections, numpy as np

class DerivedDictBase(collections.MutableMapping):
    """for conveniently adding features to a dictionary. The actual
    dictionary is in ``self.data``. Copy-paste
    and modify setitem, getitem, and delitem, if necessary"""
    def __init__(self, *args, **kwargs):
        # collections.MutableMapping.__init__(self)
        super(DerivedDictBase, self).__init__()
        # super(SolutionDict, self).__init__()  # the same
        self.data = dict(*args, **kwargs)
    def __len__(self):
        return len(self.data)
    def __contains__(self, value):
        return value in self.data
    def __iter__(self):
        return iter(self.data)
    def __setitem__(self, key, value):
        """defines self[key] = value"""
        self.data[key] = value
    def __getitem__(self, key):
        """defines self[key]"""
        return self.data[key]
    def __delitem__(self, key):
        del self.data[key]

class SolutionDict(DerivedDictBase):

    def __init__(self, *args, **kwargs):
        DerivedDictBase.__init__(self, *args, **kwargs)
        self.data_with_same_key = {}
    def key(self, x):
        try:
            return tuple(x)
        except TypeError:
            return x
    def __setitem__(self, key, value):
        """defines self[key] = value"""
        key = self.key(key)
        if key in self.data_with_same_key:
            self.data_with_same_key[key] += [self.data[key]]
        elif key in self.data:
            self.data_with_same_key[key] = [self.data[key]]
        self.data[key] = value
    def __getitem__(self, key):
        """defines self[key]"""
        return self.data[self.key(key)]
    def __delitem__(self, key):
        """remove only most current key-entry"""
        key = self.key(key)
        if key in self.data_with_same_key:
            if len(self.data_with_same_key[key]) == 1:
                self.data[key] = self.data_with_same_key.pop(key)[0]
            else:
                self.data[key] = self.data_with_same_key[key].pop(-1)
        else:
            del self.data[key]
    def truncate(self, max_len, min_iter):
        if len(self) > max_len:
            for k in list(self.keys()):
                if self[k]['iteration'] < min_iter:
                    del self[k]  # only deletes one item with k as key, should delete all?
