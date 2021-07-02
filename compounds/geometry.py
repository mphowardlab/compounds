import collections.abc
import numpy as np

class SymmetricMapping(collections.abc.MutableMapping):
    def __init__(self, mapping):
        self._data = {}
        self.update(mapping)

    def __contains__(self, key):
        return self._key(key) in self._data

    def __getitem__(self, key):
        return self._data[self._key(key)]

    def __setitem__(self, key, value):
        self._data[self._key(key)] = value

    def __delitem__(self, key):
        del self._data[self._key(key)]

    def __iter__(self):
        return iter(self._data)

    def __len__(self):
        return len(self._data)

    @staticmethod
    def _key(key):
        if isinstance(key,tuple) and len(key) >= 2:
            key_ = list(key)
            if key_[0] >= key_[-1]:
                key_[0],key_[-1] = key_[-1],key_[0]
            return tuple(key_)
        else:
            return key

## AMBER parameters
# https://doi.org/10.1021/ja00124a002
bond = SymmetricMapping({
    ('CT','CT'): 0.1526,
    ('CT','HC'): 0.1090,
    ('CT','H1'): 0.1090,
    ('CT','OS'): 0.1410,
    ('CT','OH'): 0.1410,
    ('HO','OH'): 0.0960
   })
angle = SymmetricMapping({
    ('CT','CT','OS'): np.radians(109.5),
    ('HC','CT','HC'): np.radians(109.5),
    ('H1','CT','OH'): np.radians(109.5),
    ('HC','CT','C'): np.radians(109.5),
    ('CT','OH','HO'): np.radians(108.5),
    })
