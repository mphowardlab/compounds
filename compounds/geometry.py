import numpy as np

from .collections import SymmetricMapping

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
