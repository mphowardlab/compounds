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
    ('HO','OH'): 0.0960,
    ('C','CT'): 0.1522,
    ('C','O2'): 0.1250,
    ('H1','C'): 0.1090
   })
angle = SymmetricMapping({
    ('CT','CT','OS'): np.radians(109.5),
    ('HC','CT','HC'): np.radians(109.5),
    ('H1','CT','OH'): np.radians(109.5),
    ('HC','CT','C'): np.radians(109.5),
    ('CT','OH','HO'): np.radians(108.5),
    ('CT','C','O2'): np.radians(117.00),
    ('C','CT','HC'): np.radians(109.50),
    ('CT','CT','C'): np.radians(109.5),
    ('O2','C','O2'): np.radians(126.00),
    ('H1','C','O2'): np.radians(109.5)
    })
