import numpy as np

def spin(mol, angle, axis, origin=None):
    if origin is None:
        origin = mol.center
    origin = np.asarray(origin)

    mol.translate(-origin)
    mol.rotate(angle,axis)
    mol.translate(origin)
    return mol

