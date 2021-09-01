import mbuild as mb
import numpy as np

from . import geometry
from . import moiety
from .tools import spin

class Implicit4SiteWater(mb.Compound):
    pass

class Methanol(mb.Compound):
    def __init__(self):
        super().__init__()

        b_ch = geometry.bond['CT','HC']
        b_co = geometry.bond['CT','OS']
        b_oh = geometry.bond['OH','HO']
        theta_hco = geometry.angle['H1','CT','OH']
        theta_coh = geometry.angle['CT','OH','HO']

        c = mb.Particle(name='C', element='C', pos=[0,0,0])
        h1 = mb.Particle(name='H', element='H', pos=[b_ch,0,0])
        h1.rotate(theta_hco,around=[0,0,1])
        h2 = mb.clone(h1)
        h2.rotate(2*np.pi/3, around=[1,0,0])
        h3 = mb.clone(h1)
        h3.rotate(-2*np.pi/3, around=[1,0,0])
        o = mb.Particle(name='O', element='O', pos=[b_co,0,0])
        ho = mb.Particle(name='H', element='H', pos=[b_co-b_oh,0,0])
        spin(ho,-theta_coh,[0,0,1],o.pos)

        self.add((c,h1,h2,h3,o,ho))
        self.add_bond((c,h1))
        self.add_bond((c,h2))
        self.add_bond((c,h3))
        self.add_bond((c,o))
        self.add_bond((o,ho))

class Ethanol(mb.Compound):
    def __init__(self):
        super().__init__()

        b_cc = geometry.bond['CT','CT']
        b_co = geometry.bond['CT','OS']
        b_oh = geometry.bond['OH','HO']
        theta_hco = geometry.angle['H1','CT','OH']
        theta_coh = geometry.angle['CT','OH','HO']
        theta_cco = geometry.angle['CT','CT','OS']

        # start from CH3 group
        ch3 = moiety.CH3()
        spin(ch3,0.5*(np.pi-theta_cco),[0,0,1],ch3['C'].pos)
        self.add(ch3)

        # add CH2 group
        ch2 = moiety.CH2()
        spin(ch2,0.5*np.pi,[0,1,0],ch2['C'].pos)
        ch2.translate(b_cc*np.array([np.sin(0.5*theta_cco),np.cos(0.5*theta_cco),0]))
        self.add(ch2)
        self.add_bond((ch3['C'],ch2['C']))

        # add OH group
        o = mb.Particle(name='O', element='O', pos=ch2['C'].pos+b_co*np.array([np.sin(0.5*theta_cco),-np.cos(0.5*theta_cco),0]))
        ho = mb.Particle(name='H', element='H', pos=o.pos+[b_oh,0,0])
        spin(ho,(np.pi-0.5*theta_cco)-theta_coh,[0,0,1],o.pos)
        self.add((o,ho))
        self.add_bond((ch2['C'],o))
        self.add_bond((o,ho))

class OPCWater(Implicit4SiteWater):
    """OPC water model.

    This is a 4-site water model but LAMMPS treats the dummy site
    implicitly, so we create it as a 3-site Compound.

    https://doi.org/10.1021/jz501780a

    """
    def __init__(self):
        super().__init__()

        # OPC parameters (Table 2)
        b = 0.08724 # nm
        angle = np.radians(103.6)

        # atoms
        o = mb.Particle(name='OW', element='O', pos=[0,0,0])
        h1 = mb.Particle(name='HW', element='H', pos=[b,0,0])
        h2 = mb.Particle(
            name='HW',
            element='H',
            pos=[b*np.cos(angle),b*np.sin(angle),0]
            )
        self.add((o,h1,h2))

        # bonds
        self.add_bond((o,h1))
        self.add_bond((o,h2))

    z = 0.01594 # nm
