import mbuild as mb
import numpy as np

from . import geometry
from . import moiety
from .tools import spin

class Alcohol(mb.Compound):
    """Linear alkane with primary alcohol group.

    This generates a standard alcohol:

        CH3(CH2)_{n-1}OH

    where *n* is the number of carbons. It is subclassed to give common
    synonyms (methanol, ethanol, etc.).

    """
    def __init__(self, n):
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
        self.add(ch3,'C[$]')

        # add CH2 group(s)
        last_c = ch3
        up_ch2 = b_cc*np.array([np.sin(0.5*theta_cco),np.cos(0.5*theta_cco),0])
        down_ch2 = [1,-1,1]*up_ch2
        for i in range(1,n):
            ch2 = moiety.CH2()
            spin(ch2,0.5*np.pi,[0,1,0],ch2['C'].pos)
            if i % 2 == 1:
                ch2.translate(last_c.pos+up_ch2)
            else:
                spin(ch2,np.pi,[1,0,0],ch2['C'].pos)
                ch2.translate(last_c.pos+down_ch2)
            self.add(ch2,'C[$]')
            self.add_bond((last_c['C'],ch2['C']))
            last_c = ch2

        # add OH group
        if n % 2 == 1:
            to_oh = b_co*np.array([np.sin(0.5*theta_cco),np.cos(0.5*theta_cco),0])
            spin_h = theta_coh-(np.pi-0.5*theta_cco)
        else:
            to_oh = b_co*np.array([np.sin(0.5*theta_cco),-np.cos(0.5*theta_cco),0])
            spin_h = (np.pi-0.5*theta_cco)-theta_coh
        o = mb.Particle(name='O', element='O', pos=last_c['C'].pos+to_oh)
        ho = mb.Particle(name='H', element='H', pos=o.pos+[b_oh,0,0])
        spin(ho,spin_h,[0,0,1],o.pos)
        self.add((o,ho))
        self.add_bond((last_c['C'],o))
        self.add_bond((o,ho))

class Methanol(Alcohol):
    def __init__(self):
        super().__init__(n=1)

class Ethanol(Alcohol):
    def __init__(self):
        super().__init__(n=2)

class Propanol(Alcohol):
    def __init__(self):
        super().__init__(n=3)

class Butanol(Alcohol):
    def __init__(self):
        super().__init__(n=4)

class Implicit4SiteWater(mb.Compound):
    pass

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
