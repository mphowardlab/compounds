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

class TetraPEGcenter(mb.Compound):
    """Center of the tetraPEG polymer.

    This is a carbon with four CO groups attached, each of which can
    be used to attach a PEG chain.

    """

    def __init__(self):
        super().__init__()

        # carbon at origin
        c = mb.Particle(name='C', element='C', pos=[0,0,0])
        self.add(c,'C')

        # place carbon arms and hydrogens
        c1 = mb.Particle(name='C1', element='C', pos=[0.107911, 0.095853, 0.057694])
        c1h1 = mb.Particle(name="C1H1", element="H", pos=[0.06053, 0.17198, 0.12032])
        c1h2 = mb.Particle(name="C1H2", element="H", pos=[0.15933, 0.14676, -0.02429])

        c2 = mb.Particle(name='C2', element='C', pos=[0.07170, -0.10494, -0.09039])
        c2h1 = mb.Particle(name="C2H1", element="H", pos=[0.11844, -0.05440, -0.17532])
        c2h2 = mb.Particle(name="C2H2", element="H", pos=[0.15039, -0.15501, -0.03337])

        c3 = mb.Particle(name='C3', element='C', pos=[-0.07549, -0.07157, 0.11575])
        c3h1 = mb.Particle(name="C3H1", element="H", pos=[-0.15319, -0.13690, 0.07510])
        c3h2 = mb.Particle(name="C3H2", element="H", pos=[-0.00565, -0.13355, 0.17261])

        c4 = mb.Particle(name='C4', element='C', pos=[-0.09778, 0.08787, -0.08377])
        c4h1 = mb.Particle(name="C4H1", element="H", pos=[-0.13560, 0.16946, -0.02154])
        c4h2 = mb.Particle(name="C4H2", element="H", pos=[-0.04421, 0.13187, -0.16830])

        o1 = mb.Particle(name='O_1', element='O', pos=[0.20502, 0.02589, 0.13608])
        o2 = mb.Particle(name='O_2', element='O', pos=[-0.01703, -0.20568, -0.13953])
        o3 = mb.Particle(name='O_3', element='O', pos=[-0.13654, 0.02285, 0.20438])
        o4 = mb.Particle(name='O_4', element='O', pos=[-0.20989, 0.01454, -0.13365])

        # add atoms to the compound
        self.add(c1, 'C1')
        self.add(c1h1, 'C1H1')
        self.add(c1h2, 'C1H2')
        self.add(c2, 'C2')
        self.add(c2h1, 'C2H1')
        self.add(c2h2, 'C2H2')
        self.add(c3, 'C3')
        self.add(c3h1, 'C3H1')
        self.add(c3h2, 'C3H2')
        self.add(c4, 'C4')
        self.add(c4h1, 'C4H1')
        self.add(c4h2, 'C4H2')
        self.add(o1, 'O_1')
        self.add(o2, 'O_2')
        self.add(o3, 'O_3')
        self.add(o4, 'O_4')

        # bonds to central carbon
        self.add_bond((c, c1))
        self.add_bond((c, c2))
        self.add_bond((c, c3))
        self.add_bond((c, c4))

        # bond c's to hydrogens
        self.add_bond((c1, c1h1))
        self.add_bond((c1, c1h2))
        self.add_bond((c2, c2h1))
        self.add_bond((c2, c2h2))
        self.add_bond((c3, c3h1))
        self.add_bond((c3, c3h2))
        self.add_bond((c4, c4h1))
        self.add_bond((c4, c4h2))

        # add CO groups
        self.add_bond((c1, o1))
        self.add_bond((c2, o2))
        self.add_bond((c3, o3))
        self.add_bond((c4, o4))

class Implicit4SiteWater(mb.Compound):
    pass

class OPCWater(Implicit4SiteWater):
    """OPC water model.

    This is a 4-site water model but LAMMPS treats the dummy site
    implicitly, so we create it as a 3-site Compound.

    https://doi.org/10.1021/jz501780a

    """
    def __init__(self):
        super().__init__(name="SOL")

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
