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

class Alkanoate(mb.Compound):
    """Linear alkanoate (carboxylate) anion.

    Generates HCOO- when n=0, or CH3(CH2)_{n-1}COO- when n >= 1,
    where n is the number of alkyl carbons.
    """
    def __init__(self, n):
        super().__init__()

        if n < 0:
            raise ValueError("n must be non-negative")

        b_cc = geometry.bond['CT', 'CT']
        b_ctc = geometry.bond['CT', 'C']
        b_co = geometry.bond['C', 'O2']
        b_h1c = geometry.bond['H1', 'C']
        theta_chain = geometry.angle['CT', 'CT', 'C']
        theta_oco = geometry.angle['O2', 'C', 'O2']
        theta_ccc = geometry.angle['CT', 'CT', 'C']

        # alkyl part (only for n >= 1)
        last_c = None
        if n >= 1:
            # start from CH3 group
            ch3 = moiety.CH3()
            spin(ch3, 0.5 * (np.pi - theta_ccc), [0, 0, 1], ch3['C'].pos)
            self.add(ch3, 'C[$]')
            last_c = ch3

            # add CH2 group(s)
            up_ch2 = b_cc * np.array([np.sin(0.5 * theta_ccc), np.cos(0.5 * theta_ccc), 0])
            down_ch2 = np.array([1, -1, 1]) * up_ch2
            for i in range(1, n):
                ch2 = moiety.CH2()
                spin(ch2, 0.5 * np.pi, [0, 1, 0], ch2['C'].pos)
                if i % 2 == 1:
                    ch2.translate(last_c['C'].pos + up_ch2)
                else:
                    spin(ch2, np.pi, [1, 0, 0], ch2['C'].pos)
                    ch2.translate(last_c['C'].pos + down_ch2)
                self.add(ch2, 'C[$]')
                self.add_bond((last_c['C'], ch2['C']))
                last_c = ch2

        # carboxylate group
        coo = moiety.COO()

        if n == 0:
            self.add(coo, 'COO')
            c_pos = coo['C'].pos
            h = mb.Particle(name='H', element='H', pos=c_pos + b_h1c * np.array([0.0, -1.0, 0.0]))
            
            self.add(h, 'H')
            self.add_bond((coo['C'], h))
        else:
            up_coo = b_ctc * np.array([np.sin(0.5 * theta_chain), np.cos(0.5 * theta_chain), 0])
            down_coo = np.array([1, -1, 1]) * up_coo
            step = up_coo if n % 2 == 1 else down_coo
            alpha = np.arctan2(-step[0], step[1])
            spin(coo, alpha, [0, 0, 1], coo['C'].pos)
            coo.translate(last_c['C'].pos + step - coo['C'].pos)

            self.add(coo, 'COO')
            self.add_bond((last_c['C'], coo['C']))


class Formate(Alkanoate):
    #HCOO-
    def __init__(self):
        super().__init__(0)


class Acetate(Alkanoate):
    #CH3COO-
    def __init__(self):
        super().__init__(1)


class Propanoate(Alkanoate):
    #CH3CH2COO-
    def __init__(self):
        super().__init__(2)


class Butanoate(Alkanoate):
    #CH3CH2CH2COO-
    def __init__(self):
        super().__init__(3)

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
