import mbuild as mb
import numpy as np

from . import geometry

class CH2(mb.Compound):
    def __init__(self):
        super().__init__()

        b = geometry.bond['CT','HC']
        theta = geometry.angle['HC','CT','HC']

        c = mb.Particle(name='C', element='C', pos=[0,0,0])
        h1 = mb.Particle(name='H', element='H', pos=[0,b,0])
        h1.rotate(0.5*theta,around=[0,0,1])
        h2 = mb.Particle(name='H', element='H', pos=[0,b,0])
        h2.rotate(-0.5*theta,around=[0,0,1])

        self.add(c,'C')
        self.add((h1,h2))
        self.add_bond((c,h1))
        self.add_bond((c,h2))

class CH3(mb.Compound):
    def __init__(self):
        super().__init__()

        b = geometry.bond['CT','HC']
        theta = geometry.angle['HC','CT','HC']

        c = mb.Particle(name='C', element='C', pos=[0,0,0])
        h1 = mb.Particle(name='H', element='H', pos=[b,0,0])
        h1.rotate(theta,around=[0,0,1])
        h2 = mb.clone(h1)
        h2.rotate(2*np.pi/3, around=[1,0,0])
        h3 = mb.clone(h1)
        h3.rotate(-2*np.pi/3, around=[1,0,0])

        self.add(c,'C')
        self.add((h1,h2,h3))
        self.add_bond((c,h1))
        self.add_bond((c,h2))
        self.add_bond((c,h3))

class HCCO(mb.Compound):
    def __init__(self):
        super().__init__()

        b_cc = geometry.bond['CT','CT']
        b_co = geometry.bond['CT','OS']
        b_ch = geometry.bond['CT','HC']
        theta_hcc = geometry.angle['HC','CT','C']

        c1 = mb.Particle(name='C', element='C', pos=[0,0,0])
        h = mb.Particle(name='H', element='H', pos=[0,b_ch,0])
        h.rotate(0.5*theta_hcc,around=[0,0,1])
        c2 = mb.Particle(name='C', element='C', pos=[0,b_cc,0])
        c2.rotate(-0.5*theta_hcc,around=[0,0,1])
        o = mb.Particle(name='O', element='O', pos=c2.pos+[0,b_co,0])

        self.add(c1,'C')
        self.add(c2,'C=')
        self.add(h,'H')
        self.add(o,'O=')
        self.add_bond((c1,h))
        self.add_bond((c1,c2))
        self.add_bond((c2,o))

class HCCOO(HCCO):
    def __init__(self):
        super().__init__()

        b_co = geometry.bond['CT','OS']
        theta_cco = geometry.angle['CT','CT','OS']

        co = b_co*np.array([np.sin(0.5*theta_cco), -np.cos(0.5*theta_cco), 0])
        o2 = mb.Particle(name='O', element='O', pos=self['C='].pos+co)
        self.add(o2,'O')
        self.add_bond((self['C='],o2))
