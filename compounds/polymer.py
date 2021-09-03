import mbuild as mb
import numpy as np

from . import geometry
from . import moiety
from .tools import spin

class CrosslinkedPEGDA(mb.Compound):
    def __init__(self, n):
        super().__init__()

        b_cc = geometry.bond['CT','CT']
        b_co = geometry.bond['CT','OS']

        # build the crosslink junction as a dodecagon ring
        #
        # law of cosines determines radius of dodecagon with angle 30*
        # and a C-C bond between vertices
        theta = np.radians(30)
        l = np.sqrt(0.5*b_cc**2/(1-np.cos(theta)))
        first = None
        last = None
        for i in range(12):
            if i % 2 == 0:
                if (i//2) % 2 == 0:
                    current = moiety.HCCO()
                else:
                    current = moiety.HCCOO()
                    spin(current,np.pi,[0,1,0],current['C'].pos)
                label = 'site[$]'
            else:
                current = moiety.CH2()
                label = None
            # spin and twist based on angle in circle
            spin(current,-np.pi/2,[0,1,0],current['C'].pos)
            spin(current,i*theta-np.pi/2,[0,0,1],current['C'].pos)
            # translate to vertex
            current.translate([l*np.cos(i*theta),l*np.sin(i*theta),0])
            self.add(current, label=label)
            if last is not None:
                self.add_bond((last['C'],current['C']))
            if first is None:
                first = current
            last = current
        self.add_bond((last['C'],first['C']))

        # attach PEO chains to the crosslink junction in each of cartesian directions
        L = [None,None,None]
        for i in range(3):
            site_down = self['site'][(2*i+3)%6]
            site_up = self['site'][2*i]

            # create initial PEO aligned along +z
            peo = PolyethyleneOxide(n)
            first = peo['monomer'][0]
            last = peo['monomer'][-1]

            # setup target bond displacements
            if i == 0:
                disp = np.array([b_co,0,0])
                spin(site_up['O='],np.pi,site_up['C='].pos-site_up['C'].pos,site_up['C='].pos)
            elif i == 1:
                disp = np.array([0,b_co,0])
                spin(site_up['O='],np.pi,site_up['C='].pos-site_up['C'].pos,site_up['C='].pos)
            elif i == 2:
                disp = np.array([0,0,b_co])
                #spin(site_up['O='],np.pi,site_up['C='].pos-site_up['C'].pos,site_up['C='].pos)
            port_down = site_down['O'].pos - disp
            port_up = site_up['C='].pos + disp

            # spin the PEO that will bridge each dimension (0 = x, 1 = y, 2 = z)
            rpp = port_down-port_up
            if i == 0:
                # rotate over junction
                ree = peo.Ree()
                twist = -np.arcsin((rpp[2]-ree[2])/ree[0])
                spin(peo,twist,[0,1,0],first['O'].pos)
                # rotate to port
                twist = np.arcsin((rpp[1]-ree[1])/ree[0])
                spin(peo,twist,[0,0,1],first['O'].pos)
            elif i == 1:
                # orient to +y
                spin(peo,np.pi/2.,[0,1,0],first['O'].pos)
                spin(peo,np.pi/2,[1,0,0],first['O'].pos)
                # rotate over the junction
                ree = peo.Ree()
                twist = np.arcsin((rpp[2]-ree[2])/ree[1])
                spin(peo,twist,[1,0,0],first['O'].pos)
                # rotate to port
                twist = -np.arcsin((rpp[0]-ree[0])/ree[1])
                spin(peo,twist,[0,0,1],first['O'].pos)
            elif i == 2:
                # orient to +z
                spin(peo,-np.pi/2,[0,1,0],first['O'].pos)
                # rotate over the junction
                ree = peo.Ree()
                twist = np.arcsin((rpp[0]-ree[0])/ree[2])
                spin(peo,twist,[0,1,0],first['O'].pos)
                # rotate to port
                twist = -np.arcsin((rpp[1]-ree[1])/ree[2])
                spin(peo,twist,[1,0,0],first['O'].pos)
            else:
                raise NotImplementedError('Can coordinate max of 3 chains')

            # get periodic size of box
            ree = peo.Ree()
            L[i] = np.abs(port_up[i]-port_down[i]+ree[i])

            peo.translate(port_up)
            self.add(peo,'peo[$]')
            self.add_bond((first['O'],site_up['C=']))
            peo.remove(peo['monomer'][0]['down'])

            # add ports to this compound, removing from constituent PEOs
            p = mb.Port(anchor=site_down['O'])
            p.translate_to(port_down)
            self.add(p,'down[$]')
            p = mb.Port(anchor=peo['monomer'][-1]['up'].anchor)
            p.translate_to(peo['monomer'][-1]['up'].pos)
            self.add(p,'up[$]')
            peo.remove(peo['monomer'][-1]['up'])
        self.L = np.array(L)

class IdealPEGDANetwork(mb.Compound):
    def __init__(self, n, num_junctions):
        super().__init__()

        if isinstance(num_junctions, int):
            num_junctions = (num_junctions,num_junctions,num_junctions)
        else:
            num_junctions = tuple(num_junctions)
            if len(num_junctions) != 3:
                raise TypeError('Number of crosslinks must be int or 3-tuple')

        # assemble the junctions, then bond them
        xs = []
        for i in range(num_junctions[0]):
            xs.append([])
            for j in range(num_junctions[1]):
                xs[i].append([])
                for k in range(num_junctions[2]):
                    x = CrosslinkedPEGDA(n)
                    x.translate(x.L*np.array([i,j,k]))
                    self.add(x)
                    xs[i][j].append(x)

        for i in range(num_junctions[0]):
            for j in range(num_junctions[1]):
                for k in range(num_junctions[2]):
                    up = xs[i][j][k]
                    # +x
                    down = xs[(i+1)%num_junctions[0]][j][k]
                    self.add_bond((up['up'][0].anchor,down['down'][0].anchor))
                    up.remove(up['up'][0])
                    down.remove(down['down'][0])
                    # +y
                    down = xs[i][(j+1)%num_junctions[1]][k]
                    self.add_bond((up['up'][1].anchor,down['down'][1].anchor))
                    up.remove(up['up'][1])
                    down.remove(down['down'][1])
                    # +z
                    down = xs[i][j][(k+1)%num_junctions[2]]
                    self.add_bond((up['up'][2].anchor,down['down'][2].anchor))
                    up.remove(up['up'][2])
                    down.remove(down['down'][2])

        # try to shift polymer to keep it mostly in the bounding box
        self.translate(-np.amin(self.xyz,axis=0))

        self.box = mb.Box(xs[0][0][0].L*num_junctions)
        self.periodicity = [True,True,True]

class PolyethyleneOxide(mb.Compound):
    class Monomer(mb.Compound):
        def __init__(self):
            super().__init__()

            # oxygen at origin
            o = mb.Particle(name='O', element='O', pos=[0,0,0])
            self.add(o,'O')

            # orient backbone down x axis
            b_cc = geometry.bond['CT','CT']
            b_co = geometry.bond['CT','OS']
            theta_cco = geometry.angle['CT','CT','OS']
            oc1 = b_co*np.array([np.sin(0.5*theta_cco), np.cos(0.5*theta_cco), 0])
            c1c2 = b_cc*np.array([np.sin(0.5*theta_cco), -np.cos(0.5*theta_cco), 0])

            # first CH2, translated
            c1 = moiety.CH2()
            spin(c1,np.pi/2,[0,1,0],c1['C'].pos)
            c1.translate(oc1)
            self.add(c1,'C1')

            # second CH2, spun around and translated
            c2 = moiety.CH2()
            spin(c2,np.pi/2,[0,1,0],c2['C'].pos)
            c2.translate(oc1+c1c2)
            spin(c2,np.pi,[1,0,0],c2['C'].pos)
            self.add(c2,'C2')

            # bond
            self.add_bond((o,c1['C']))
            self.add_bond((c1['C'],c2['C']))

            # ports to next monomer
            self.add(mb.Port(anchor=o),'down')
            self['down'].translate(0.5*oc1*[-1,1,0])
            self.add(mb.Port(anchor=c2['C']), 'up')
            self['up'].translate(0.5*oc1)

    def __init__(self, n):
        super().__init__()

        # polymerize n monomers
        last = None
        for i in range(n):
            # manually make monomers overlap (get wrong dihedral with force_overlap)
            current = self.Monomer()
            if i % 2 == 1:
                spin(current,np.pi,[1,0,0],current['O'].pos)
            self.add(current,label='monomer[$]')
            if last is not None:
                c2o = 2*(last['up'].pos-last['C2']['C'].pos)
                current.translate(last['C2']['C'].pos+c2o)
                self.add_bond((last['C2']['C'],current['O']))
                self.remove(last['up'])
                self.remove(current['down'])

            last = current

    def Ree(self):
        return self['monomer'][-1]['C2']['C'].pos-self['monomer'][0]['O'].pos
