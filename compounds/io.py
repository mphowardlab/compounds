from . import molecule
from .collections import SymmetricMapping

def save_lammps(compound, filename, cutoff, **kwargs):
    if '.lammps' not in filename and '.lmp' not in filename:
        raise ValueError('File must be LAMMPS data format .lammps or .lmp')

    if 'atom_style' in kwargs or 'unit_style' in kwargs:
        raise KeyError('atom_style and unit_style are fixed as full and real')

    compound.save(filename, atom_style='full', unit_style='real', **kwargs)

    # extract types from written file
    sections = (('Masses','atom',1),
                ('Pair Coeffs','pair',0),
                ('Bond Coeffs','bond',2),
                ('Angle Coeffs','angle',3),
                ('Dihedral Coeffs','dihedral',4))
    type_map = {}
    potential_map = {'pair': None, 'bond': None, 'angle': None, 'dihedral': None}
    h_types = []
    with open(filename) as f:
        line = f.readline()
        section = 0
        while len(line) > 0 and section < len(sections):
            if sections[section][0] in line:
                section_key = sections[section][1]
                section_count = sections[section][2]

                # extract potential type from coeffs
                header = line.split()
                if section_key in potential_map:
                    potential_map[section_key] = header[-1]

                # extract types for each potential
                if section_count > 0:
                    type_map[section_key] = SymmetricMapping()
                f.readline()
                line = f.readline()
                while len(line) > 0:
                    if line == '\n':
                        break
                    result = line.strip().split()

                    if section_count > 1:
                        result_key = tuple(result[-section_count:])
                        type_map[section_key][result_key] = result[0]
                    elif section_count == 1:
                        result_key = result[-1]
                        type_map[section_key][result_key] = result[0]

                    if section_key == 'atom':
                        m = float(result[1])
                        if np.isclose(m,1.008,atol=1e-4,rtol=0):
                            h_types.append(result[-1])

                    line = f.readline()

                section += 1
            else:
                line = f.readline()

    # write forcefield stub
    with open(filename + '.ff','w') as f:
        f.write('# generated using mbuild {} / foyer {}\n\n'.format(mb.__version__,foyer.__version__))

        f.write('atom_style full\n')
        f.write('units real\n')

        f.write('\n# types\n')
        for section in type_map:
            for key,value in type_map[section].items():
                if not isinstance(key,str):
                    key = '_'.join(key)
                key = key.replace(' ','_')
                f.write('variable {}_{} equal {}\n'.format(section,key,value))

        f.write('\n# potentials\n')
        for p,sec in potential_map.items():
            if p == 'pair':
                # check for tip4p pairs
                tip4p_z = None
                for c in compounds_:
                    if isinstance(c,molecule.Implicit4SiteModel):
                        if tip4p_z is None:
                            tip4p_z = 10*c.z # A
                        else:
                            raise TypeError('Only 1 tip4p model currently supported')

                kspace_style = None
                if tip4p_z is not None:
                    f.write('pair_style lj/cut/tip4p/long {o} {h} {b} {a} {z} {cut}\n'.format(
                        o=type_map['atom']['opc_OW'],
                        h=type_map['atom']['opc_HW'],
                        b=type_map['bond']['opc_OW','opc_HW'],
                        a=type_map['angle']['opc_HW','opc_OW','opc_HW'],
                        z=tip4p_z,
                        cut=cutoff))
                    kspace_style = "tip4p"
                else:
                    f.write('pair_style lj/cut/coul/long {cut}\n'.format(cut=cutoff))
                    kspace_style = "pppm"
                f.write('pair_modify mix geometric tail yes\n')
                f.write('special_bonds lj/coul 0.0 0.0 0.5\n')

                if kspace_style is not None:
                    f.write('kspace_style {} 1.0e-4\n'.format(kspace_style))
            else:
                if sec is not None:
                    f.write('{}_style {}\n'.format(p,sec))

    # write shake constraints on hydrogen (deduce using mass)
    with open(filename + '.ff.shake','w') as f:
        f.write('# generated using mbuild {} / foyer {}\n\n'.format(mb.__version__,foyer.__version__))

        shake_bonds = []
        for (i,j),b in type_map['bond'].items():
            if i in h_types or j in h_types:
                shake_bonds.append(b)
        fix_shake = 'fix shake all shake 1.0e-4 20 0 b ' + ' '.join(shake_bonds)
        if 'opc' in potential_map['pair']:
            fix_shake += 'a ' + type_map['angle']['opc_HW','opc_OW','opc_HW']
        f.write(fix_shake + '\n')
