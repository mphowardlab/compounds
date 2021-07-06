import foyer
import mbuild as mb
import numpy as np

from . import molecule
from .collections import SymmetricMapping

def save_lammps(compound, filename, cutoff, kspace_tol=1e-4, shake_tol=1e-4, **kwargs):
    if '.lammps' not in filename and '.lmp' not in filename:
        raise ValueError('File must be LAMMPS data format .lammps or .lmp')

    if 'atom_style' in kwargs or 'unit_style' in kwargs:
        raise KeyError('atom_style and unit_style are fixed as full and real')

    # make a forcefield object
    ff_kwargs = {}
    if 'forcefield_name' in kwargs:
        ff_kwargs['name'] = kwargs['forcefield_name']
    if 'forcefield_files' in kwargs:
        ff_kwargs['forcefield_files'] = kwargs['forcefield_files']
    ff = foyer.Forcefield(**ff_kwargs)

    # use mbuild to save to lammps format, though
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
                tip4p_model = None
                tip4p_type = None
                for c in compound.to_networkx():
                    if isinstance(c,molecule.Implicit4SiteWater):
                        if tip4p_type is None:
                            tip4p_model = {
                                'z': type(c).z,
                                'OW': None,
                                'HW': None
                                }
                            s = ff.apply(c)
                            for a in s.atoms:
                                if tip4p_model['OW'] is None and a.name == 'OW':
                                    tip4p_model['OW'] = a.type
                                elif tip4p_model['HW'] is None and a.name == 'HW':
                                    tip4p_model['HW'] = a.type
                            if tip4p_model['OW'] is None or tip4p_model['HW'] is None:
                                raise TypeError('Unable to deduce all implicit 4-site types')
                            tip4p_type = type(c)
                        elif type(c) != tip4p_type:
                            raise TypeError('Only 1 tip4p model currently supported')

                kspace_style = None
                if tip4p_model is not None:
                    ow = tip4p_model['OW']
                    hw = tip4p_model['HW']
                    f.write('pair_style lj/cut/tip4p/long {o} {h} {b} {a} {z} {cut}\n'.format(
                        o=type_map['atom'][ow],
                        h=type_map['atom'][hw],
                        b=type_map['bond'][ow,hw],
                        a=type_map['angle'][hw,ow,hw],
                        z=tip4p_model['z'],
                        cut=cutoff))
                    kspace_style = "pppm/tip4p"
                else:
                    f.write('pair_style lj/cut/coul/long {cut}\n'.format(cut=cutoff))
                    kspace_style = "pppm"

                if kspace_style is not None:
                    f.write('kspace_style {} {}\n'.format(kspace_style,kspace_tol))

                # set combining rule
                comb_rule = ff.combining_rule
                if not isinstance(comb_rule,str):
                    comb_rule = set(ff.combining_rule)
                    if len(comb_rule) != 1:
                        raise TypeError('Conflicting combining rules in forcefield.')
                    comb_rule = comb_rule[0]
                # lammps calls 'lorentz' the 'arithmetic' mixing rule
                if comb_rule == 'lorentz':
                    comb_rule = 'arithmetic'
                if comb_rule not in ('geometric','arithmetic','sixthpower'):
                    raise ValueError('Combining rule not recognized')
                f.write('pair_modify mix {} tail yes\n'.format(comb_rule))

                # set 1-4 exclusions
                f.write('special_bonds lj 0.0 0.0 {lj} coul 0.0 0.0 {coul}\n'.format(
                    lj=ff.lj14scale,
                    coul=ff.coulomb14scale))
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
        fix_shake = 'fix shake all shake {} 20 0 b '.format(shake_tol) + ' '.join(shake_bonds)
        if tip4p_model is not None:
            ow = tip4p_model['OW']
            hw = tip4p_model['HW']
            fix_shake += ' a ' + type_map['angle'][hw,ow,hw]
        f.write(fix_shake + '\n')
