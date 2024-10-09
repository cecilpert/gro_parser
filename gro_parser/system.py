import re
from .base_logger import logger
import logging
from .atom import Bond, Angle, Dihedral
from .residue import Residue
from .utils import add_to_gro_number, get_relation_params
import json
from pathlib import Path
import os

source_path = Path(__file__).resolve()
source_dir = source_path.parent


logger.setLevel(logging.INFO)

GRO_LINE_REGEX = '^([\d ]{5})([\w ]{5})([\w+|\- ]{5})([\d ]{5})([\d. -]{8})([\d. -]{8})([\d. -]{8})(([\d. -]{8})([\d. -]{8})([\d. -]{8}))?'

ITP_DESC = json.load(open(f'{source_dir}/itp_description.json'))

class GroSystem:
    def __init__(self, gro_file: str, top_file:str = None, itp_file: str = None, name = None):
        self.gro_file = gro_file
        self.top_file = top_file
        self.itp_file = itp_file
        self.name = name
        self.box_vectors = [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
        self._residues = []
        self._index_by_residue_name = {}
        self._index_by_resnumber = {}
        self._index_by_atomnumber = {}
        self._top_includes = []
        self._residues_stack = []
        self._index_residues_stack_by_name = {}
        self.itp_info = {}
        self.itp_headers = {}
       
        self.is_martini = False
        
        self._max_x = None
        self._max_y = None
        self._max_z = None
        self._min_x = None
        self._min_y = None
        self._min_z = None
        self.mol_name_in_itp = None

        logger.info(f'Load gromacs system with gro_parser')
        logger.debug(f'Read gro file : {gro_file}')
        self._parse(gro_file)
        if top_file:
            self._parse_top(top_file)
        if itp_file:
            self._parse_itp(itp_file)

        

    @property
    def residues(self):
        """
        Provides an iterator over the Residue objects in the molecular system.
        """
        return iter(self._residues)
    
    @property
    def residue_names(self):
        """
        Retrieve a set of unique residue names present in the molecular system.'''
        """
        return set([res.name for res in self.residues])

    @property
    def atoms(self):
        """
        Returns a list of all Atom objects in the molecular system.
        """
        return [a for res in self._residues for a in res.atoms]
    
    @property
    def atoms_with_itp(self):
        return [a for res in self._residues for a in res.atoms if a.in_itp]

    @property
    def max_x(self):
        """
        Returns the maximum x-coordinate among all atoms in the molecular system and store it to avoid re computing if it's asked several time
        """
        if self._max_x == None:
            self._max_x = max([atom.coordinates[0] for atom in self.atoms])
        return self._max_x
    
    @property
    def min_x(self):
        """
        Returns the minimum x-coordinate among all atoms in the molecular system and store it to avoid re computing if it's asked several time
        """
        if self._min_x == None:
            self._min_x = min([atom.coordinates[0] for atom in self.atoms])
        return self._min_x
    
    @property
    def max_y(self):
        """
        Returns the maximum y-coordinate among all atoms in the molecular system and store it to avoid re computing if it's asked several time
        """
        if self._max_y == None:
            self._max_y = max([atom.coordinates[1] for atom in self.atoms])
        return self._max_y
    
    @property
    def min_y(self):
        """
        Returns the minimum y-coordinate among all atoms in the molecular system and store it to avoid re computing if it's asked several time
        """
        if self._min_y == None:
            self._min_y = min([atom.coordinates[1] for atom in self.atoms])
        return self._min_y
    
    @property
    def max_z(self):
        """
        Returns the maximum z-coordinate among all atoms in the molecular system and store it to avoid re computing if it's asked several time
        """
        if self._max_z == None:
            self._max_z = max([atom.coordinates[2] for atom in self.atoms])
        return self._max_z
    
    @property
    def min_z(self):
        """
        Returns the minimum z-coordinate among all atoms in the molecular system and store it to avoid re computing if it's asked several time
        """
        if self._min_z == None:
            self._min_z = min([atom.coordinates[2] for atom in self.atoms])
        return self._min_z

    @property
    def itp_names(self):
        return [include.lstrip('"#include ').rstrip('"') for include in self._top_includes] 
    
    @property
    def angles_by_funct(self):
        to_return = {}
        for angle in self.angles:
            if not angle.funct in to_return:
                to_return[angle.funct] = []
            to_return[angle.funct].append(angle)
        return to_return
    
    @property
    def bonds_by_funct(self):
        to_return = {}
        for b in self.bonds:
            if not b.funct in to_return:
                to_return[b.funct] = []
            to_return[b.funct].append(b)
        return to_return
    
    @property
    def dihedrals_by_funct(self):
        to_return = {}
        for d in self.dihedrals:
            if not d.funct in to_return:
                to_return[d.funct] = []
            to_return[d.funct].append(d)
        return to_return
    
    @property
    def angles(self):
        angles = []
        for a in self.atoms:
            for angle in a.angles:
                if angle not in angles:
                    angles.append(angle)
        return angles
    
    @property
    def bonds(self):
        bonds = []
        for a in self.atoms:
            for b in a.bonds:
                if b not in bonds:
                    bonds.append(b)
        return bonds

    @property
    def dihedrals(self):
        dihedrals = []
        for a in self.atoms:
            for d in a.dihedrals:
                if d not in dihedrals:
                    dihedrals.append(d)
        return dihedrals
    
    def _parse(self, gro_file):
        with open(gro_file) as f:
            #Read the system name from the first line if not provided by user
            gro_name = f.readline().strip()
            if not self.name:
                self.name = gro_name
            f.readline() #Skip second line
            prev_resnum = -1 # Initialize the previous residue number to track when we go to the next residue
            idx = 0 # Initialize an index to assign to each residue
            prev_resname = ""
            #Loop through lines from line 3
            for l in f:
                line_match = re.search(GRO_LINE_REGEX, l)
                # Extract information from the line using regular expressions
                if line_match:
                    resnum = int(line_match.group(1).strip())
                    resname = line_match.group(2).strip()
                    atomname = line_match.group(3).strip()
                    atomnumber = int(line_match.group(4).strip())
                    coordinates = [float(line_match.group(i).strip()) for i in range(5,8)]
                    velocity_group = line_match.group(8)
                    velocities = []
                    if velocity_group:
                        velocities = [line_match.group(i).strip() for i in range(9,12)]

                    if resnum != prev_resnum:
                        # New residue: create a new Residue object and add it to the system
                        if prev_resname != resname:
                            new_stack = True
                        else: 
                            new_stack = False
                        residue = self.add_residue_at_the_end(resname, resnum, idx, new_stack)
                        idx += 1

                    else:
                        # Same residue: check if the residue name matches the previous one
                        if resname != residue.name:
                            raise GroParsingException(f'A residue changed its name : resid {idx} from {residue.name} to {resname}')
                    
                    prev_resnum = resnum # Update the previous residue number for the next iteration
                    prev_resname = resname
                    #Add the atoms to the residue
                    atom = residue.add_atom(atomname, atomnumber, coordinates, velocities)
                    self._register_atom_in_index(atom)

                else:
                     # If gro is correctly formated, it's the last line, check that
                    next_line = f.readline()
                    if next_line:
                        raise GroParsingException(f'A line is not gro compatible : {l}')
                    
                    #Check if it's the box vectors
                    last_line_values = l.split()
                    if len(last_line_values) < 3:
                        raise GroParsingException(f'Last line (box vectors) is not gro compatible : you need at least 3 numbers')
                    for i, value in enumerate(last_line_values):
                        try :
                            float_value = float(value)
                            self.box_vectors[i] = float_value
                        except ValueError:
                            raise GroParsingException('Last line (box vectors) is not gro compatible : not numbers')

                    # Last checking, verify if the box is GROMACS compatible (v1(y)=v1(z)=v2(z)=0). If not, print a warning message.
                    if not (self.box_vectors[3] == self.box_vectors[4] == self.box_vectors[5] == 0):
                        logger.warn('WARNING : your box is not gromacs compatible, v1(y), v1(z) and v2(z) needs to be zero')

        self._create_residue_stack_index()

    def _parse_top(self, top_file):
        with open(top_file) as f:
            for l in f:
                if l.startswith('#include'):
                    self._top_includes.append(l.rstrip())
                    if 'martini' in l:
                        self.is_martini = True
        if self.is_martini:
            logger.debug('Your force field seems to be martini')

    def _parse_itp(self, itp_file):
        current_cat_line = 0
        with open(itp_file) as f:
            current_category = 'header'
            line_idx = 0
            for l in f:
                line_idx += 1
                l = l.strip().rstrip('\n')
                if l != '':
                    if l.startswith('['):
                        current_cat_line = line_idx
                        current_category = l.replace('[', '').replace(']', '').strip()
                    else:
                        if current_category not in self.itp_info:
                            self.itp_info[current_category] = []
                            self.itp_headers[current_category] = []
                        if l.startswith(';') and line_idx == current_cat_line + 1:
                            for col in l.lstrip(';').split():
                                self.itp_headers[current_category].append(col)
                        
                        self.itp_info[current_category].append(l.rstrip('\n'))
        
        self._load_name_from_itp()
        self._try_to_populate_atom_with_itp()
        
        
        self._register_itp_relation('bonds')
        self._register_itp_relation('angles')
        self._register_itp_relation('dihedrals')

        #self._top_includes.append(f'#include "{os.path.abspath(itp_file)}"')

    def _load_name_from_itp(self):
        # Temporary : can only load if the itp has one moleculetype
        moleculetype = [line for line in self.itp_info['moleculetype'] if not line.startswith(';')]
        self.mol_name_in_itp = moleculetype[0].split()[0]

    def _try_to_populate_atom_with_itp(self):
        if not self.itp_info['atoms']:
            logger.warn('No atoms category in itp file, impossible to load itp atoms information into system')
            return
        atoms_info = [l.split() for l in self.itp_info['atoms'] if not l.startswith(';')]
        # if not atoms_info[0].startswith(';'):
        #     logger.warn('No header line to identify columns in itp atoms part, impossible to load itp atoms information into system')
        #     return
        header_cats = ITP_DESC['atoms']
        atom_idx_idx = header_cats.index('atom_idx')
        atom_type_idx = header_cats.index('atom_type')
        charge_group_number_idx = header_cats.index('charge_group_number')
        charge_idx = header_cats.index('charge')
        mass_idx = header_cats.index('mass')
        
        # try:
        #     itp_atom_name_idx = header_cats.index('atom')
        # except ValueError:
        #     logger.warn('No atom column in itp atoms part, impossible to load itp atoms information into system')
        #     return

        # if not 'atom' in header_cats:
        #     logger.warn('No atom column in itp atoms part, impossible to load itp atoms information into system')
        #     return
        # atoms_line = atoms_info[1:]

        # if len(atoms_info) != len(self.atoms):
        #     logger.warn('Not the same number of atoms in system and itp file, impossible to load itp atoms information into system')
        #     return
        # for idx, atom in enumerate(self.atoms):
        #     itp_atom_name = atoms_line[idx].split()[itp_atom_name_idx]
        #     if not itp_atom_name == atom.name:
        #         logger.warn(f'Correspondance error with the {idx}th atom between itp and gro. Names are different : {itp_atom_name} in itp, {atom.name} in gro. Impossible to load itp atoms information into system')
        #         return

        
        for itp_atom in atoms_info:
            try:
                atom_idx = int(itp_atom[atom_idx_idx])
            except ValueError:
                logger.error("First column of [ atoms ] should be atom number. It's not a int. Check your itp file")
                return
            

            system_atom = self.get_atom_by_number(atom_idx)
            if not check_itp_atom(itp_atom, system_atom):
                logger.error('# Inconsistency between gro and itp. Impossible to load atoms part')
                return
            
            #more information to add to system atoms
            atom_type = itp_atom[atom_type_idx]
            system_atom.atom_type = atom_type

            try:
                charge_group_number = itp_atom[charge_group_number_idx]
            except ValueError:
                raise ItpParsingException(f'[atoms] : charge_group_number (col {charge_group_number_idx}) should be a int.')
            
            system_atom.charge_group_number = charge_group_number

            try:
                charge = float(itp_atom[charge_idx])
            except ValueError:
                raise ItpParsingException(f'[atoms] : charge_group_number (col {charge_idx}) should be a float.')

            system_atom.charge = charge
            try: 
                mass = float(itp_atom[mass_idx])
                system_atom.mass = mass
            except ValueError:
                print(f'[atoms] : no mass column')
            
            

            system_atom.in_itp = True
            system_atom.residue.in_itp = True


            
        logger.info("itp file has been used to register more informations on system's atoms")
                            
    def _register_itp_relation(self, part_to_check):
        if part_to_check not in ['angles', 'bonds', 'dihedrals']:
            logger.warn(f'{part_to_check} is not supported yet to be registered into system')
        
        if not part_to_check in self.itp_info:
            logger.warn(f'No {part_to_check} part in itp. Impossible to register {part_to_check} into system.')
            return
       
        
        header = ITP_DESC[part_to_check]
        lines = [l for l in self.itp_info[part_to_check] if not l.startswith(';')]
        
        if part_to_check == "bonds":
            atom_ids = ["i", "j"]
        
        if part_to_check == "angles":
            atom_ids = ["i", "j", "k"]
        
        if part_to_check == "dihedrals":
            atom_ids = ["i", "j", "k", "l"]

        
        col_idx = [header['default'].index(col) for col in atom_ids]
        funct_idx = header['default'].index('funct')
        
        for line in lines:
            if line.startswith(';'):
                register_comment = True
            else: 
                register_comment = False

            line = line.lstrip(';').split()
            
            atoms = []
            for atom_idx in col_idx:
                try:
                    atom = self.get_atom_by_number(int(line[atom_idx]))
                    atoms.append(atom)
                except AtomNotFound:
                    logger.warn(f'Atom {line[atom_idx]} not found in system. Impossible to register this {part_to_check} into system.')
                    
                    return
            

            funct = line[funct_idx]
            other_params = get_relation_params(header, funct, line, part_to_check)

            if part_to_check == 'bonds':
                link_obj = Bond(*atoms, funct, other_params, comment = register_comment)
                for atom in atoms:
                    atom.register_bond(link_obj)
            
            if part_to_check == "angles":
                link_obj = Angle(*atoms, funct, other_params, comment = register_comment)
                for atom in atoms:
                    atom.register_angle(link_obj)

            if part_to_check == "dihedrals":
                link_obj = Dihedral(*atoms, funct, other_params, comment = register_comment)
                for atom in atoms:
                    atom.register_dihedral(link_obj)

        
        logger.info(f'itp file has been used to register {part_to_check} into system')

    def add_itps(self,itp_paths):
        for itp in itp_paths:
            self._top_includes.append(f'#include "{os.path.abspath(itp)}"')

    def add_residue_at_the_end(self, resname, resnum, idx, new_stack):
        #Create and add the residue
        residue = Residue(resnum, resname, idx, self)
        self._residues.append(residue)
        if new_stack:
            self._residues_stack.append([residue])
        else:
            self._residues_stack[-1].append(residue)

         # Update the internal dictionaries for efficient residue retrieval
        if resname not in self._index_by_residue_name:
            self._index_by_residue_name[resname] = []
        self._index_by_residue_name[resname].append(residue)

        # Ensure that the residue number is unique and not already present in the system
        if resnum not in self._index_by_resnumber:
            self._index_by_resnumber[resnum] = []
            
        self._index_by_resnumber[resnum].append(residue)
        
        return residue   

    def insert_residue(self, residue):
        logger.debug(f'Insert <{residue}> in {residue.idx}th position')

        # Calculate the number of atoms in the new residue to shift atom numbers
        nb_atoms_to_shift = len(residue.atoms)

        # Shift existing residues and their atom numbers to accommodate the new residue
        for res in self._residues[residue.idx:]:
            self._clear_from_index_resnum(res)
            res.number = add_to_gro_number(res.number, 1)
            self._add_to_index_resnum(res)
            res.idx = res.idx + 1

            for atom in res.atoms:
                atom.number = add_to_gro_number(atom.number, nb_atoms_to_shift)

        # Insert the new residue into the molecular system list
        tmp_new_residues = self._residues[:residue.idx] + [residue] + self._residues[residue.idx:]
        self._residues = tmp_new_residues

        # Add to last residue stack
        if not residue.name in self._index_residues_stack_by_name:
            self._index_residues_stack_by_name[residue.name] = [[]] 
        self._index_residues_stack_by_name[residue.name][-1].append(residue)

        # Update the internal dictionaries to include the new residue for efficient retrieval
        if not residue.name in self._index_by_residue_name:
            self._index_by_residue_name[residue.name] = []
        self._index_by_residue_name[residue.name].append(residue)
        self._add_to_index_resnum(residue)
        
        self._index_by_atomnumber = {}
        for a in self.atoms:
            self._register_atom_in_index(a)

    def delete_residue(self, residue):
        self._residues.remove(residue)
        self._index_by_residue_name[residue.name].remove(residue)
        for residue_stack in self._index_residues_stack_by_name[residue.name]:
            if residue in residue_stack:
                residue_stack.remove(residue)
        
    def get_residue_by_idx(self, idx) : 
        """
        Retrieve a Residue object from the molecular system based on its index.
        """
        try:
            return self._residues[idx]
        except IndexError:
            raise ResidueNotFound(f'in {idx}th position')
    
    def get_residue_by_number(self, number):
        """
        Retrieve a Residue object from the molecular system based on its number.
        """
        try:
            return self._index_by_resnumber[number]
        except KeyError:
            raise ResidueNotFound(number)
    
    def get_residue_by_name(self, name):
        """
        Retrieve a list of Residue objects from the molecular system based on residue name.
        """
        try: 
            return self._index_by_residue_name[name]
        except KeyError:
            raise ResidueNotFound(name)      

    def get_atom_by_number(self, number):
        try:
            return self._index_by_atomnumber[number]
        except KeyError: 
            raise AtomNotFound(number)

    def _add_to_index_resnum(self, res):
        if res.number not in self._index_by_resnumber:
            self._index_by_resnumber[res.number] = []
        self._index_by_resnumber[res.number].append(res)

    def _clear_from_index_resnum(self, res):
        self._index_by_resnumber[res.number].remove(res)

    def _create_residue_stack_index(self):
        for residue_stack in self._residues_stack:
            name = residue_stack[0].name
            if name not in self._index_residues_stack_by_name:
                self._index_residues_stack_by_name[name] = []
            self._index_residues_stack_by_name[name].append(residue_stack)           

    def select_residues_from_name(self, name, select_func, *kwargs):
        residues_to_look_at = self.get_residue_by_name(name)
        to_return = []
        for res in residues_to_look_at:
            if select_func(res, *kwargs):
                to_return.append(res)
        return to_return
    
    def select_residues(self, select_func, *kwargs):
        to_return = []
        for res in self.residues:
            if select_func(res, *kwargs):
                to_return.append(res)
        return to_return

    def change_itp_molname(self, name):
        self.mol_name_in_itp = name

    def write_gro(self, gro_path):
        with open(gro_path, 'w') as o:
            o.write(f'{self.name}\n') # Write the system name as the first line in the .gro file
            o.write(f' {len(self.atoms)}\n') # Write the total number of atoms as the second line
            for res in [res for res_stack in self._residues_stack for res in res_stack]: # Iterate through all residues in the system
                for atom in res.atoms: # Iterate through all atoms in each residue
                    # Format the atom data for the .gro file
                    _resnum_char = len(str(res.number))
                    resnumber = ' ' * (5-_resnum_char) + str(res.number)
                    _resname_char = len(res.name)
                    resname = res.name + ' ' * (5 - _resname_char)
                    _atomname_char = len(atom.name)
                    atomname = ' ' * (5-_atomname_char) + atom.name
                    _atomnum_char = len(str(atom.number))
                    atomnum = ' ' * (5 - _atomnum_char) + str(atom.number)
                    coords_str = ''
                    for coord in atom.coordinates:
                        str_coord = "{:.3f}".format(coord)
                        _coord_char = len(str_coord)
                        coords_str += ' ' * (8 - _coord_char) + str_coord
                    velocities_str = ''
                    if atom.velocities:
                        for vel in atom.velocities:
                            _vel_char = len(str(vel))
                            velocities_str += ' ' * (8 - _vel_char) + str(vel)
                    # Write the formatted atom data to the .gro file
                    o.write(f'{resnumber}{resname}{atomname}{atomnum}{coords_str}{velocities_str}\n')

            o.write(f' {" ".join([str(val) for val in self.box_vectors])}\n') # Write box vectors at the end
        logger.info(f'System gro file written in {gro_path}')
    
    def write_top(self, top_path):
        self.redo_index_resname()
        self.redo_residue_stack_from_index_resname()
        with open(top_path, 'w') as o:
            for include in self._top_includes:
                o.write(f'{include}\n')
            o.write('\n')
            o.write(f'[ system ]\n; name\n{self.name}\n\n')
            o.write(f'[ molecules ]\n; name number\n')
            itp_mol_written = False
            for residue_stack in self._residues_stack:
                if residue_stack[0].name == 'ION' and self.is_martini:
                    logger.debug('Your system is martini with ion, ION will be replaced by atom name on the top file')
                    atom_names = set([atom.name for res in residue_stack for atom in res.atoms])
                    if len(atom_names) == 1:
                        o.write(f'{atom_names.pop()} {len(residue_stack)}\n')  
                    else:
                        ion_count = {}
                        for res in residue_stack:
                            for atom in res.atoms:
                                if atom.name not in ion_count:
                                    ion_count[atom.name] = 0
                                ion_count[atom.name] += 1
                        for ion, nb in ion_count.items():
                            o.write(f'{ion} {nb} \n')
                else:
                    if residue_stack[0].in_itp:
                        if itp_mol_written:
                            continue
                        else: 
                            o.write(f'{self.mol_name_in_itp} 1\n') 
                            itp_mol_written = True
                            #WARNING TEMP
                    else:
                        o.write(f'{residue_stack[0].name} {len(residue_stack)}\n')  
        logger.info(f'System top file written in {top_path}') 

    def write_itp(self, itp_path):
        should_be_handled_part = ['atoms', 'angles', 'bonds', 'dihedrals']

        with open(itp_path, 'w') as itp:
            for part in self.itp_info:
                if part == 'header':
                    itp.write('\n'.join(self.itp_info['header']) + '\n')
                    itp.write('; system loaded with gro_parser\n;\n')
                elif part == 'moleculetype':
                    itp.write(f'[ {part} ]\n')
                    itp.write(f'{self.mol_name_in_itp} 1\n\n')
                elif part in should_be_handled_part:
                    itp.write(f'[ {part} ]\n')
                    if part == 'atoms':
                        header = ITP_DESC[part]
                        header_str = "\t".join(header)
                        itp.write('; ' + header_str + "\n")
                        for atom in self.atoms_with_itp:
                            itp.write(f'{atom.number}\t{atom.atom_type}\t{atom.residue.number}\t{atom.residue.name}\t{atom.name}\t{atom.charge_group_number}\t{atom.charge}')
                            if atom.mass != None:
                                itp.write(f'\t{atom.mass}')
                            itp.write('\n')
                        itp.write('\n')
                    else:
                        header_default = ITP_DESC[part]['default']
                        header_default_str = "\t".join(header_default)

                        if part == 'bonds':
                            browsed_relations = self.bonds_by_funct.items()
                        elif part == 'angles':
                            browsed_relations = self.angles_by_funct.items()
                        elif part == "dihedrals":
                            browsed_relations = self.dihedrals_by_funct.items()
                            

                        written_relation = []
                        for rel_funct, relations in browsed_relations:
                            header_funct = ITP_DESC[part][str(rel_funct)]
                            header_funct_str = "\t".join(header_funct)
                            itp.write(f";{header_default_str}\t{header_funct_str}\n")
                            for rel in relations:
                                if rel not in written_relation: 
                                    if part == "bonds":
                                        line = f'{rel.atom1.number}\t{rel.atom2.number}\t{rel.funct}'
                                    elif part == "angles":
                                        line = f'{rel.atom1.number}\t{rel.atom2.number}\t{rel.atom3.number}\t{rel.funct}'
                                    elif part == "dihedrals":
                                        line = f'{rel.atom1.number}\t{rel.atom2.number}\t{rel.atom3.number}\t{rel.atom4.number}\t{rel.funct}'
                                    for other_param in header_funct['params']:
                                        line += f'\t{rel.other_params[other_param]}'
                                    if rel.comment:
                                        line = ";" + line
                                    if rel.comment_str:
                                        line = line + " ;" + rel.comment_str
                                    itp.write(line + "\n")
                                    written_relation.append(rel)     
                        itp.write('\n')      
                else:
                    itp.write(f'[ {part} ]\n')
                    itp.write('\n'.join(self.itp_info[part]) + '\n\n')

        logger.info(f'System itp file written in {itp_path}')

    def make_change_on_residues(self, change_func, *kwargs):
        for res in self.residues:
            change_func(res, *kwargs)

    def select_atoms(self, select_func, *kwargs):
        to_return = []
        for atom in self.atoms:
            if select_func(atom, *kwargs):
                to_return.append(atom)
        return to_return
    
    def _register_atom_in_index(self, atom):
        self._index_by_atomnumber[atom.number] = atom
    
    def _clear_from_index_atomnum(self, atom):
        del self._index_by_atomnumber[atom.number]

    def redo_index_resname(self):
        self._index_by_residue_name = {}
        for res in self.residues:
            if res.name not in self._index_by_residue_name: 
                self._index_by_residue_name[res.name] = []
            self._index_by_residue_name[res.name].append(res)

    def redo_residue_stack_from_index_resname(self):
        self._residues_stack = []
        for residues in self._index_by_residue_name.values():
            self._residues_stack.append(residues)



def check_itp_atom(itp_line, system_atom):
    atom_desc = ITP_DESC['atoms']
    atom_name_idx = atom_desc.index('atom_name')
    residue_name_idx = atom_desc.index('residue_name')
    residue_number_idx = atom_desc.index('residue_number')
    
    atom_name = str(itp_line[atom_name_idx])
    if atom_name != system_atom.name:
        logger.error(f'itp atom {itp_line} is not consistent with system_atom {system_atom} : different name')
        return False
    
    residue_name = str(itp_line[residue_name_idx])
    if residue_name != system_atom.residue.name:
        logger.error(f'itp atom {itp_line} is not consistent with system_atom {system_atom} : different residue name')
        return False
    
    try:
        residue_number = int(itp_line[residue_number_idx])
    except ValueError:
        raise ItpParsingException(f'residue_number (col {residue_number_idx}) is not a int')
    
    if residue_number != system_atom.residue.number:
        logger.error(f'itp atom {itp_line} is not consistent with system_atom {system_atom} : different residue number')
        return False

    return True

    

class GroParsingException(Exception):
    pass

class ResidueNotFound(Exception):
    def __init__(self, residue):
        message = f'Residue {residue} not found'
        super().__init__(message)


class AtomNotFound(Exception):
    def __init__(self, number):
        message = f'Atom {number} not found'
        super().__init__(message)

class ItpParsingException(Exception):
    def __init__(self, message):
        super().__init__(message)