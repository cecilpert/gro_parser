import re
from .base_logger import logger
import logging
from .atom import Bond, Angle, Dihedral
from .residue import Residue
from .utils import add_to_gro_number

logger.setLevel(logging.INFO)

GRO_LINE_REGEX = '^([\d ]{5})([\w ]{5})([\w ]{5})([\d ]{5})([\d. -]{8})([\d. -]{8})([\d. -]{8})(([\d. -]{8})([\d. -]{8})([\d. -]{8}))?'

class GroSystem:
    def __init__(self, gro_file: str, top_file:str = None, itp_file: str = None):
        self.gro_file = gro_file
        self.top_file = top_file
        self.itp_file = itp_file
        self.name = ''
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
        logger.info(f'Load gromacs system with gro_parser')
        logger.debug(f'Read gro file : {gro_file}')
        self._parse(gro_file)
        logger.debug(f'Read top file : {top_file}') #Just to keep the includes lines
        self.is_martini = False
        if top_file:
            self._parse_top(top_file)
        if itp_file:
            self._parse_itp(itp_file)
        self._max_x = None
        self._max_y = None
        self._max_z = None
        self._min_x = None
        self._min_y = None
        self._min_z = None

        

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

    def _parse(self, gro_file):
        with open(gro_file) as f:
            #Read the system name from the first line
            self.name = f.readline().strip()
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
        
        self._try_to_populate_atom_with_itp()
        
        self._register_itp_relation('bonds', ['i', 'j'])
        self._register_itp_relation('angles', ['i', 'j', 'k'])
        self._register_itp_relation('dihedrals', ['i', 'j', 'k', 'l'])

        self._top_includes.append(f'#include "{itp_file}"')


    def _try_to_populate_atom_with_itp(self):
        if not self.itp_info['atoms']:
            logger.warn('No atoms category in itp file, impossible to load itp atoms information into system')
            return
        atoms_info = self.itp_info['atoms']
        if not atoms_info[0].startswith(';'):
            logger.warn('No header line to identify columns in itp atoms part, impossible to load itp atoms information into system')
            return
        header_cats = atoms_info[0].lstrip(';').split()
        try:
            itp_atom_name_idx = header_cats.index('atom')
        except ValueError:
            logger.warn('No atom column in itp atoms part, impossible to load itp atoms information into system')
            return

        if not 'atom' in header_cats:
            logger.warn('No atom column in itp atoms part, impossible to load itp atoms information into system')
            return
        atoms_line = atoms_info[1:]
        if len(atoms_line) != len(self.atoms):
            logger.warn('Not the same number of atoms in system and itp file, impossible to load itp atoms information into system')
            return
        for idx, atom in enumerate(self.atoms):
            itp_atom_name = atoms_line[idx].split()[itp_atom_name_idx]
            if not itp_atom_name == atom.name:
                logger.warn(f'Correspondance error with the {idx}th atom between itp and gro. Names are different : {itp_atom_name} in itp, {atom.name} in gro. Impossible to load itp atoms information into system')
                return
        
        for idx, atom in enumerate(self.atoms):
            atom.itp_info = {}
            for idx_header, header in enumerate(header_cats):
                atom.itp_info[header] = atoms_line[idx].split()[idx_header]

        logger.info("itp file has been used to register more informations on system's atoms")
                            
    def _register_itp_relation(self, part_to_check, col_to_check):
        if part_to_check not in ['angles', 'bonds', 'dihedrals']:
            logger.warn(f'{part_to_check} is not supported yet to be registered into system')
        
        if not part_to_check in self.itp_info:
            logger.warn(f'No {part_to_check} part in itp. Impossible to register {part_to_check} into system.')
            return
        if not self.itp_info[part_to_check][0].startswith(';'):
            logger.warn(f'No header line to identify columns in itp {part_to_check} part. Impossible to register {part_to_check} into system.')
            return
        
        header = self.itp_info[part_to_check][0].lstrip(';').split()
        lines = self.itp_info[part_to_check][1:]

        if part_to_check == "bonds":
            if not 'length' in header:
                logger.warn('No length column in bonds part. Impossible to register bonds into system.')
                return
            length_idx = header.index('length')

        for col in col_to_check:
            if not col in header:
                logger.warn(f'No {col} column in bonds part. Impossible to register {part_to_check} into system.')
                return
            
        col_idx = [header.index(col) for col in col_to_check]
        
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
            
            if part_to_check == 'bonds':
                link_obj = Bond(*atoms, float(line[length_idx]), comment = register_comment)
                for atom in atoms:
                    atom.register_bond(link_obj)
            
            if part_to_check == "angles":
                link_obj = Angle(*atoms, comment = register_comment)
                for atom in atoms:
                    atom.register_angle(link_obj)

            if part_to_check == "dihedrals":
                link_obj = Dihedral(*atoms, comment = register_comment)
                for atom in atoms:
                    atom.register_dihedral(link_obj)

            for idx, col in enumerate(header):
                link_obj.itp_info[col] = line[idx]
        
        logger.info(f'itp file has been used to register {part_to_check} into system')

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
        with open(top_path, 'w') as o:
            for include in self._top_includes:
                o.write(f'{include}\n')
            o.write('\n')
            o.write(f'[ system ]\n; name\n{self.name}\n\n')
            o.write(f'[ molecules ]\n; name number\n')
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
                    o.write(f'{residue_stack[0].name} {len(residue_stack)}\n')  
        logger.info(f'System top file written in {top_path}') 

    def write_itp(self, itp_path):
        should_be_handled_part = ['atoms', 'bonds', 'angles', 'dihedrals']

        with open(itp_path, 'w') as itp:
            for part in self.itp_info:
                if part == 'header':
                    itp.write('\n'.join(self.itp_info['header']) + '\n')
                    itp.write('; system loaded with gro_parser\n;\n')
                elif part in should_be_handled_part:
                    itp.write(f'[ {part} ]\n')
                    header = self.itp_headers[part]
                    itp.write(f'; {"\t".join(header)}\n')
                    
                    if part == 'atoms':
                        for atom in self.atoms:
                            itp.write("\t".join(atom.itp_info[col] for col in header) + '\n')
                    
                    if part == 'bonds':
                       written_bonds = []
                       for atom in self.atoms:
                           for bond in atom.bonds:
                                if bond not in written_bonds: 
                                    line = "\t".join(bond.itp_info[col] for col in header) + '\n'
                                    if bond.comment : 
                                        line = ";" + line
                                    itp.write(line)
                                    written_bonds.append(bond)
                    if part == 'angles':
                        written_angles = []
                        for atom in self.atoms:
                           for angle in atom.angles:
                                if angle not in written_angles: 
                                    line = "\t".join(angle.itp_info[col] for col in header) + '\n'
                                    if angle.comment:
                                        line = ";" + line
                                    itp.write(line)
                                    written_angles.append(angle)

                    if part == "dihedrals":
                        written_dihedrals = []
                        for atom in self.atoms:
                           for d in atom.dihedrals:
                                if d not in written_dihedrals: 
                                    line = "\t".join(str(d.itp_info.get(col, ' ')) for col in header)
                                    if d.comment:
                                        line = ";" + line
                                    if d.comment_str:
                                        line = line + " ;" + d.comment_str
                                    itp.write(line + "\n")
                                    written_dihedrals.append(d)
                            
                else:
                    itp.write(f'[ {part} ]\n')
                    itp.write('\n'.join(self.itp_info[part]) + '\n')

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