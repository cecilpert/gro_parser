import re

GRO_LINE_REGEX = '^([\d ]{5})([\w ]{5})([\w ]{5})([\d ]{5})([\d. -]{8})([\d. -]{8})([\d. -]{8})([\d. -]{8})([\d. -]{8})([\d. -]{8})'

class GroSystem:
    def __init__(self, gro_file: str, top_file:str):
        
        self.name = ''
        self.box_vectors = [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
        self._residues = []
        self._index_by_residue_name = {}
        self._index_by_resnumber = {}
        self._top_includes = []
        print(f'Read gro file : {gro_file}')
        self._parse(gro_file)
        print(f'Read top file : {top_file}') #Just to keep the includes lines
        self._parse_top(top_file)

    @property
    def residues(self):
        return iter(self._residues)

    @property
    def atoms(self):
        return [a for res in self._residues for a in res.atoms]

    def _parse(self, gro_file):
        with open(gro_file) as f:
            self.name = f.readline().strip()
            self.atom_number = int(f.readline().strip())
            prev_resnum = -1
            idx = 0
            for l in f:
                line_match = re.search(GRO_LINE_REGEX, l)
                
                if line_match:
                    resnum = int(line_match.group(1).strip())
                    resname = line_match.group(2).strip()
                    atomname = line_match.group(3).strip()
                    atomnumber = int(line_match.group(4).strip())
                    coordinates = [float(line_match.group(i).strip()) for i in range(5,8)]
                    velocities = [line_match.group(i).strip() for i in range(8,11)]

                    

                    if resnum != prev_resnum:
                        # New residue
                        residue = self.add_residue_at_the_end(resname, resnum, idx)
                        idx += 1

                    else:
                        if resname != residue.name:
                            raise GroParsingException(f'A residue changed its name : resid {resid} from {residue.name} to {resname}')
                    
                    prev_resnum = resnum
                    residue.add_atom(atomname, atomnumber, coordinates, velocities)
                    

                else:
                     # If gro is correctly formated, it's the last line, check that
                    next_line = f.readline()
                    if next_line:
                        raise GroParsingException(f'A line is not gro compatible : {l}')
                    
                    #Check it's the box vectors
                    last_line_values = l.split()
                    if len(last_line_values) < 3:
                        raise GroParsingException(f'Last line (box vectors) is not gro compatible : you need at least 3 numbers')
                    for i, value in enumerate(last_line_values):
                        try :
                            float_value = float(value)
                            self.box_vectors[i] = float_value
                        except ValueError:
                            raise GroParsingException('Last line (box vectors) is not gro compatible : not numbers')

                    #Last checking, check if the box is gromacs compatible ( v1(y)=v1(z)=v2(z)=0 ). If not, print a warning. 
                    if not (self.box_vectors[3] == self.box_vectors[4] == self.box_vectors[5] == 0):
                        print('WARNING : your box is not gromacs compatible, v1(y), v1(z) and v2(z) needs to be zero')

    def _parse_top(self, top_file):
        with open(top_file) as f:
            for l in f:
                if l.startswith('#include'):
                    self._top_includes.append(l.rstrip())

    def add_residue_at_the_end(self, resname, resnum, idx):
        residue = Residue(resnum, resname, idx, self)
        self._residues.append(residue)

        if resname not in self._index_by_residue_name:
            self._index_by_residue_name[resname] = []
        self._index_by_residue_name[resname].append(residue)

        if resnum in self._index_by_resnumber:
            raise GroParsingException('Two residues with the same number')
        self._index_by_resnumber[resnum] = residue
        
        return residue   

    def insert_residue(self, residue):
        print(f'Insert <{residue}> in {residue.idx}th position')

        #Shift part
        nb_atoms_to_shift = len(residue.atoms)
        for res in self._residues[residue.idx:]:
            res.number = res.number + 1
            self._index_by_resnumber[res.number] = res
            res.idx = res.idx + 1
            for atom in res.atoms:
                atom.number = atom.number + nb_atoms_to_shift

        tmp_new_residues = self._residues[:residue.idx] + [residue] + self._residues[residue.idx:]
        self._residues = tmp_new_residues
        self._index_by_residue_name[residue.name].append(residue)
        self._index_by_resnumber[residue.number] = residue
        
    def get_residue_by_idx(self, idx) : 
        try:
            return self._residues[idx]
        except IndexError:
            raise ResidueNotFound(f'in {idx}th position')
    
    def get_residue_by_number(self, number):
        try:
            return self._index_by_resnumber[number]
        except KeyError:
            raise ResidueNotFound(number)
    
    def get_residue_by_name(self, name):
        try: 
            return self._index_by_residue_name[name]
        except KeyError:
            raise ResidueNotFound(name)      

    def write_gro(self, gro_path):
        with open(gro_path, 'w') as o:
            o.write(f'{self.name}\n')
            o.write(f' {len(self.atoms)}\n')
            for res in self.residues:
                for atom in res.atoms:
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
                    for vel in atom.velocities:
                        _vel_char = len(str(vel))
                        velocities_str += ' ' * (8 - _vel_char) + str(vel)
                    o.write(f'{resnumber}{resname}{atomname}{atomnum}{coords_str}{velocities_str}\n')
            o.write(f' {" ".join([str(val) for val in self.box_vectors])}')
        print(f'System gro file written in {gro_path}')
    
    def write_top(self, top_path):
        with open(top_path, 'w') as o:
            for include in self._top_includes:
                o.write(f'{include}\n')
            o.write('\n')
            o.write(f'[ system ]\n; name\n{self.name}\n\n')
            o.write(f'[ molecules ]\n; name number\n')
            for resname, residues in self._index_by_residue_name.items():
                o.write(f'{resname} {len(residues)}\n')  
        print(f'System top file written in {top_path}') 

class Residue:
    def __init__(self, number: int, name: str, idx, system):
        self.system = system
        self.name = name
        self.number = number
        self.idx = idx
        self.atoms = []

    def __repr__(self):
        return f"{self.name}{self.number} with {len(self.atoms)} atoms"
    
    def add_atom(self, name, number, coordinates, velocities):
        self.atoms.append(Atom(name, number, coordinates, velocities, self))
    
    def copy(self):
        print(f'I want to copy <{self}>')
        residues_stack = self.system.get_residue_by_name(self.name)
        print(f'Will be copied at the end of {len(residues_stack)} {self.name} stack')
        last = residues_stack[-1]
        last_atom = last.atoms[-1]
        new_number = last.number + 1
        new_idx = last.idx + 1
        new_residue = Residue(new_number, self.name, new_idx, self.system)
        new_atom_number = last_atom.number + 1
        for atom in self.atoms:
            new_coords = atom.coordinates[:]
            new_coords[0] = round(new_coords[0] + 0.21, 3) #vdw radius
            new_residue.add_atom(atom.name, new_atom_number, new_coords, atom.velocities)
            new_atom_number += 1
        
        self.system.insert_residue(new_residue)
        



        








            


class Atom:
    def __init__(self, name: str, number: int, coordinates, velocities, residue):
        self.name = name
        self.number = number
        self.coordinates = coordinates
        self.velocities = velocities
        self.residue = residue

    def __repr__(self):
        return f'{self.name} {self.number}'

class GroParsingException(Exception):
    pass

class ResidueNotFound(Exception):
    def __init__(self, residue):
        message = f'Residue {residue} not found'
        super().__init__(message)