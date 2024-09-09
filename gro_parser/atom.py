import copy

class ItpPart:
    def __init__(self, type, comment: bool = False):
        self.itp_info = {}
        self.type = type
        self.comment = comment
        self.comment_str = None

    def add_itp_column(self, column, value):
        self.itp_info[column] = value
        system = self.atom1.residue.system
        if column not in system.itp_headers[self.type]:
            system.itp_headers[self.type].append(column)

    def add_a_comment(self, comment):
        
        if self.comment_str and comment != self.comment_str:
            self.comment_str += f' {comment}'
        else:
            self.comment_str = comment
    
class Atom:
    def __init__(self, name: str, number: int, coordinates, velocities, residue):
        self.name = name
        self.number = number
        self.coordinates = coordinates
        self.velocities = velocities
        self.residue = residue
        self.atom_type = None #from itp
        self.charge_group_number = None #from itp
        self.charge = None #from itp
        self.mass = None #optionnal, from itp
        self.bonds = []
        self.angles = []
        self.dihedrals = []
        self.in_itp = False
        
    def __repr__(self):
        return f'{self.name} {self.number}'
    
    @property
    def x(self):
        return self.coordinates[0]
    
    @property
    def y(self):
        return self.coordinates[1]
    
    @property
    def z(self):
        return self.coordinates[2]
    
    @property
    def type(self):
        return self.atom_type
    
    def set_z(self, z):
        self.coordinates[2] = z

    def set_x(self, x):
        self.coordinates[0] = x

    def set_y(self, y):
        self.coordinates[1] = y

    def register_bond(self, bond):
        self.bonds.append(bond)
    
    def register_angle(self, angle):
        self.angles.append(angle)
    
    def register_dihedral(self, dihedral):
        self.dihedrals.append(dihedral)
    
    def delete_all_angles(self):
        copied_angles = self.angles[:]
        for angle in copied_angles:
            angle.delete()

    def delete_one_angle(self, angle):
        self.angles.remove(angle)

    def delete_all_dihedrals(self):
        copied_dihedrals = self.dihedrals[:]
        for d in copied_dihedrals:
            d.delete()

    def delete_one_dihedral(self, d):
        self.dihedrals.remove(d) 

    def delete_all_bonds(self):
        copied_bonds = self.bonds[:]
        for bond in copied_bonds:
            bond.delete()

    def comment_all_dihedrals(self):
        for d in self.dihedrals:
            d.comment = True

    def delete_all_relations(self):
        self.delete_all_bonds()
        self.delete_all_angles()
        self.delete_all_dihedrals()

    def change_atomtype(self, atomtype):
        self.atom_type = atomtype

    def copy(self):
        new_atom = copy.deepcopy(self)
        new_atom.residue = self.residue
        return new_atom
    
    def set_name(self, name):
        self.name = name
        
    def set_number(self, number):
        self.number = number
    
    def set_charge_group_number(self, number):
        self.charge_group_number = number

    def increase_x(self, by):
        self.coordinates[0] += by
    
    def increase_y(self, by):
        self.coordinates[1] += by

    def get_neighbors_through_bond(self):
        neighbors = []
        for b in self.bonds:
            if b.atom1 == self:
                neighbors.append(b.atom2)
            else:
                neighbors.append(b.atom1)
        return neighbors
    
    def get_bond(self, other_atom):
        for b in self.bonds:
            if b.atom1 == other_atom or b.atom2 == other_atom:
                return b
        

class Bond(ItpPart):
    def __init__(self, atom1:Atom, atom2: Atom, funct, other_params = {}, comment: bool = False):
        super().__init__('bonds', comment)
        self.atom1 = atom1
        self.atom2 = atom2
        self.funct = str(funct)
        self.other_params = other_params

    def __repr__(self):
        return f'Bond btw {self.atom1.name} {self.atom2.name}'
    
    def set_length(self, length):
        if not 'length' in self.other_params:
            raise InvalidBondChange(self.funct, 'length')
        
        self.other_params['length'] = length

    def delete(self):
        self.atom1.bonds.remove(self)
        self.atom2.bonds.remove(self)


    @property
    def length(self):
        if not 'length' in self.other_params:
            raise InvalidParameter(self.funct, 'length', 'bond')
        return float(self.other_params['length'])


class Angle(ItpPart):
    def __init__(self, atom1: Atom, atom2: Atom, atom3: Atom, funct = None, other_params = {}, comment: bool = False):
        super().__init__('angles', comment)
        self.atom1 = atom1 
        self.atom2 = atom2
        self.atom3 = atom3
        self.funct = str(funct)
        self.other_params = other_params

    def __repr__(self):
        return f'Angle ({self.funct}) btw {self.atom1.name} {self.atom2.name} {self.atom3.name}'
    
    def delete(self):
        self.atom1.delete_one_angle(self)
        self.atom2.delete_one_angle(self)
        self.atom3.delete_one_angle(self)

    def set_angle_value(self, value):
        if not 'angle' in self.other_params:
            raise InvalidAngleChange(self.funct, 'angle')
        self.other_params['angle'] = value

    
class Dihedral(ItpPart):
    def __init__(self, atom1: Atom, atom2: Atom, atom3: Atom, atom4: Atom, funct, other_params = {}, comment : bool = False):
        super().__init__('dihedrals', comment)
        self.atom1 = atom1 
        self.atom2 = atom2
        self.atom3 = atom3
        self.atom4 = atom4
        self.funct = str(funct)
        self.other_params = other_params
       
    def __repr__(self):
        return f'Dihedral btw {self.atom1.name} {self.atom2.name} {self.atom3.name} {self.atom4.name}'
    
    def delete(self):
        self.atom1.delete_one_dihedral(self)
        self.atom2.delete_one_dihedral(self)
        self.atom3.delete_one_dihedral(self)
        self.atom4.delete_one_dihedral(self)


class MissingItpDescription(Exception):
    def __init__(self, missing_thing):
        message = f'{missing_thing} is not described in itp file'
        super().__init__(message)

class InvalidAngleChange(Exception):
    def __init__(self, funct, param):
        message = f"Impossible to modify {param} on angle funct {funct}, this parameter doesn't exist"
        super().__init__(message)

class InvalidBondChange(Exception):
    def __init__(self, funct, param):
        message = f"Impossible to modify {param} on bond funct {funct}, this parameter doesn't exist"
        super().__init__(message)

class InvalidParameter(Exception):
    def __init__(self, funct, param, type):
        message = f"Impossible to get {param} on {type} funct {funct}, this parameter doesn't exist"
        super().__init__(message)


# ITP_PARAMETERS = {
#     'dihedrals' :  {
#         (4, 'periodic improper dihedrals') : ['angle', 'k', 'multiplicity']
#     }
# }