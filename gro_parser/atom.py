class Atom:
    def __init__(self, name: str, number: int, coordinates, velocities, residue):
        self.name = name
        self.number = number
        self.coordinates = coordinates
        self.velocities = velocities
        self.residue = residue
        self.bonds = []
        self.angles = []
        self.dihedrals = []
        self.itp_info = {}

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
    
    def set_z(self, z):
        self.coordinates[2] = z

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

    def comment_all_dihedrals(self):
        for d in self.dihedrals:
            d.comment = True

    def change_atomtype(self, atomtype):
        if 'type' not in self.itp_info:
            raise MissingItpDescription('atom type')
        self.itp_info['type'] = atomtype

class Bond:
    def __init__(self, atom1:Atom, atom2: Atom, length: float, comment: bool = False):
        self.atom1 = atom1
        self.atom2 = atom2
        self.length = length
        self.itp_info = {}
        self.comment = comment

    def __repr__(self):
        return f'Bond btw {self.atom1.name} {self.atom2.name}, length {self.length}'
    
    def set_length(self, length):
        self.length = length
        self.itp_info['length'] = str(length)


class Angle:
    def __init__(self, atom1: Atom, atom2: Atom, atom3: Atom, comment: bool = False):
        self.atom1 = atom1 
        self.atom2 = atom2
        self.atom3 = atom3
        self.itp_info = {}
        self.comment = comment

    def __repr__(self):
        return f'Angle btw {self.atom1.name} {self.atom2.name} {self.atom3.name}'
    
    def delete(self):
        self.atom1.delete_one_angle(self)
        self.atom2.delete_one_angle(self)
        self.atom3.delete_one_angle(self)

    def set_angle_value(self, value):
        self.itp_info['angle'] = str(value)
    
class Dihedral:
    def __init__(self, atom1: Atom, atom2: Atom, atom3: Atom, atom4: Atom, comment : bool = False):
        self.atom1 = atom1 
        self.atom2 = atom2
        self.atom3 = atom3
        self.atom4 = atom4
        self.itp_info = {}
        self.comment = comment
    
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