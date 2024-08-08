from .atom import Atom
from .base_logger import logger
from .utils import add_to_gro_number
import numpy as np

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
        atom = Atom(name, number, coordinates, velocities, self)
        self.atoms.append(atom)
        self.system._register_atom_in_index(atom)
        return atom
    
    def add_atom_object(self, atom_object):
        self.atoms.append(atom_object)
        self.system._register_atom_in_index(atom_object)

    @property
    def min_x(self):
        """
        Property: min_x

        Get the minimum X-coordinate among all atoms in the residue.
        """
        return min([atom.coordinates[0] for atom in self.atoms])
    
    @property
    def max_x(self):
        """
        Property: max_x

        Get the maximum X-coordinate among all atoms in the residue.
        """
        return max([atom.coordinates[0] for atom in self.atoms])
    
    @property
    def min_y(self):
        """
        Property: min_y

        Get the minimum Y-coordinate among all atoms in the residue.
        """
        return min([atom.coordinates[1] for atom in self.atoms])
    
    @property
    def max_y(self):
        """
        Property: max_y

        Get the maximum Y-coordinate among all atoms in the residue.
        """
        return max([atom.coordinates[1] for atom in self.atoms])
    
    @property
    def min_z(self):
        """
        Property: min_z

        Get the minimum Z-coordinate among all atoms in the residue.
        """
        return min([atom.coordinates[2] for atom in self.atoms])
    
    @property
    def max_z(self):
        """
        Property: max_z

        Get the maximum Z-coordinate among all atoms in the residue.
        """
        return max([atom.coordinates[2] for atom in self.atoms])
    
    def copy(self, offset=0.21):
        """
        Method: copy

        Create a copy of the Residue object.

        Returns:
        - None

        This method creates a new Residue object that is an updated copy of the current Residue object.
        The new Residue object has the same 'name' and 'system', but the 'number' and 'idx' attributes
        are updated based to put it at the end of the stack of residues with same name.

        Additionally, all Atom objects in the current residue are duplicated in the new residue. The
        coordinates of each atom are adjusted by adding 0.21 to the X-coordinate. This adjustment is
        meant to represent the van der Waals (vdW) radius of the atom.

        The new Residue object is not connected to the original Residue or its Atoms, meaning that
        any modifications made to the new Residue or its Atoms will not affect the original Residue
        or its Atoms, and vice versa.

        Example Usage:
        Suppose you have a 'residue' instance of the Residue class representing a molecular residue
        in a molecular system. You can create a copy of the residue like this:
        residue.copy()
        """
        logger.debug(f'I want to copy <{self}>')
        residues_stack = self.system.get_residue_by_name(self.name)
        logger.debug(f'Will be copied at the end of {len(residues_stack)} {self.name} stack')
        # Determine the 'number' and 'idx' for the new Residue based on the last residue with the same name.
        last = residues_stack[-1]
        last_atom = last.atoms[-1]

        new_number = add_to_gro_number(last.number, 1)
        new_idx = last.idx + 1

        # Create a new Residue object with updated 'number' and 'idx'.
        new_residue = Residue(new_number, self.name, new_idx, self.system)
        new_atom_number = add_to_gro_number(last_atom.number, 1)

        # Duplicate all Atom objects in the current residue and associate them with the new residue, shift their number 
        # The new Atom objects are not connected to the original Atom objects.
        for atom in self.atoms:
            new_coords = atom.coordinates[:]
            new_coords[0] = round(new_coords[0] + offset, 3) # Adjust X-coordinate by adding 0.21 (vdW radius).
            # Create a new Atom object and add it to the new residue's 'atoms' list.
            new_residue.add_atom(atom.name, new_atom_number, new_coords, atom.velocities)
            new_atom_number = add_to_gro_number(new_atom_number, 1)
        
        # Insert the new Residue into the molecular system using the 'insert_residue' method.
        self.system.insert_residue(new_residue)
        return new_residue

    def change_coordinates(self, new_coords):
        """
        Method: change_coordinates

        Update the coordinates of all atoms in the residue to new coordinates.

        Parameters:
        - new_coords (list): A list of new X, Y, Z coordinates to set for all atoms in the residue.

        Returns:
        - None

        Example Usage:
        residue = Residue(number=1, name="ALA", idx=0, system=my_system)
        new_coordinates = [1.0, 2.0, 3.0]
        residue.change_coordinates(new_coordinates)
        # All atoms in the residue will have their coordinates set to [1.0, 2.0, 3.0].
        """
        for atom in self.atoms:
            atom.coordinates = new_coords
    
    def change_velocity(self, new_velocity):
        """
        Method: change_velocity

        Update the velocities of all atoms in the residue to new velocities.

        Parameters:
        - new_velocity (list): A list of new velocity components (e.g., VX, VY, VZ) to set for all atoms in the residue.

        Returns:
        - None

        This method iterates through all atoms in the residue and updates their velocities to the new velocities
        provided. It allows for changing the motion state of the entire residue in the molecular system.

        Example Usage:
        residue = Residue(number=1, name="ALA", idx=0, system=my_system)
        new_velocity = [0.1, 0.2, 0.3]
        residue.change_velocity(new_velocity)
        # All atoms in the residue will have their velocities set to [0.1, 0.2, 0.3].
        """
        for atom in self.atoms:
            atom.velocities = new_velocity

    def delete(self):
        logger.debug(f'I want to delete <{self}>')
        self.system.delete_residue(self)

    def change_name(self, new_name):
        self.name = new_name
        self.system.redo_index_resname()
        self.system.redo_residue_stack_from_index_resname()

    def get_coordinates(self):
        coords = []
        for atom in self.atoms:
            coords.append(atom.coordinates)
        return np.array(coords)