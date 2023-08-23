import re

GRO_LINE_REGEX = '^([\d ]{5})([\w ]{5})([\w ]{5})([\d ]{5})([\d. -]{8})([\d. -]{8})([\d. -]{8})(([\d. -]{8})([\d. -]{8})([\d. -]{8}))?'

class GroSystem:
    def __init__(self, gro_file: str, top_file:str):
        """
        Initialize a molecular system based on GROMACS .gro and .top files.

        Parameters:
        - gro_file (str): Path to the .gro file containing atomic coordinates and velocities.
        - top_file (str): Path to the .top file containing information about the molecular topology.

        Attributes:
        - name (str): Name of the system.
        - box_vectors (list): A list of nine elements representing the simulation box vectors.
        - _residues (list): A private list to store Residue objects representing the chemical residues in the system.
        - _index_by_residue_name (dict): A private dictionary to quickly access residues by their names.
        - _index_by_resnumber (dict): A private dictionary to quickly access residues by their numbers.
        """
        self.name = ''
        self.box_vectors = [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
        self._residues = []
        self._index_by_residue_name = {}
        self._index_by_resnumber = {}
        self._top_includes = []
        self._residues_stack = []
        self._index_residues_stack_by_name = {}
        print(f'Read gro file : {gro_file}')
        self._parse(gro_file)
        print(f'Read top file : {top_file}') #Just to keep the includes lines
        self.is_martini = False
        self._parse_top(top_file)
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
        """
            Private Method: _parse

            Parse the information from the .gro file and populate the molecular system.

            Parameters:
            - gro_file (str): Path to the .gro file containing atomic coordinates and velocities.

            This private method reads the contents of the .gro file and extracts information
            about the atoms, residues, and box vectors of the molecular system. It populates the
            internal attributes such as _residues, _index_by_residue_name, and _index_by_resnumber
            to efficiently store and access data related to the residues and atoms. The method
            ensures correct association of atoms with their respective residues, maintains residue
            order, and handles the box vectors information. It also performs necessary validation
            to ensure the file is in GROMACS-compatible format. This method is called during the
            initialization of the GroSystem class to read and set up the molecular system from the
            provided .gro file.
        """

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
                    residue.add_atom(atomname, atomnumber, coordinates, velocities)
                    

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
                        print('WARNING : your box is not gromacs compatible, v1(y), v1(z) and v2(z) needs to be zero')

        self._create_residue_stack_index()

    def _parse_top(self, top_file):
        """
            Private Method: _parse_top

            Parse the .top file and store the include lines.

            Parameters:
            - top_file (str): Path to the .top file containing information about the molecular topology.

            This private method reads the contents of the .top file and extracts include lines, if any,
            from the file. Include lines in the .top file reference other topology files. These lines
            are stored in the _top_includes list, which helps keep track of any include statements that
            may need to be processed later. This method is called during the
            initialization of the GroSystem class to handle include lines from the provided .top file.
        """
        with open(top_file) as f:
            for l in f:
                if l.startswith('#include'):
                    self._top_includes.append(l.rstrip())
                    if 'martini' in l:
                        self.is_martini = True
        if self.is_martini:
            print('Your force field seems to be martini')

    def add_residue_at_the_end(self, resname, resnum, idx, new_stack):
        """
            Method: add_residue_at_the_end

            Add a new Residue to the molecular system at the end of the residue list.

            Parameters:
            - resname (str): Name of the new residue.
            - resnum (int): Number of the new residue.
            - idx (int): Index to assign to the new residue.

            Returns:
            - Residue: The newly created Residue object.
        """

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
        """
        Method: insert_residue

        Insert a new Residue into the molecular system at a specific position.

        Parameters:
        - residue (Residue): The Residue object to be inserted.

        This method inserts a new Residue object into the molecular system at a specific position.
        The method takes a Residue object as input and places it at the designated position (given by its idx attribute), updating
        the system's residue list, residue indices, and atom numbers accordingly. Residues and atoms
        following the insertion position are shifted to accommodate the new residue. The method ensures
        that the inserted residue maintains its correct position, and other residues and atoms are updated
        accordingly. The new residue is placed between two existing residues, at the end of the stack of residues with the same name,  and its index is updated to
        reflect its correct position within the system. This method is useful when inserting a new residue
        obtained from another system or when modifying the system's topology with new components.
        """
        print(f'Insert <{residue}> in {residue.idx}th position')

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
        self._index_residues_stack_by_name[residue.name][-1].append(residue)

        # Update the internal dictionaries to include the new residue for efficient retrieval
        self._index_by_residue_name[residue.name].append(residue)
        self._add_to_index_resnum(residue)

        
        
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

    # def select_residues(self, select_func, *kwargs):
    #     to_return = []
    #     for res in self.residues:
    #         if select_func(res, *kwargs):     
    #             to_return.append(res)
    #     return to_return   

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
        """
        Method: select_residues_from_name

        Select specific residues from the molecular system based on their name and filter criteria.

        Parameters:
        - name (str): Name of the residues to be filtered and selected.
        - select_func (function): A filtering function that takes a Residue object and optional
                                arguments (kwargs) as input and returns a boolean value.
        - *kwargs: Optional additional arguments to be passed to the select_func.

        Returns:
        - list: A list of Residue objects that meet the filtering criteria.

        This method allows selection of specific Residue objects from the molecular system based on
        their shared name (residue name). The input parameter 'name' corresponds to the name of the
        residues to be selected. The method also requires a filtering function 'select_func', which
        takes a Residue object as input and optional keyword arguments (kwargs). The filtering function
        should return True for the Residue objects that match the desired criteria and False for others.

        The method iterates through all Residue objects with the specified name and applies the filtering
        function to each of them. Residues that meet the criteria are included in the returned list,
        providing access to their properties and atoms.

        Example Usage:
        Suppose you want to select all residues with the name "ALA" whose atoms have a maximum value of 10.0 in the X-coordinate. 
        You can define a filtering function like this:
        def custom_filter(residue, max_x_coord):
            for atom in residue.atoms:
                if atom.x > max_x_coord:
                    return False
            return True

        Then, you can call the method like this:
        selected_residues = system.select_residues_from_name("ALA", custom_filter, max_x_coord=10)

        Raises:
        - ResidueNotFound: If no residues with the specified name are found in the system.
        """

        residues_to_look_at = self.get_residue_by_name(name)
        to_return = []
        for res in residues_to_look_at:
            if select_func(res, *kwargs):
                to_return.append(res)
        return to_return

    def write_gro(self, gro_path):
        """
        Method: write_gro

        Write the molecular system's coordinates to a GROMACS .gro file.

        Parameters:
        - gro_path (str): Path to the output GROMACS .gro file.

        Example Usage:
        Suppose you have a 'system' instance of the GroSystem class representing the molecular system
        with residues and atoms. You can write the system's coordinates to a GROMACS .gro file like this:
        system.write_gro("output.gro")
        """
        with open(gro_path, 'w') as o:
            o.write(f'{self.name}\n') # Write the system name as the first line in the .gro file
            o.write(f' {len(self.atoms)}\n') # Write the total number of atoms as the second line
            for res in self.residues: # Iterate through all residues in the system
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
        print(f'System gro file written in {gro_path}')
    
    def write_top(self, top_path):
        """
        Method: write_top

        Write the molecular system's topology to a GROMACS .top file.

        Parameters:
        - top_path (str): Path to the output GROMACS .top file.

        This method writes the molecular system's topology to a GROMACS .top file. The input
        parameter 'top_path' specifies the file path where the topology data will be saved.

        This method write a simple topology file including the following informations : the #includes lines readed from the top file you provide to construct the system,the system name and the
        number of molecules of each type (residue).

        Example Usage:
        Suppose you have a 'system' instance of the GroSystem class representing the molecular system
        with multiple residue types. You can write the system's topology to a GROMACS .top file like this:
        system.write_top("output.top")
        """

        with open(top_path, 'w') as o:
            for include in self._top_includes:
                o.write(f'{include}\n')
            o.write('\n')
            o.write(f'[ system ]\n; name\n{self.name}\n\n')
            o.write(f'[ molecules ]\n; name number\n')
            for residue_stack in self._residues_stack:
                if residue_stack[0].name == 'ION' and self.is_martini:
                    print('Your system is martini with ion, ION will be replaced by atom name on the top file')
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
        print(f'System top file written in {top_path}') 

class Residue:
    def __init__(self, number: int, name: str, idx, system):
        """
        Class: Residue

        Represents a molecular residue in a molecular system.

        Parameters:
        - number (int): The unique number assigned to the residue.
        - name (str): The name of the residue (e.g., "ALA", "GLY").
        - idx: An index or position of the residue in the molecular system.
        - system (GroSystem): A reference to the GroSystem containing the residue.

        Attributes:
        - number (int): The unique number assigned to the residue.
        - name (str): The name of the residue (e.g., "ALA", "GLY").
        - idx: An index or position of the residue in the molecular system.
        - system (GroSystem): A reference to the GroSystem containing the residue.
        - atoms (list): A list of Atom objects representing atoms in the residue.

        Example Usage:
        Creating a Residue object:
        residue = Residue(number=1, name="ALA", idx=0, system=my_system)
        """
        
        self.system = system
        self.name = name
        self.number = number
        self.idx = idx
        self.atoms = []

    def __repr__(self):
        """
        Method: __repr__

        Get the string representation of the Residue object.

        Returns:
        - str: A string representation of the Residue object.

        Example Usage:
        residue = Residue(number=1, name="ALA", idx=0, system=my_system)
        print(residue)  # Output: "ALA1 with 10 atoms"
        """
        return f"{self.name}{self.number} with {len(self.atoms)} atoms"
    
    def add_atom(self, name, number, coordinates, velocities):
        """
        Method: add_atom

        Add an Atom object to the residue.

        Parameters:
        - name (str): The name of the atom (e.g., "CA", "CB").
        - number (int): The unique number assigned to the atom.
        - coordinates (list): The X, Y, Z coordinates of the atom.
        - velocities (list): The velocities of the atom (optional).

        Example Usage:
        residue = Residue(number=1, name="ALA", idx=0, system=my_system)
        residue.add_atom(name="CA", number=1, coordinates=[0.0, 1.0, 2.0], velocities=[0.1, 0.2, 0.3])
        """
        self.atoms.append(Atom(name, number, coordinates, velocities, self))

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
    
    def copy(self):
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
        print(f'I want to copy <{self}>')
        residues_stack = self.system.get_residue_by_name(self.name)
        print(f'Will be copied at the end of {len(residues_stack)} {self.name} stack')
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
            new_coords[0] = round(new_coords[0] + 0.21, 3) # Adjust X-coordinate by adding 0.21 (vdW radius).
            # Create a new Atom object and add it to the new residue's 'atoms' list.
            new_residue.add_atom(atom.name, new_atom_number, new_coords, atom.velocities)
            new_atom_number = add_to_gro_number(new_atom_number, 1)
        
        # Insert the new Residue into the molecular system using the 'insert_residue' method.
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
    
    @property
    def x(self):
        return self.coordinates[0]
    
    @property
    def y(self):
        return self.coordinates[1]
    
    @property
    def z(self):
        return self.coordinates[2]
    
def add_to_gro_number(number, to_add):
    """
    Function: add_to_gro_number

    Add a value to a GROMACS-style atom or residue number while handling rollover.

    Parameters:
    - number (int): The current atom or residue number.
    - to_add (int): The value to add to the current number.

    Returns:
    - int: The new atom or residue number after addition, considering rollover.

    This function is used to add a specified value to a GROMACS-style atom or residue number, while taking
    into account the rollover behavior. In GROMACS, atom and residue numbers are typically limited to five digits
    (ranging from 0 to 99999). When adding a value to the current number, if the result exceeds 99999, the
    number rolls over to 0 and continues counting. This function ensures that the rollover behavior is properly
    handled and the new number stays within the valid range.

    Example Usage:
    current_number = 99998
    added_value = 5
    new_number = add_to_gro_number(current_number, added_value)
    # In this case, new_number would be 3, as it rolls over from 99998 to 99999 and then to 1.

    Notes:
    This function assumes that the maximum atom or residue number in GROMACS is 99999.
    """
    new_number = number + to_add
    if new_number >= 99999 :
        new_number = number - 99999
    return new_number

class GroParsingException(Exception):
    pass

class ResidueNotFound(Exception):
    def __init__(self, residue):
        message = f'Residue {residue} not found'
        super().__init__(message)