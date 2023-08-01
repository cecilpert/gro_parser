# gro_parser

The gro_parser module is a Python implementation that allows reading and writing GROMACS coordinate (gro) and topology (top) files. It provides a set of classes to represent a molecular system, including residues and atoms, and supports basic operations such as adding and copying residues.

## Installation

```bash
git clone https://github.com/cecilpert/gro_parser.git
pip install -e gro_parser
```

## Usage

### Loading a system and accessing informations
```python
from gro_parser import load_system

gro_file_path = 'path/to/your.gro'
top_file_path = 'path/to/your.top'

system = load_system(gro_file_path, top_file_path)

# Browse all residues 
for residue in system.residues:
    print(residue)

# Get a specific residue by index or residue number 
residue = system.get_residue_by_idx(0) 
residue = system.get_residue_by_number(1)

# Get all residues with a specific name
residues_with_name = system.get_residue_by_name('RESIDUENAME')
```

### Manipulating the system 
```python
# Copy a residue. Example : I want to copy the first DOPC of my system. It will insert it at the end of DOPC stack
residue_to_copy = system.get_residue_by_name('DOPC')[0]
residue_to_copy.copy()
```

### Writing the system to gro and top files
```python
# Write the gro file corresponding to the system
system.write_gro('path/to/output.gro')

# Write the top file corresponding to the system
system.write_top('path/to/output.top')
```
