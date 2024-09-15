from ase.io import read, write
from ase.atoms import Atoms
import numpy as np

# Step 1: Read the LAMMPS dump file
# Read all structures from the LAMMPS dump file
atoms_list = read('lmp.dump', index=':', format='lammps-dump-text')

# Step 2: Group atoms by element
# Define your custom order of elements (example: Aluminum and Oxygen)
element_order = ['H', 'He', 'Li', 'Be', 'B', 'C', 'N']  # Replace with your desired elements
atoms = atoms_list[0]
element_order = element_order[:len(np.unique(atoms.get_atomic_numbers()))]

# Loop through each structure and write it to a POSCAR file
for i, atoms in enumerate(atoms_list):
    sorted_atoms = atoms.copy()
    sorted_atoms.arrays['positions'] = atoms.arrays['positions'][atoms.get_atomic_numbers().argsort()]
    #write(f'POSCAR_{i}', sorted_atoms, format='vasp')
    # Sort atoms according to the defined element order
    sorted_indices = []
    for element in element_order:
        sorted_indices += [i for i, atom in enumerate(atoms) if atom.symbol == element]
    
    # Create a new Atoms object with atoms in the desired order
    sorted_atoms.arrays['positions'] = atoms.arrays['positions'][atoms.get_atomic_numbers().argsort()]
    
    # Step 3: Write to POSCAR format with the element order explicitly specified
    write('POSCAR', sorted_atoms, format='vasp', direct=True, vasp5=True, sort=False)
    
    # Step 4: Update the POSCAR manually to include element counts if needed
    # Count the number of atoms for each element in the specified order
    element_counts = [sum(1 for atom in sorted_atoms if atom.symbol == element) for element in element_order]
    
    # Modify the POSCAR header to include the correct element types and counts
    with open('POSCAR', 'r') as file:
        lines = file.readlines()
    
    # Manually update the POSCAR header with element types and their counts
    lines[0] = 'mds POSCAR by SLUSCHI\n'  # Line 6: Element names
    lines[5] = ' '.join(element_order) + '\n'  # Line 6: Element names
    lines[6] = ' '.join(map(str, element_counts)) + '\n'  # Line 7: Element counts
    
    # Write the updated POSCAR file
    with open(f'POSCAR_{i}', 'w') as file:
        file.writelines(lines)

