from ase.io import read, write

# Read all structures from the LAMMPS dump file
atoms_list = read('lmp.dump', index=':', format='lammps-dump-text')

# Loop through each structure and write it to a POSCAR file
for i, atoms in enumerate(atoms_list):
    write(f'POSCAR_{i}', atoms, format='vasp')
