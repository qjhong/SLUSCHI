from ase.io import read, write
from ase.build import make_supercell

# Read the VASP POSCAR file
atoms = read('POSCAR')
# Define a transformation matrix for the supercell (e.g., 2x2x2)
supercell_matrix = [[2, 0, 0],
                    [0, 2, 0],
                    [0, 0, 2]]

# Build the supercell
supercell = make_supercell(atoms, supercell_matrix)
# Write the supercell to a LAMMPS data file
write('lammps.data', supercell, format='lammps-data')
