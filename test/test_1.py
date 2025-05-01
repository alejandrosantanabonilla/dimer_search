import numpy as np
from dimer_search import *
import glob
import os
from ase.build import molecule
from ase.io import write, read
from ase import Atoms

print("Running MoleculeProcessor Example...")

yaw = np.deg2rad([0, 180])
pitch = np.deg2rad([0, 180])
roll = np.deg2rad([45, 0])

translations = [
    [5.0, 0.0, 0.0],
    [0.0, 0.0, 0.0]
]

tblite_config = {
    "electronic_temperature": 5500.0,
    "max_iterations": 300,
}

mh_config = {
    "T0": 1200.0,
    "Ediff0": 0.6,
    "fmax": 0.1
}

# --- Create the initial system (Two Water Molecules) ---
# Create the first water molecule
water1 = molecule('H2O')

# Save the initial system to mol.xyz so the processor can read it
input_filename = "mol.xyz"
write(input_filename, water1)
print(f"Created initial system with {len(water1)} atoms and saved to '{input_filename}'.")


# --- Test Case 1: Process without relaxation ---
print("\n--- Test Case 1: No Relaxation ---")
processor1 = MoleculeProcessor(input_file="mol.xyz", output_file="initial_assembly.xyz")
result1 = processor1.process_molecules(yaw_rad=yaw, pitch_rad=pitch, roll_rad=roll, translation_vector=translations,relax_molecule=False)
if result1:
    print("Test Case 1 completed. Check 'initial_assembly.xyz'.")
