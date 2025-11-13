import numpy as np
from dimer_search import *
import glob
import os
from ase.build import molecule
from ase.io import write, read
from ase import Atoms


print("Running Horizontal Dimer Example...")

# --- Define Orientation and Position ---
# To arrange in the "same direction", all rotation angles for the second molecule should be 0.
# Angles for [Molecule 1, Molecule 2]
yaw = np.deg2rad([0, 50])      # No rotation around Z-axis
pitch = np.deg2rad([180, 0])   # No rotation around Y-axis
roll = np.deg2rad([-180, 0])    # No rotation around X-axis

# To arrange "horizontally", we place one at the origin and translate the other along an axis (e.g., X).
# Positions for [Molecule 1, Molecule 2]
translations = [
    [0.0, 0.0, 0.0],  # Molecule 1 at the origin
    [0.0, 0.0, -14.5]   # Molecule 2 translated 5 Angstroms along X-axis
]

# --- Configs (kept from original example, though not used by mock) ---
tblite_config = {
    "electronic_temperature": 5500.0,
    "max_iterations": 300,
}

mh_config = {
    "T0": 1200.0,
    "Ediff0": 0.6,
    "fmax": 0.1
}

# --- Define Input File ---
# This script assumes a file named 'mol.xyz' already exists in the same directory.
input_filename = "mol.xyz"

if not os.path.exists(input_filename):
    print(f"Warning: Input file '{input_filename}' not found.")
    print("Please create a 'mol.xyz' file. The script will likely fail.")
else:
    print(f"Using existing '{input_filename}' as the base molecule.")


# --- Process the Dimer ---
print("\n--- Creating Horizontal Dimer (No Relaxation) ---")
output_filename = "horizontal_dimer.xyz"
processor = MoleculeProcessor(input_file=input_filename, output_file=output_filename)

# Call process_molecules with relax_molecule=False to just get the geometry
result = processor.process_molecules(
        yaw_rad=yaw, pitch_rad=pitch, roll_rad=roll, translation_vector=translations,
        relax_molecule=True,                # Turn relaxation ON
        restart_relax=False,                # Start fresh
        totalsteps=45,                      # Short relaxation run
        tblite_params=tblite_config,        # Pass tblite params  
        mh_params= mh_config,               # Pass MH params
        output_filename_prefix="relax_run"  # Prefix for relax outputs
    )

if result:
    print(f"\nTest Case completed. Check '{output_filename}'.")
    print("The resulting file 'horizontal_dimer.xyz' should contain two molecules,")
    print("both with the same orientation, separated by 5 Angstroms along the X-axis.")
else:
    print("\nDimer processing failed.")
