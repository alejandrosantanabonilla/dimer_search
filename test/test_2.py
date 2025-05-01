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

# --- Test Case 2: Process WITH relaxation (fresh start) ---
print("\n--- Test Case 2: With Relaxation (Fresh Start) ---")
# Clean up potential relaxation files from previous runs
for f in ['minima.traj', 'qn*.traj', 'hop.log', 'relax_run_*.xyz', 'relax_run_*.png']:
    for ff in glob.glob(f):
        if os.path.exists(ff): os.remove(ff)

processor2 = MoleculeProcessor(input_file="mol.xyz", output_file="initial_assembly_for_relax.xyz")
result2 = processor2.process_molecules(
        yaw_rad=yaw, pitch_rad=pitch, roll_rad=roll, translation_vector=translations,
        relax_molecule=True,                # Turn relaxation ON
        restart_relax=False,                # Start fresh
        totalsteps=50,                      # Short relaxation run
        tblite_params=tblite_config,        # Pass tblite params  
        mh_params= mh_config,               # Pass MH params
        output_filename_prefix="relax_run"  # Prefix for relax outputs
    )

if result2:
    print("Test Case 2 completed. Check 'initial_assembly_for_relax.xyz', 'relax_run_minima.xyz', etc.")
