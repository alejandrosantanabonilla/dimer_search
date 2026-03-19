import os
import copy
import numpy as np
from ase.io import read, write
from tblite.ase import TBLite
from dimer_search import *

print("==================================================")
print("  STEP 1: RIGID-BODY GRID SEARCH PRE-SCREENING    ")
print("==================================================")

# 1. Load the two different molecules
if not os.path.exists('ptb7.xyz') or not os.path.exists('cnnpvn.xyz'):
    raise FileNotFoundError("Make sure 'ptb7.xyz' and 'cnnpvn.xyz' are in the same directory!")

mol1_base = read('ptb7.xyz')
mol2_base = read('cnnpvn.xyz')

# Center them at the origin for predictable math
mol1_base.center()
mol1_base.set_center_of_mass([0, 0, 0])
mol2_base.center()

# 2. Define the Search Grid 
distance = 8.5  # Ångstroms
positions = [
    [distance, 0, 0], [-distance, 0, 0],
    [0, distance, 0], [0, -distance, 0],
    [0, 0, distance], [0, 0, -distance]
]

# 90-degree rotations to check different molecular faces
rotations = [
    (0, 0, 0), (90, 0, 0), (0, 90, 0), 
    (0, 0, 90), (180, 0, 0), (0, 180, 0)
]

results = []
calc = TBLite(method="GFN2-xTB", electronic_temperature=5500.0)
mol1_len = len(mol1_base)

print(f"Testing {len(positions) * len(rotations)} rigid-body configurations...\n")

calc_count = 0
for pos_idx, pos in enumerate(positions):
    for rot_idx, rot in enumerate(rotations):
        
        # Reset geometries for this step
        mol1 = copy.deepcopy(mol1_base)
        mol2 = copy.deepcopy(mol2_base)
        
        # Orient and translate Molecule 2
        mol2.euler_rotate(phi=rot[0], theta=rot[1], psi=rot[2], center='COM')
        mol2.set_center_of_mass(pos)
        
        # Combine
        dimer = mol1 + mol2
        
        # Collision Check (Throw out if atoms are closer than 1.5 Å)
        dist_matrix = dimer.get_all_distances()
        cross_distances = dist_matrix[:mol1_len, mol1_len:]
        
        if np.min(cross_distances) < 1.5:
            print(f"[{calc_count}] Skipped: Collision detected.")
            calc_count += 1
            continue 
            
        # Single Point Energy Calculation
        dimer.calc = calc
        try:
            energy = dimer.get_potential_energy()
            print(f"[{calc_count}] Success: E = {energy:.3f} eV")
            results.append({'energy': energy, 'dimer': dimer})
        except Exception as e:
            print(f"[{calc_count}] Failed SCF calculation. Error: {e}")
            
        calc_count += 1

# 3. Extract the Best Result
if not results:
    raise RuntimeError("All configurations collided! Increase the 'distance' variable.")

results.sort(key=lambda x: x['energy'])
best_dimer = results[0]['dimer']
best_energy = results[0]['energy']

print(f"\nGrid Search Complete! Lowest Energy Found: {best_energy:.3f} eV")

# Save this best configuration to feed into the Processor
best_input_file = "best_grid_start.xyz"
write(best_input_file, best_dimer)
print(f"Saved optimal starting point to '{best_input_file}'.")


print("\n==================================================")
print("  STEP 2: BASIN HOPPING / MINIMA HOPPING RELAX    ")
print("==================================================")

# Since best_grid_start.xyz ALREADY contains the correctly positioned dimer, 
# we tell the processor to apply ZERO further rigid-body movements.
yaw = np.deg2rad([0])    
pitch = np.deg2rad([0])  
roll = np.deg2rad([0])   
translations = [[0.0, 0.0, 0.0]]

# Set up configurations for the Processor
tblite_config = {
    "method": "GFN2-xTB",  
    "electronic_temperature": 500.0,
    "max_iterations": 300,
}

mh_config = {
    "T0": 400.0,
    "Ediff0": 0.2,
    "fmax": 0.1
}

output_filename = "final_relaxed_heterodimer.xyz"
processor = MoleculeProcessor(input_file=best_input_file, output_file=output_filename)

print("Starting molecular processor relaxation...")

# Run the Relaxation
result = processor.process_molecules(
        yaw_rad=yaw, pitch_rad=pitch, roll_rad=roll, 
        translation_vector=translations,
        relax_molecule=True,                # Run the Basin/Minima Hopping
        restart_relax=False,                # Start fresh from the grid search guess
        totalsteps=100,                     # Number of hopping steps
        tblite_params=tblite_config,          
        mh_params=mh_config,               
        output_filename_prefix="relax_run"  
    )

if result:
    print(f"\n==================================================")
    print(f"SUCCESS! The final optimized dimer is saved as '{output_filename}'.")
    print("==================================================")
else:
    print("\nProcessing failed during the relaxation step. Check your logs.")
99999999999999999999999999
