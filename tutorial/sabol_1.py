import os
import copy
import numpy as np
from ase.io import read, write
from tblite.ase import TBLite
from scipy.stats import qmc  # Requires: pip install scipy
from ase.geometry import distance
from dimer_search import MoleculeProcessor # Assuming this is in dimer_search.py

# --- Configuration ---
NUM_SAMPLES = 50           # Number of Sobol configurations to test
R_RANGE = [6.5, 7.5, 10.0, 12.0]      # Distance range for COM-to-COM exploration (Angstroms)
COLLISION_THRESHOLD = 2.2  # Minimum allowed distance between any two atoms (Angstroms)
XTB_METHOD = "GFN2-xTB"
E_TEMP = 5000.0            # High electronic temp helps convergence for bad initial guesses

print("==================================================")
print("  STEP 1: 6D SOBOL SAMPLING (RIGID-BODY SEARCH)   ")
print("==================================================")

# 1. Load and prepare molecules
if not os.path.exists('ptb7.xyz') or not os.path.exists('cnnpvn.xyz'):
    raise FileNotFoundError("Check that 'ptb7.xyz' and 'cnnpvn.xyz' exist in the script folder.")

mol1_base = read('ptb7.xyz')
mol2_base = read('cnnpvn.xyz')

# Standardize: Center both at (0,0,0)
mol1_base.center()
mol1_base.set_center_of_mass([0, 0, 0])
mol2_base.center()
mol2_base.set_center_of_mass([0, 0, 0])

mol1_len = len(mol1_base)
calc = TBLite(method=XTB_METHOD, electronic_temperature=E_TEMP)

# 2. Setup 6D Sobol Sampler
# Dim 0: Distance r, Dim 1-2: Position on Sphere, Dim 3-5: Euler Rotations
sampler = qmc.Sobol(d=6, scramble=True)
sample_points = sampler.random(n=NUM_SAMPLES)

results = []

print(f"Sampling {NUM_SAMPLES} configurations in 6D space...")

for i, sample in enumerate(sample_points):
    # Map Sobol [0,1] to physical parameters
    # A. Distance
    r = R_RANGE[0] + sample[0] * (R_RANGE[1] - R_RANGE[0])
    
    # B. Position on the sphere (Uniform surface sampling)
    phi_pos = sample[1] * 2 * np.pi
    cos_theta_pos = 2 * sample[2] - 1 
    sin_theta_pos = np.sqrt(1 - cos_theta_pos**2)
    pos = [
        r * sin_theta_pos * np.cos(phi_pos),
        r * sin_theta_pos * np.sin(phi_pos),
        r * cos_theta_pos
    ]

    # C. Euler Rotations for Mol2
    rot_phi = sample[3] * 360
    rot_theta = sample[4] * 180
    rot_psi = sample[5] * 360

    # Apply transformations to copies
    mol1 = mol1_base.copy()
    mol2 = mol2_base.copy()
    
    mol2.euler_rotate(phi=rot_phi, theta=rot_theta, psi=rot_psi, center='COM')
    mol2.set_center_of_mass(pos)
    
    dimer = mol1 + mol2
    
    # D. Collision Check
    dists = dimer.get_all_distances()
    # Check only distances between atoms of mol1 and atoms of mol2
    cross_dists = dists[:mol1_len, mol1_len:]
    min_d = np.min(cross_dists)
    
    if min_d < COLLISION_THRESHOLD:
        print(f"[{i:03d}] Skipped: Collision ({min_d:.2f} A)")
        continue
    
    # E. Energy Calculation
    dimer.calc = calc
    try:
        energy = dimer.get_potential_energy()
        print(f"[{i:03d}] Success: r={r:.2f} A, E={energy:.4f} eV")
        results.append({'energy': energy, 'dimer': dimer})
    except Exception as e:
        print(f"[{i:03d}] SCF Error: {e}")

# 3. Process Results
if not results:
    raise RuntimeError("Zero valid configurations found. Decrease COLLISION_THRESHOLD or increase R_RANGE.")

results.sort(key=lambda x: x['energy'])
best_dimer = results[0]['dimer']
best_energy = results[0]['energy']

input_for_processor = "best_sobol_guess.xyz"
write(input_for_processor, best_dimer)

print(f"\nPhase 1 Complete. Best Pre-screened Energy: {best_energy:.4f} eV")
print(f"Optimal guess saved to '{input_for_processor}'.")

print("\n==================================================")
print("  STEP 2: LOCAL OPTIMIZATION / BASIN HOPPING      ")
print("==================================================")

# These are placeholders since the dimer is already positioned
yaw = np.deg2rad([0])
pitch = np.deg2rad([0])
roll = np.deg2rad([0])
translations = [[0.0, 0.0, 0.0]]

tblite_config = {
    "method": XTB_METHOD,
    "electronic_temperature": 300.0, # Use physical temp for final relax
    "max_iterations": 500,
}

mh_config = {
    "T0": 400.0,
    "Ediff0": 0.1,
    "fmax": 0.05 # Tighter force convergence
}

output_final = "final_optimized_dimer.xyz"
processor = MoleculeProcessor(input_file=input_for_processor, output_file=output_final)

print("Starting deep relaxation of the best guess...")

success = processor.process_molecules(
    yaw_rad=yaw, pitch_rad=pitch, roll_rad=roll,
    translation_vector=translations,
    relax_molecule=True,
    restart_relax=False,
    totalsteps=100, # More steps for better exploration of the local minimum
    tblite_params=tblite_config,
    mh_params=mh_config,
    output_filename_prefix="step2_relax"
)

if success:
    print(f"\nDONE! Final structure available in '{output_final}'.")
else:
    print("\nRelaxation failed. Check the log files for SCF or force issues.")
