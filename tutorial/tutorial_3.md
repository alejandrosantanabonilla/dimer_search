# Tutorial: Hierarchical Dimer Conformational Search 

Complicated Potential Energy Surfaces due to many different possibilities.

# Installation Guide: `dimer_search` on UCL Young Cluster

This guide provides step-by-step instructions to properly install the `dimer_search` package on the UCL Young cluster. It ensures that the correct compilers, Python versions, and math libraries are loaded before building the package in an isolated virtual environment.

## Step 1: Clean and Load Required Modules
Before creating your environment, purge any conflicting modules and load the specific dependencies required for the quantum chemistry engines to compile correctly.

```bash
# Clear any currently loaded modules to avoid conflicts
module purge

# Load Python 3.11
module load python/3.11.4-gnu-10.2.0

# Load the GNU compilers and OpenBLAS (crucial for underlying math/LAPACK operations)
module load compilers/gnu/10.2.0 
module load openblas/0.3.13-openmp/gnu-10.2.0 
```

## Step 2: Clone the Repository
Navigate to your desired working directory and clone the package from GitHub.

```bash
# Clone the repository
git clone https://github.com/alejandrosantanabonilla/dimer_search.git
```

## Step 3: Create and Activate a Virtual Environment
Create an isolated Python environment so your package dependencies do not interfere with the global cluster modules.

```bash
# Create a virtual environment named 'env_dim'
python3 -m venv env_dim

# Activate the environment
source env_dim/bin/activate
```
*(Note: You must run `source env_dim/bin/activate` every time you log into the cluster to work on this project).*

## Step 4: Install the Package
Move into the cloned directory and install the package in "editable" mode (`-e`). This allows you to modify the source code later without needing to reinstall the package.

```bash
# Navigate into the source code folder
cd dimer_search/

# Install the package and its dependencies
pip install -e .
```

## Step 5: Verify the Installation
Run one of the built-in tutorial scripts to ensure the package and its dependencies (like `tblite`) are communicating correctly.

```bash
# Navigate to the tutorial folder
cd tutorial/

# Run the tutorial script
python3 tutorial_3.py 
```

## Step 6: Optional Clean-Up
If the tutorial generates standard output files that you don't need to keep, you can clean them up using the following command:

```bash
# Remove relaxation trajectories and generated geometries
rm -r md00001.* qn0000* final_relaxed_heterodimer.xyz minima.traj best_grid_start.xyz
```


## 1. The Solution: A Two-Step Hierarchical Workflow
To efficiently find the true global minimum basin, we separate the macro-movements from the micro-relaxations:

1. **Step 1: Rigid-Body Grid Search (Macro-exploration)** We freeze the internal degrees of freedom of both molecules. We anchor Molecule 1 at the origin and 
   systematically translate and rotate Molecule 2 around it. We calculate a fast single-point energy for each configuration, instantly discarding any setups that 
   result in steric clashes (atoms closer than 1.5 Å).
2. **Step 2: Basin Hopping Relaxation (Micro-exploration)**
   We take the absolute lowest-energy configuration from the Grid Search and feed it into our `MoleculeProcessor`. We then apply gentle Minima Hopping to allow 
   the side chains to bend and the $\pi$-$\pi$ stacking distance to optimize, settling into the true local minimum.

---

## 2. Tuning the Grid Search Variables
The success of this method depends on how you define your initial grid. In the script below, you can modify three core components to expand or shrink your search space:

* **`distance`**: The starting center-of-mass separation (in Ångstroms). 
    * *Tip:* Start with `8.0`. If all configurations fail due to collisions, increase it to `10.0` or `12.0`.
* **`positions`**: The coordinates where Molecule 2 will be placed relative to Molecule 1. 
    * *Default:* A 6-point octahedral grid (±X, ±Y, ±Z). 
    * *Advanced:* You can add diagonals like `[distance, distance, 0]` for a denser search.
* **`rotations`**: The Euler angles `(phi, theta, psi)` applied to Molecule 2. 
    * *Default:* 90° increments to test the top, bottom, and side faces.
    * *Advanced:* Add 45° increments or 180° flips depending on your molecule's symmetry.

---

## 3. The Complete End-to-End Script (`dimer_workflow.py`)

Ensure `ptb7.xyz` and `cnnpvn.xyz` are in the same directory as this script.

```python
import os
import copy
import numpy as np
from ase.io import read, write
from tblite.ase import TBLite
from dimer_search import *

print("==================================================")
print("  STEP 1: RIGID-BODY GRID SEARCH PRE-SCREENING    ")
print("==================================================")

# --- 1. Load and Prep Molecules ---
if not os.path.exists('ptb7.xyz') or not os.path.exists('cnnpvn.xyz'):
    raise FileNotFoundError("Missing 'ptb7.xyz' or 'cnnpvn.xyz' in the current directory!")

mol1_base = read('ptb7.xyz')
mol2_base = read('cnnpvn.xyz')

# Center molecules for predictable rotation and translation math
mol1_base.center()
mol1_base.set_center_of_mass([0, 0, 0])
mol2_base.center()

# --- 2. Define the Search Grid ---
distance = 8.0  # Adjust if molecules are too large and constantly collide
positions = [
    [distance, 0, 0], [-distance, 0, 0],
    [0, distance, 0], [0, -distance, 0],
    [0, 0, distance], [0, 0, -distance]
]

# (yaw, pitch, roll) in degrees
rotations = [
    (0, 0, 0), (90, 0, 0), (0, 90, 0), 
    (0, 0, 90), (180, 0, 0), (0, 180, 0)
]

results = []
calc = TBLite(method="GFN2-xTB", electronic_temperature=300.0) # 300K prevents SCF instability
mol1_len = len(mol1_base)

print(f"Testing {len(positions) * len(rotations)} rigid-body configurations...\n")

calc_count = 0
for pos_idx, pos in enumerate(positions):
    for rot_idx, rot in enumerate(rotations):
        
        # Isolate instances for this loop iteration
        mol1 = copy.deepcopy(mol1_base)
        mol2 = copy.deepcopy(mol2_base)
        
        # Apply spatial transformations to Molecule 2
        mol2.euler_rotate(phi=rot[0], theta=rot[1], psi=rot[2], center='COM')
        mol2.set_center_of_mass(pos)
        
        dimer = mol1 + mol2
        
        # --- Collision Filter ---
        dist_matrix = dimer.get_all_distances()
        cross_distances = dist_matrix[:mol1_len, mol1_len:]
        
        # Discard if any atoms from Mol 1 are within 1.5 A of Mol 2
        if np.min(cross_distances) < 1.5:
            print(f"[{calc_count}] Skipped: Steric collision detected.")
            calc_count += 1
            continue 
            
        # --- Single Point Energy Calculation ---
        dimer.calc = calc
        try:
            energy = dimer.get_potential_energy()
            print(f"[{calc_count}] Success: E = {energy:.3f} eV")
            results.append({'energy': energy, 'dimer': dimer})
        except Exception as e:
            print(f"[{calc_count}] Failed SCF calculation. Error: {e}")
            
        calc_count += 1

# --- 3. Extract and Save the Best Result ---
if not results:
    raise RuntimeError("All configurations collided! Increase the 'distance' variable.")

# Sort by lowest energy
results.sort(key=lambda x: x['energy'])
best_dimer = results[0]['dimer']
best_energy = results[0]['energy']

print(f"\nGrid Search Complete! Lowest Energy Found: {best_energy:.3f} eV")

best_input_file = "best_grid_start.xyz"
write(best_input_file, best_dimer)
print(f"Saved optimal starting point to '{best_input_file}'.")


print("\n==================================================")
print("  STEP 2: BASIN HOPPING / MINIMA HOPPING RELAX    ")
print("==================================================")

# The grid search already positioned the molecules perfectly.
# We instruct the processor to apply ZERO further rigid-body movements.
yaw = np.deg2rad([0])    
pitch = np.deg2rad([0])  
roll = np.deg2rad([0])   
translations = [[0.0, 0.0, 0.0]]

# --- Processor Configurations ---
tblite_config = {
    "method": "GFN2-xTB",  
    "electronic_temperature": 300.0, # Keep at room temp to avoid eigenvalue crashes
    "max_iterations": 300,
}

mh_config = {
    "T0": 400.0,     # Keep MD thermal kicks gentle for large organic molecules
    "Ediff0": 0.2,   
    "fmax": 0.1
}

output_filename = "final_relaxed_heterodimer.xyz"
processor = MoleculeProcessor(input_file=best_input_file, output_file=output_filename)

print("Starting molecular processor relaxation...")

# Run the Basin Hopping Relaxation
result = processor.process_molecules(
        yaw_rad=yaw, pitch_rad=pitch, roll_rad=roll, 
        translation_vector=translations,
        relax_molecule=True,                
        restart_relax=False,                
        totalsteps=100,                     
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
```
---
-ubleshooting & Pro-Tips
* **`info=386` or Eigenvalue Errors:** If the relaxation step crashes with a math error, your `T0` (Minima Hopping Temperature) or `electronic_temperature` (TBLite smearing) is too high. Lower `T0` from `400.0` to `300.0` or even `200.0` to prevent the molecules from tearing themselves apart during the MD step.
* **Visualizing the Search:** Once complete, open the `best_grid_start.xyz` and the `final_relaxed_heterodimer.xyz` in a visualizer like **VMD**, **Avogadro**, or **ASE GUI**. You will be able to see exactly how the side-chains and $\pi$-faces rearranged during the hopping process to achieve the final minimum.
