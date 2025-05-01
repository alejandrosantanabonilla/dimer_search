# Dimer Search via Basin Hopping with RDKit, ASE, and TBLite

This repository provides Python code for processing molecule dimers, including rotation, joining, and relaxation
using RDKit, ASE, and TBLite.

## Code Structure

This Python code defines a MoleculeProcessor class using the ASE library to:

1. Load a molecule from an input file (e.g., XYZ format).
2. Align the molecule by moving its center of mass to the origin and rotating its principal axes to align with the Cartesian axes.
3. Generate an assembly by creating multiple copies of the aligned molecule, applying specified rotations (yaw, pitch, roll) and translations to each copy.
4. Write the initial, unrelaxed assembly to an output file.
5. Optionally relax the final assembly using an external relax function (imported from utils.py), with support for restarting the relaxation from a previous state.

The main process_molecules method orchestrates these steps.

## Installation

It's highly recommended to use a virtual environment to manage dependencies.  Here's how:

1.  **Create a virtual environment:**

```bash
python3 -m venv dimer_env  # Creates a virtual environment named "dimer_env"
```

2. **Activate the virtual environment:**


```bash
source dimer_env/bin/activate  # On Linux/macOS
```

3. **Clone the repository:**

```bash
git clone https://github.com/alejandrosantanabonilla/dimer_search.git
cd dimer_search
```

4. **Install dimer_search**

```bash
pip install .
```

# Tutorial: Using MoleculeProcessor for Generating and Relaxing Molecular Assemblies

This tutorial explains how to use the `MoleculeProcessor` class, presumably from a library like `dimer_search`, to generate molecular assemblies based on specified orientations and translations. It also covers how to perform structural relaxation (e.g., energy minimization or conformational search) on these assemblies and utilize the restart feature for long calculations. We will walk through the three distinct test cases provided in the example script.

## Prerequisites

* A working Python environment.
* Necessary libraries installed: `numpy` and the library containing `MoleculeProcessor` (referred to as `dimer_search` in the script).
* An input coordinate file named `mol.xyz` located in the same directory as the script. This file should contain the atomic coordinates of the single molecule you wish to assemble into dimers or larger structures.

## Initial Setup

The script begins by importing required libraries and defining the parameters that control how the molecular assemblies are generated and, optionally, relaxed.

```python
import numpy as np
from dimer_search import *
import glob
import os
from ase.build import molecule
from ase.io import write, read
from ase import Atoms

print("Running MoleculeProcessor Example...")

# --- Create the initial system (Two Water Molecules) ---
# Create the first water molecule
water1 = molecule('H2O')

# Create a second water molecule and translate it slightly
water2 = molecule('H2O')

water2.translate([0.0, 0.0, 5.0]) # Example translation

# Combine them into a single Atoms object
initial_system = water1 + water2 # ASE allows direct addition

# Save the initial system to mol.xyz so the processor can read it
input_filename = "mol.xyz"
write(input_filename, initial_system)
print(f"Created initial system with {len(initial_system)} atoms and saved to '{input_filename}'.")

print("Running MoleculeProcessor Example...")

# --- Define Orientations ---
# Yaw, Pitch, and Roll angles specify the rotation applied to the second molecule
# relative to the first one. Multiple angles can be provided to generate
# different relative orientations. Angles are converted to radians.
yaw = np.deg2rad([0, 180])
pitch = np.deg2rad([0, 180])
roll = np.deg2rad([0, 0])

# --- Define Translations ---
# These vectors specify the displacement (shift) of the second molecule's
# center relative to the first molecule's center after rotation.
# Each translation vector corresponds to a set of (yaw, pitch, roll) angles.
translations = [
    [10.5, 0.0, 8.0], # Corresponds to yaw[0], pitch[0], roll[0]
    [1.0, 0.0, 9.0]   # Corresponds to yaw[1], pitch[1], roll[1]
]

```

This initial block sets up two distinct relative placements (defined by rotation and translation) for creating dimer structures from the base molecule in mol.xyz. It also prepares
configuration dictionaries (tblite_config, mh_config) that will only be used if the relaxation feature is turned on in the subsequent steps.

## Test Case 1: Generate Assembly Without Relaxation

### Goal

To generate the initial dimer structures based on the defined orientations and translations without performing any subsequent energy minimization or structural relaxation. This is 
useful for quickly creating and visualizing the starting geometries.

```python
# --- Test Case 1: Process without relaxation ---
print("\n--- Test Case 1: No Relaxation ---")

# Initialize the processor with input and the desired output file for the assemblies
processor1 = MoleculeProcessor(input_file="mol.xyz", output_file="initial_assembly.xyz")

# Process the molecules
result1 = processor1.process_molecules(
    yaw_rad=yaw,
    pitch_rad=pitch,
    roll_rad=roll,
    translation_vector=translations,
    relax_molecule=False # Key parameter: Relaxation is turned OFF
)

# Check if the process completed successfully
if result1:
    print("Test Case 1 completed. Check 'initial_assembly.xyz'.")
```

### Explanation
1. An instance of MoleculeProcessor is created. It's configured to read from mol.xyz and write the generated structures to initial_assembly.xyz.
2. The process_molecules method is called with the orientation angles (yaw_rad, pitch_rad, roll_rad) and translation_vector defined earlier.
3. The crucial parameter here is relax_molecule=False. This instructs the processor to:
    - Read the base molecule from mol.xyz.
    - For each specified placement (combination of rotation and translation):
      - Create a copy of the base molecule.
      - Apply the rotation (yaw, pitch, roll) and translation to the copy.
      - Combine the original molecule and the transformed copy to form a dimer.
4. Write the coordinates of all generated dimer structures sequentially into the specified output_file (initial_assembly.xyz).
   No energy calculations or geometry changes are performed.

### Expected Output
A single file named initial_assembly.xyz. This file will contain the atomic coordinates of the generated dimer structures, concatenated one after another in standard XYZ format. 
The number of structures will match the number of orientation/translation pairs provided (two in this example).

## Test Case 2: Generate Assembly WITH Relaxation (Fresh Start)

### Goal
To generate the initial dimer structures (as in Test Case 1) and then perform a structural relaxation or energy minimization on each generated structure. This test case starts the relaxation process from scratch, ignoring any previous relaxation runs.

```python
# --- Test Case 2: Process WITH relaxation (fresh start) ---
print("\n--- Test Case 2: With Relaxation (Fresh Start) ---")

# Clean up potential relaxation output files from any previous runs
# This ensures that restart_relax=False truly starts fresh.
for f in ['minima.traj', 'qn*.traj', 'hop.log', 'relax_run_*.xyz', 'relax_run_*.png']:
    for ff in glob.glob(f):
        if os.path.exists(ff): os.remove(ff)

# Initialize the processor. Note: output_file still refers to the *initial* assembly.
processor2 = MoleculeProcessor(input_file="mol.xyz", output_file="initial_assembly_for_relax.xyz")

# Process the molecules, enabling relaxation
result2 = processor2.process_molecules(
    yaw_rad=yaw, pitch_rad=pitch, roll_rad=roll, translation_vector=translations,
    relax_molecule=True,          # Key parameter: Relaxation is turned ON
    restart_relax=False,         # Key parameter: Do NOT restart, start fresh
    totalsteps=50,               # Number of steps for the relaxation simulation/optimization
    tblite_params=tblite_config, # Pass the calculator configuration
    mh_params= mh_config,        # Pass the relaxation algorithm configuration
    output_filename_prefix="relax_run" # Prefix for relaxation-specific output files
)

if result2:
    print("Test Case 2 completed. Check 'initial_assembly_for_relax.xyz', 'relax_run_minima.xyz', etc.")
```

### Explanation

1. Cleanup: The code first attempts to delete common output files associated with the relaxation process (like .traj, .log files). This step is crucial to guarantee that setting restart_relax=False results in a completely new relaxation run, unaffected by any previous attempts.
2. An instance of MoleculeProcessor is created. The output_file (initial_assembly_for_relax.xyz) will contain the initial, unrelaxed structures generated before relaxation begins.
3. The process_molecules method is called with several important parameters enabled for relaxation:
   - relax_molecule=True: This activates the relaxation procedure for each dimer generated.
   - restart_relax=False: This explicitly tells the relaxation algorithm not to look for or use any existing state/restart files. Each relaxation starts from the initial generated geometry.
   - totalsteps=50: This defines the duration or extent of each relaxation run (e.g., maximum number of optimization steps, molecular dynamics steps, or search steps). Here, it's set to a small value for a short demonstration run.
   - tblite_params, mh_params: The configuration dictionaries defined earlier are passed to configure the energy calculator and the relaxation algorithm, respectively.
   - output_filename_prefix="relax_run": This string is used as a base name for various output files generated during the relaxation phase. This helps organize outputs when multiple relaxation runs might occur.

### Expected Output

1. initial_assembly_for_relax.xyz: Contains the coordinates of the initial, unrelaxed dimer structures (should be identical in content to initial_assembly.xyz from Test Case 1).
2. Files starting with the prefix relax_run: These files contain the results and trajectory of the relaxation process itself. The exact files depend on the dimer_search library's implementation, but common examples include:
   - relax_run_minima.xyz: Often contains the lowest-energy structures found during the relaxation for each starting dimer.
   - minima.traj or qn*.traj: Trajectory files storing intermediate geometries or identified local minima. These are often crucial for the restart functionality.
   - hop.log: A log file detailing the steps, energies, acceptance criteria, etc., of the relaxation/search algorithm.
   - Other potential files like images (.png) or intermediate structure files (.xyz).

## Test Case 3: Generate Assembly WITH Relaxation (Attempting Restart)

### Goal

To generate the initial assemblies and perform relaxation, similar to Test Case 2, but this time attempting to continue or restart the relaxation process from where a previous run (specifically, Test Case 2 in this sequence) left off. This is extremely useful for extending long simulations or recovering from interruptions 
without starting over.

```python
# --- Test Case 3: Process WITH relaxation (attempt restart) ---
# This test assumes Test Case 2 ran successfully and created the necessary
# restart files (e.g., minima.traj, qn*.traj, hop.log with the 'relax_run' prefix).
print("\n--- Test Case 3: With Relaxation (Attempt Restart) ---")

# Initialize the processor
processor3 = MoleculeProcessor(input_file="mol.xyz", output_file="initial_assembly_for_restart.xyz")

# Process molecules, attempting to restart the relaxation
result3 = processor3.process_molecules(
    yaw_rad=yaw, pitch_rad=pitch, roll_rad=roll, translation_vector=translations,
    relax_molecule=True,          # Key parameter: Relaxation is turned ON
    restart_relax=True,          # Key parameter: Signal to attempt RESTART
    totalsteps=240,              # Set the *new* target total number of steps
    tblite_params=tblite_config, # Pass calculator configuration
    mh_params= mh_config,        # Pass relaxation algorithm configuration
    output_filename_prefix="relax_run" # MUST match the prefix from the run to be restarted
)

if result3:
    print("Test Case 3 completed (Restart likely occurred).")
else:
    print("Test Case 3 failed (Perhaps expected if restart files were missing or previous run failed).")

# Clean up dummy input file (if created by the library)
if os.path.exists("dummy_input.xyz"): os.remove("dummy_input.xyz")
print("\nExample execution finished.")
```

### Explanation

1. Dependency: This test case explicitly relies on Test Case 2 having been executed successfully beforehand. It needs the output files generated by Test Case 2 (specifically those associated with the relax_run prefix) to perform the restart.
2. An instance of MoleculeProcessor is created. Again, output_file (initial_assembly_for_restart.xyz) will just contain the initial structures.
3. The process_molecules method is called with parameters configured for a restarted relaxation:
    - relax_molecule=True: Relaxation is still enabled.
    - restart_relax=True: This is the key flag that activates the restart mechanism. The processor and its underlying relaxation algorithm will now look for existing state files.
    - totalsteps=240: This specifies the target total number of steps for the simulation. If Test Case 2 ran for 50 steps and the restart is successful, this run will perform approximately 190 additional steps (240 total - 50 previous = 190 new). If restart fails and it starts fresh, it will run for 240 steps.
    - output_filename_prefix="relax_run": This must be identical to the prefix used in the run you intend to continue (Test Case 2). This tells the processor which set of output files to look for and update.

### Expected Output

1. initial_assembly_for_restart.xyz: Contains the initial, unrelaxed dimer structures.
2. Updated files starting with relax_run: If the restart was successful, the existing files (relax_run_minima.xyz, minima.traj, qn*.traj, hop.log, etc.) will be appended to or updated. The simulation history within these files will now reflect a total run length approaching the new totalsteps value (240). If the restart failed (e.g., files were missing), these files might be overwritten with a new run of 240 steps, or the process might error out depending on the library's implementation.
3. The console output provides a hint about whether the restart likely succeeded based on the completion status.

## Detailed Explanation of the Restart Option (restart_relax=True) 

The restart_relax=True parameter enables the continuation of a potentially long and computationally expensive relaxation or conformational search process.

### How it Works:

1. Signal to Resume: When process_molecules is called with restart_relax=True, it signals to the internal relaxation algorithm that it should attempt to resume a previous calculation rather than starting anew.
2. File Identification: The algorithm uses the output_filename_prefix (e.g., "relax_run") to identify the relevant set of files from a previous run in the current working directory.
3. State File Check: It searches for specific files known to store the state of the simulation. Common candidates include:
   - Trajectory Files (.traj): Files like minima.traj or qn_*.traj often store the sequence of atomic coordinates visited, potentially including the very last configuration reached.
   - Log Files (.log): Files like hop.log usually record the progress, including the number of steps already performed, energies, acceptance rates, and possibly other algorithm-specific state information.
   - Dedicated Restart Files: Some algorithms might write a specific checkpoint or restart file containing a more complete snapshot of the simulation state (coordinates, velocities, optimizer state, random number generator state, step count, etc.).
4. State Restoration: If valid and compatible restart files are found, the algorithm reads the necessary information to restore its internal state as closely as possible to where the previous run left off. This typically involves loading the last coordinates and retrieving the number of steps already completed.
5. Continuation: The relaxation process then resumes from this restored state. It continues executing steps until the cumulative number of steps reaches the new totalsteps value specified in the current process_molecules call.

### Why Use Restart?

1. Long Simulations: Break down simulations that might take days or weeks into manageable segments. Run for 12 hours, save state, restart.
2. Resource Limits: Overcome wall-time limits on computing clusters. Run until the time limit, then submit a new job that restarts from the previous state.
3. Interruption Recovery: Recover from unexpected interruptions like power outages, system crashes, or accidental process termination without losing all progress.
4. Progressive Exploration: Run an initial short exploration, analyze the results (e.g., check energies in relax_run_minima.xyz), and then decide to continue the simulation for longer if needed.

### Requirements for Successful Restart:

1. A previous run using relax_molecule=True must have completed at least partially and successfully generated the necessary output/state files.
2. The output_filename_prefix used in the restart call must exactly match the prefix used in the run you want to continue.
3. The required state/output files from the previous run must be present in the directory where the restart script is executed.
4. The relax_molecule=True flag must also be set in the restart call.
5. Ideally, other simulation parameters (tblite_params, mh_params) should remain consistent, although some libraries might tolerate minor changes.

## Running the Example

1. Save the Python code provided into a file (e.g., run_dimer_tests.py).
2. Create or obtain a valid molecular structure file named mol.xyz and place it in the same directory.
3. Ensure your Python environment has numpy and the dimer_search library installed.
4. Open a terminal or command prompt, navigate to the directory containing the script and mol.xyz.
5. Execute the script: python run_dimer_tests.py


Monitor the console output, which will indicate the progress through the three test cases. After execution, inspect the generated files (initial_assembly*.xyz, relax_run*.*) to see the results of each step. Compare 
initial_assembly.xyz and initial_assembly_for_relax.xyz (they should be the same). Examine the relax_run files after Test Case 2 and see how they are potentially updated after Test Case 3.
This tutorial provides a comprehensive guide to using the MoleculeProcessor for generating molecular assemblies, performing




