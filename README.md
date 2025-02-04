# Dimer Search

A Python package for creating and manipulating molecular dimers.  This package provides tools for rotating fragments, joining them, and performing geometry optimizations.

## Installation

```bash
pip install dimer_manipulator
```
## Usage

This is an example for using this code:

```python
  from dimer_manipulator import DimerProcessor  # Import from the package

# 1. Create a DimerProcessor instance, providing the PDB file and atom indices.
#    Make sure "dimer.pdb" is in the same directory as your script, or provide the full path.
processor = DimerProcessor(pdb_filename='dimer.pdb', atom_indices=[0, 1, 2])  # Atom indices for the plane

# 2. Rotate and join the fragments.  This creates the initial joined dimer structure
#    and writes it to "dimer_initial.xyz" (by default, because of how the class is defined).
joined_molecule = processor.rotate_and_join(z_angle=35)  # Rotate by 35 degrees

# 3. Relax the joined molecule using TBLite and MinimaHopping.  This optimizes the
#    geometry and writes the trajectory to "dimer_minima.xyz" and generates a plot
#    "dimer_summary.png".
relaxed_molecule = processor.relax(method="GFN2-xTB", Ediff0=0.5, T0=3000, totalsteps=50)  # Customize the relaxation

# 4. Now you have the relaxed molecule in the `relaxed_molecule` variable.
```
