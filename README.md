# Dimer Search

A Python package for creating and manipulating molecular dimers.  This package provides tools for rotating fragments, joining them, and performing geometry optimizations.

## Installation

```bash
pip install dimer_search
```
## Usage

This is an example for using this code:

```python
from dimer_search import DimerProcessor  # Import from the package

# Example usage:
processor = DimerProcessor(pdb_filename='dimer.pdb', atom_indices=[0, 1, 2])
joined_molecule = processor.rotate_and_join(z_angle=35)

relaxed_molecule = processor.relax(method="GFN1-xTB", Ediff0=0.5, T0=3000, totalsteps=5)
print("Relaxation complete. Files written.")

relaxed_molecule = processor.relax()  # Using default values
print("Relaxation complete. Files written.")
```
