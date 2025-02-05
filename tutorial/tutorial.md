# Dimer Processing with RDKit, ASE, and TBLite

This repository provides Python code for processing molecule dimers, including rotation, joining, and relaxation using RDKit, ASE, and TBLite.

## Code Structure

The code consists of three main classes:

*   `MolFrag`: Handles molecule fragmentation.
*   `CreateDimer`: Creates and manipulates dimers (two molecules joined together).
*   `DimerProcessor`: Orchestrates the entire dimer processing workflow, including rotation, joining, and relaxation.

## Installation

It's highly recommended to use a virtual environment to manage dependencies.  Here's how:

1.  **Create a virtual environment:**

```bash
python3 -m venv dimer_env  # Creates a virtual environment named "dimer_env"
```

2. **Activate the virtual environment:**


```bash
source dimer_env/bin/activate  # On Linux/macOS
dimer_env\Scripts\activate  # On Windows
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

## Usage

1. **MolFrag Class**
Purpose: This class is used to break a molecule into fragments.

Initialization:

```python
from rdkit import Chem
mol = Chem.MolFromSmiles('c1ccccc1c2ccccc2') # Example molecule (biphenyl)
mol_frag = MolFrag(mol, output_format='pdb', write_files=True)
```

This creates a MolFrag object.  output_format specifies the format for writing fragment files ('pdb', 'xyz', or 'mol'). write_files=True will write the fragment files to disk.

## get_mol_frags_blocks() function:


```Python
pdb_blocks = mol_frag.get_mol_frags_blocks(pdb_filename_prefix="biphenyl_frag_")
print(pdb_blocks)
```

This function splits the molecule into fragments and returns a list of PDB blocks (strings) for each fragment. It also writes the fragment files to disk according to the specified format and prefix.

2. **CreateDimer Class**

Purpose: This class handles the creation and manipulation of dimers.

Initialization:

```Python
mol1 = Chem.MolFromSmiles('c1ccccc1')  # Benzene
mol2 = Chem.MolFromSmiles('c1ncccc1')  # Pyridine
dimer = CreateDimer(mol1, mol2)
```

This creates a CreateDimer object with two RDKit molecule objects.

## rotate_molecule_around_plane() function:

```Python
rotated_mol2 = dimer.rotate_molecule_around_plane(dimer.mol2, atom_indices=[0, 1, 2], z_angle=30)
```

This rotates mol2 by 30 degrees around the Z-axis, where the Z-axis is defined by atoms with indices 0, 1, and 2.  Note that the atom_indices are 0-based.  You can also specify x_angle and y_angle for rotations around the X and Y axes, respectively.

3. **DimerProcessor Class**

Purpose: This class manages the entire dimer processing workflow.

Initialization:

```Python
processor = DimerProcessor("dimer.pdb", atom_indices=[0, 1, 2], output_filename_prefix="my_dimer")
```

This initializes a DimerProcessor object.  dimer.pdb is the input PDB file containing the molecule to be fragmented. atom_indices specifies the atoms that define the plane of rotation.  output_filename_prefix is used to name output files.

## rotate_and_join() function:

```Python
joined_molecule = processor.rotate_and_join(z_angle=45)
```

This rotates the second fragment by 45 degrees around the Z-axis (defined by the provided atom_indices) and joins the two fragments.  The resulting joined molecule (ASE Atoms object) is written to my_dimer_initial.xyz.

## relax() function:

```Python
relaxed_molecule = processor.relax(method="GFN2-xTB", Ediff0=0.5, T0=1500.0, totalsteps=20)
```

This relaxes the joined molecule using the GFN2-xTB method with specified parameters for MinimaHopping. The relaxed structure is written to my_dimer_minima.xyz.  A summary plot of the MinimaHopping trajectory is saved as my_dimer_summary.png.

## Example Usage

The file **dimer.pdb** can be found in the folder **test**.

```python
from dimer_search import *  

# Example usage:
processor = DimerProcessor(pdb_filename='dimer.pdb', atom_indices=[0, 1, 2])
joined_molecule = processor.rotate_and_join(z_angle=35)

relaxed_molecule = processor.relax(method="GFN1-xTB", Ediff0=0.5, T0=3000, totalsteps=2)
print("Relaxation complete. Files written.")

relaxed_molecule = processor.relax()  # Using default values
print("Relaxation complete. Files written.")
```

This script demonstrates the complete workflow: loading a molecule from a PDB file, rotating one fragment, joining the fragments, and then relaxing the joined molecule.  Remember to install the necessary libraries.
