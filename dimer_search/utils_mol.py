import numpy as np
from scipy.spatial.transform import Rotation as R

from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import AllChem, rdMolAlign

import parmed as pmd
from io import StringIO

class CreateDimer:
    """
    Class for creating and manipulating molecule dimers.
    This class focuses on applying transformations to the second molecule (mol2) 
    while keeping the first molecule (mol1) unchanged.
    """
    def __init__(self, mol1, mol2):
        """
        Initializes CreateDimer object with two RDKit molecules.

        Args:
            mol1: The first RDKit molecule object.
            mol2: The second RDKit molecule object.
        """
        self.mol1 = mol1
        self.mol2 = mol2

    def rotate_molecule_around_plane(self, mol, atom_indices, x_angle=0, y_angle=0, z_angle=0):
        """
        Rotates the given molecule around a plane defined by the provided atom indices.
        The plane's normal vector is calculated using SVD.
        Additional rotations around X, Y, and Z axes can also be applied using scipy.spatial.transform.Rotation.
        The function does not modify the original molecule and returns a new molecule with the rotated conformation.

        Args:
            mol: The RDKit molecule object to rotate.
            atom_indices: A list of atom indices (0-based) defining the plane.
            x_angle: Rotation angle around the X-axis in degrees.
            y_angle: Rotation angle around the Y-axis in degrees.
            z_angle: Rotation angle around the Z-axis in degrees.

        Returns:
            A new RDKit molecule object with the rotated conformation.
        """

        # Create a copy of the molecule to avoid modifying the original
        mol = Chem.Mol(mol)

        # Get coordinates of the molecule
        conf = mol.GetConformer()
        coords = conf.GetPositions()

        # Select atoms defining the plane
        plane_atoms = coords[atom_indices]

        # Compute the best-fit plane using SVD
        centroid = np.mean(coords, axis=0)  # Use the centroid of the entire molecule
        plane_atoms -= centroid
        _, _, vh = np.linalg.svd(plane_atoms)
        normal_vector = vh[-1, :]  # Normal vector of the plane

        # Create a rotation object using scipy.spatial.transform.Rotation
        # Correct the rotation axis by inverting the normal vector
        #rotation = R.from_rotvec(-np.pi * normal_vector)  # Rotate by 180 degrees around the inverted normal vector

        # Apply additional rotations around X, Y, and Z axes
        rotation = R.from_euler('x', x_angle, degrees=True) * \
                   R.from_euler('y', y_angle, degrees=True) * \
                   R.from_euler('z', z_angle, degrees=True)

        # Apply rotation to the molecule's coordinates
        new_coords = rotation.apply(coords - centroid) + centroid

        # Set the new coordinates for the molecule
        for i in range(len(coords)):
            conf.SetAtomPosition(i, new_coords[i])

        return mol

    
class MolFrag:
    """
    Class for handling molecule fragments.
    """
    def __init__(self, mol, output_format='pdb', write_files=False):
        """
        Initializes MolFrag object with an RDKit molecule.

        Args:
            mol: An RDKit molecule object.
            output_format: The desired output format ('pdb', 'xyz', or 'mol').
            write_files: Whether to write the fragments to files (True/False).
        """
        self.mol = mol
        self.output_format = output_format
        self.write_files = write_files

    def get_mol_frags_blocks(self, pdb_filename_prefix="mol_"):
        """
        Splits a molecule into fragments, generates 3D coordinates,
        and returns a list of PDB blocks for each fragment.

        Args:
            pdb_filename_prefix: Prefix for the PDB filenames.

        Returns:
            A list of PDB blocks (strings) for each fragment.
        """

        mol_frags = Chem.GetMolFrags(self.mol, asMols=True)
        pdb_blocks = []

        # Define a dictionary to map output_format to the corresponding Chem function
        format_to_function = {
            'pdb': Chem.MolToPDBFile,
            'xyz': Chem.MolToXYZFile,
            'mol': Chem.MolToMolFile
        }

        for i, frag in enumerate(mol_frags):
            # Generate PDB block (this remains the same)
            pdb_block = Chem.MolToPDBBlock(frag)  
            pdb_blocks.append(pdb_block)

            if self.write_files:
                # Get the appropriate function from the dictionary
                write_function = format_to_function[self.output_format]  # Directly access using the key

                # Determine file extension
                file_extension = f".{self.output_format}"

                # Call the function to write the file
                write_function(frag, f"{pdb_filename_prefix}{i+1}{file_extension}")

        return pdb_blocks




