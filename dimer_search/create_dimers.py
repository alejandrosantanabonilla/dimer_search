from rdkit import Chem
import numpy as np
from utils_mol import MolFrag, CreateDimer

from ase import Atoms
from ase.io import write
from ase import Atoms
from ase.optimize import BFGS
from ase.optimize.minimahopping import MinimaHopping
from tblite.ase import TBLite
from ase.units import Bohr,Rydberg,kJ,kB,fs,Hartree,mol,kcal
from ase.io import read, write
from ase.optimize.minimahopping import MHPlot

class DimerProcessor:
    def __init__(self, pdb_filename, atom_indices, output_filename_prefix="dimer"):
        self.pdb_filename = pdb_filename
        self.atom_indices = atom_indices
        self.output_filename_prefix = output_filename_prefix
        self.mol = Chem.MolFromPDBFile(self.pdb_filename, removeHs=False)
        self.mol_frag = MolFrag(self.mol)
        self.mol1, self.mol2 = self._get_fragments()
        self.dimer = CreateDimer(self.mol1, self.mol2)
        self.joined_molecule = None

    def _get_fragments(self):
        pdb_blocks = self.mol_frag.get_mol_frags_blocks()
        return Chem.MolFromPDBBlock(pdb_blocks[0], removeHs=False), Chem.MolFromPDBBlock(pdb_blocks[1], removeHs=False)

    def rdkit_mol_to_ase_atoms(self, rdkit_mol: Chem.Mol) -> Atoms:
        """Convert an RDKit molecule to an ASE Atoms object."""
        ase_atoms = Atoms(
            numbers=[atom.GetAtomicNum() for atom in rdkit_mol.GetAtoms()],
            positions=rdkit_mol.GetConformer().GetPositions()
        )
        return ase_atoms

    def rotate_and_join(self, z_angle):
        """Rotates the second fragment, joins, and optionally writes to file."""
        rotated_mol2 = self.dimer.rotate_molecule_around_plane(
            self.dimer.mol2,
            atom_indices=self.atom_indices,
            z_angle=z_angle
        )

        ase_mol1 = self.rdkit_mol_to_ase_atoms(self.mol1)
        ase_rotated_mol2 = self.rdkit_mol_to_ase_atoms(rotated_mol2)
        self.joined_molecule = ase_mol1 + ase_rotated_mol2

        output_filename = f"{self.output_filename_prefix}_initial.xyz"
        write(output_filename, self.joined_molecule)
        return self.joined_molecule

    def relax(self, method="GFN2-xTB", Ediff0=1.0, T0=2000.0, totalsteps=10):
        """Relaxes the joined molecule using TBLite and MinimaHopping with options.

        Args:
            method (str): The TBLite method to use (default: "GFN2-xTB").
            Ediff0 (float): The energy difference threshold for MinimaHopping (default: 1.0).
            T0 (float): The initial temperature for MinimaHopping (default: 2000.0).
            totalsteps (int): The total number of MinimaHopping steps (default: 10).

        Returns:
            ase.Atoms: The relaxed ASE Atoms object.
        """
        if self.joined_molecule is None:
            raise ValueError("Molecule must be joined first. Call rotate_and_join().")

        print(f"Relaxing molecule with: method={method}, Ediff0={Ediff0}, T0={T0}, totalsteps={totalsteps}")  # Print parameters

        calculator = TBLite(method=method)
        self.joined_molecule.set_calculator(calculator)

        hop = MinimaHopping(atoms=self.joined_molecule, Ediff0=Ediff0, T0=T0)
        hop(totalsteps=totalsteps)

        mhplot = MHPlot()
        mhplot.save_figure(f"{self.output_filename_prefix}_summary.png")

        traj_filename = "minima.traj"
        traj = read(traj_filename)
        write(f"{self.output_filename_prefix}_minima.xyz", traj, format="xyz")

        return self.joined_molecule


