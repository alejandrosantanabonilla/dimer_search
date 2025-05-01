import ase
from ase.atoms import Atoms
import numpy as np
from ase.io import read, write
from scipy.spatial.transform import Rotation as R

try:
    from .utils import relax
except ImportError:
    # Fallback if running as a script and utils.py is in the same folder
    try:
        from utils import relax
    except ImportError:
        print("ERROR: Cannot find the 'relax' function. Ensure 'utils.py' is in the same directory or Python path.")
        # Define a dummy relax function to allow the script to load
        def relax(*args, **kwargs):
            print("Dummy relax function called. 'utils.py' not found.")
            return kwargs.get('atoms', None) # Return something simple

class MoleculeProcessor:
    """
    A class to process molecular structures, including loading, aligning, rotating, translating,
    and optionally relaxing molecules.

    Attributes:
        input_file (str): Path to the input file containing the molecular structure.
        output_file (str): Path to the output file where the processed structure will be saved.
        original_atoms (Atoms): The original molecular structure loaded from the input file.
        aligned_atoms (Atoms): The molecular structure after alignment and translation.
    """
    
    def __init__(self, input_file, output_file):
        """
        Initializes the MoleculeProcessor with input and output file paths.

        Args:
            input_file (str): Path to the input file containing the molecular structure.
            output_file (str): Path to the output file where the processed structure will be saved.
        """
        
        self.input_file = input_file
        self.output_file = output_file
        self.original_atoms = None
        self.aligned_atoms = None
        print(f"MoleculeProcessor initialized for input '{input_file}', output '{output_file}'.")

    def load_atoms(self):
        """
        Loads the molecular structure from the input file.

        Returns:
            bool: True if the file was successfully loaded, False otherwise.

        Raises:
            FileNotFoundError: If the input file does not exist.
            Exception: For other errors during file reading.
        """
        
        print(f"Loading atoms from '{self.input_file}'...")
        try:
            self.original_atoms = read(self.input_file)
            print(f"Successfully loaded {len(self.original_atoms)} atoms.")
            return True
        except FileNotFoundError:
            print(f"Error: File '{self.input_file}' not found.")
            self.original_atoms = None
            return False
        except Exception as e:
            print(f"An error occurred while reading the file: {e}")
            self.original_atoms = None
            return False

    def align_and_translate(self, atoms):
        """
        Aligns the molecule's principal axes with Cartesian axes and translates the center of mass (COM) to the origin.

        Args:
            atoms (Atoms): The molecular structure to align and translate.

        Returns:
            Atoms: A new `Atoms` object with the aligned and translated structure.
            None: If the alignment or translation fails.

        Notes:
            - For single-atom structures, only the COM translation is performed.
            - If the alignment fails, the method returns None.
        """

        if not isinstance(atoms, Atoms) or len(atoms) == 0:
             print("Warning: Cannot align empty atoms object.")
             return None # Return None to indicate failure
        print("Aligning molecule and translating COM to origin...")
        if len(atoms) == 1:
            print("Warning: Cannot align principal axes for a single atom. Only translating COM.")
            new_atoms = atoms.copy()
            new_atoms.positions -= new_atoms.get_center_of_mass() # Use positions directly
            print("Single atom translated to origin.")
            return new_atoms

        new_atoms = atoms.copy()
        # 1. Translate COM to origin
        try:
             com = new_atoms.get_center_of_mass()
             new_atoms.positions -= com # Use positions directly
        except Exception as e:
            print(f"Warning: Could not calculate COM or translate: {e}")
            # Decide if you want to proceed without centering or return None
            return None # Indicate failure

        # 2. Align principal axes
        try:
            inertia = new_atoms.get_inertia_tensor()
            eigvals, eigvecs = np.linalg.eigh(inertia) # eigvecs columns are eigenvectors
            # Ensure a right-handed coordinate system for the eigenvectors
            if np.linalg.det(eigvecs) < 0.0:
                eigvecs[:, -1] *= -1.0 # Flip the last eigenvector

            # ASE's rotate method needs the vectors defining the rotation.
            # We want to rotate FROM eigvecs TO the identity matrix (Cartesian axes).
            # The 'v' argument in rotate should be the vectors *in their current orientation*
            # The 'a' argument specifies *what they should become*.
            # We rotate the eigenvectors (v=eigvecs) to align with Cartesian axes (a=[(1,0,0),(0,1,0),(0,0,1)])
            new_atoms.rotate(v=eigvecs, a=[(1.0, 0.0, 0.0), (0.0, 1.0, 0.0), (0.0, 0.0, 1.0)])
            print("Molecule aligned to principal axes and centered.")
        except Exception as e:
             print(f"Warning: Could not compute or apply principal axes alignment: {e}")
             # Return centered atoms even if alignment fails
             print("Returning centered (but possibly not aligned) atoms.")
             # Depending on requirements, you might return None here if alignment is critical

        return new_atoms


    def rotate_and_translate_copies(self, yaw, pitch, roll, translation_vector):
        """
        Rotates and translates multiple copies of the aligned molecule.

        Args:
            yaw (list[float]): List of yaw angles (rotation around Z-axis) in radians.
            pitch (list[float]): List of pitch angles (rotation around Y-axis) in radians.
            roll (list[float]): List of roll angles (rotation around X-axis) in radians.
            translation_vector (list[list[float]]): List of translation vectors for each copy.

        Returns:
            Atoms: A combined `Atoms` object containing all rotated and translated copies.
            None: If the operation fails.

        Notes:
            - The lengths of `yaw`, `pitch`, `roll`, and `translation_vector` must match.
            - If alignment was not performed, this method will fail.
        """

        if self.aligned_atoms is None:
            print("Error: Cannot rotate/translate copies because molecule alignment failed or wasn't performed.")
            return None
        if not (len(yaw) == len(pitch) == len(roll) == len(translation_vector)):
             print("Error: Lengths of yaw, pitch, roll, and translation_vector arrays must match.")
             return None

        num_copies = len(yaw)
        print(f"Creating {num_copies} rotated and translated copies...")
        all_atoms_list = [] # Build a list of Atoms objects

        for i in range(num_copies):
            print(f"Processing copy {i+1}/{num_copies}...")
            new_atoms_copy = self.aligned_atoms.copy()

            # Create rotation object from Euler angles (ZYX convention)
            # Ensure angles are in degrees if SciPy expects degrees, or radians if it expects radians (check docs - R.from_euler expects radians by default)
            # Assuming yaw, pitch, roll provided to process_molecules are in RADIANS
            try:
                rotation = R.from_euler('zyx', [yaw[i], pitch[i], roll[i]])
                # Note: new_atoms_copy is already centered from align_and_translate
                # Apply rotation around the origin (current COM)
                new_atoms_copy.positions = rotation.apply(new_atoms_copy.positions)
                # Apply translation
                new_atoms_copy.translate(translation_vector[i])
                all_atoms_list.append(new_atoms_copy)
            except Exception as e:
                 print(f"Error during rotation/translation for copy {i+1}: {e}")
                 # Decide whether to skip this copy or stop entirely
                 return None # Stop if any copy fails

        # Combine all individual Atoms objects into one large Atoms object
        if not all_atoms_list:
            print("Warning: No copies were successfully generated.")
            return Atoms() # Return empty Atoms object

        final_assembly = Atoms()
        for atoms_obj in all_atoms_list:
            final_assembly.extend(atoms_obj)

        print(f"Generated final assembly with {len(final_assembly)} atoms.")
        return final_assembly

    def process_molecules(self,
                          yaw_rad,
                          pitch_rad,
                          roll_rad,
                          translation_vector,
                          relax_molecule: bool = False, # Control if relaxation is done
                          restart_relax: bool = False, # Control if relaxation should *try* to restart
                          tblite_params=None,
                          mh_params=None,
                          totalsteps=20,
                          output_filename_prefix="last"):
        
        """
        Processes the molecules by loading, aligning, rotating, translating, and optionally relaxing them.

        Args:
            yaw_rad (list[float]): List of yaw angles (rotation around Z-axis) in radians.
            pitch_rad (list[float]): List of pitch angles (rotation around Y-axis) in radians.
            roll_rad (list[float]): List of roll angles (rotation around X-axis) in radians.
            translation_vector (list[list[float]]): List of translation vectors for each copy.
            relax_molecule (bool): Whether to perform relaxation on the final structure.
            restart_relax (bool): Whether to restart relaxation from previous calculations.
            tblite_params (dict): Parameters for the TBLite calculator.
            mh_params (dict): Parameters for the Monte Carlo Metropolis-Hastings algorithm.
            totalsteps (int): Maximum number of steps for relaxation.
            output_filename_prefix (str): Prefix for output files generated during relaxation.

        Returns:
            Atoms: The final processed molecular structure (relaxed if requested).
            None: If any step in the process fails.

        Notes:
            - The method writes the initial unrelaxed structure to the output file.
            - If relaxation is requested, the relaxed structure is returned.
        """

        print("\n=== Starting Molecule Processing ===")
        # 1. Load
        if not self.load_atoms():
            print("Processing failed: Could not load input file.")
            return None # Indicate failure

        # 2. Align original molecule
        self.aligned_atoms = self.align_and_translate(self.original_atoms)
        if self.aligned_atoms is None:
             print("Processing failed: Alignment step failed.")
             return None

        # 3. Create rotated/translated assembly
        rotated_atoms = self.rotate_and_translate_copies(yaw_rad, pitch_rad, roll_rad, translation_vector)

        if not isinstance(rotated_atoms, Atoms) or len(rotated_atoms) == 0:
            print("Processing failed: No atoms generated after rotation/translation.")
            return None

        # 4. Write the initial (unrelaxed) assembly
        print(f"Writing initial generated structure to '{self.output_file}'...")
        try:
            write(self.output_file, rotated_atoms)
            print(f"Successfully wrote initial structure with {len(rotated_atoms)} atoms.")
        except Exception as e:
            print(f"Error writing initial structure to '{self.output_file}': {e}")
            # Decide if this error is critical
            return None # Stop if writing initial file fails

        # 5. Optional Relaxation
        final_atoms = rotated_atoms # Default to the unrelaxed structure
        if relax_molecule:
            print(f"\n--- Relaxation Requested (restart={restart_relax}) ---")
            try:
                # Call the relax function from utils.py
                # Pass the 'restart_relax' flag from this method to 'restart_calc' in relax()
                final_atoms = relax(
                    atoms=rotated_atoms.copy(), # Pass a copy to avoid modifying rotated_atoms if relax fails midway
                    output_filename_prefix=output_filename_prefix,
                    tblite_params=tblite_params,
                    mh_params=mh_params,
                    totalsteps=totalsteps,
                    restart_calc=restart_relax # Pass the flag here
                )
                print(f"\n--- Relaxation process finished successfully ---")
                # Note: relax() function handles writing its own specific outputs (_minima.xyz, _summary.png etc)
            except FileNotFoundError as e:
                 # This error is raised by relax() if restart_calc=True but no files found
                 print(f"\n--- Relaxation Error ---")
                 print(f"Relaxation failed: {e}")
                 print("Cannot proceed with processing.")
                 return None # Stop processing if restart failed when requested
            except Exception as e:
                 # Catch other potential errors during relaxation (calculator errors, MH errors)
                 print(f"\n--- Relaxation Error ---")
                 print(f"An unexpected error occurred during the relaxation process: {e}")
                 print("Cannot proceed with processing.")
                 return None # Stop processing if relaxation had a critical error

        print("\n=== Molecule Processing Finished ===")

        # Return the final state (relaxed if relaxation was done, otherwise the initial rotated structure)
        return final_atoms


