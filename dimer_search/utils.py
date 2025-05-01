import os
import glob
import ase
from ase.io import read, write
from ase.io.trajectory import Trajectory
from ase.optimize.minimahopping import MinimaHopping
from ase.atoms import Atoms
import logging
logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)

logger.info("TBLite found and imported.")
logger.warning("Warning: tblite not found.")

# Assuming MHPlot is available
try:
    from ase.optimize.minimahopping import MHPlot
    print("MHPlot found and imported.")
except ImportError:
    print("Warning: MHPlot not found. Plotting will be skipped.")
    MHPlot = None # Placeholder


def relax(
    atoms,
    tblite_params=None,
    mh_params=None,
    totalsteps=20,
    output_filename_prefix="last",
    restart_calc: bool = False # Parameter to control restart
):
    """
    Relaxes the molecule using TBLite and MinimaHopping with options.

    Args:
        atoms (ase.Atoms): Initial atoms object (used if restart_calc is False
                           or if restart fails when requested).
        tblite_params (dict, optional): Parameters for the TBLite calculator.
        mh_params (dict, optional): Parameters for MinimaHopping.
        totalsteps (int): Number of MinimaHopping steps to perform in this run.
        output_filename_prefix (str): Prefix for output files (plot, final XYZ).
        restart_calc (bool): If True, attempt to restart from existing
                             'minima.traj' or 'qn*.traj' files. If False (default),
                             start from the provided 'atoms' object.
    """
    # --- Parameter Handling (Defaults & Merging) ---
    if tblite_params is None: tblite_params = {}
    if mh_params is None: mh_params = {}

    # Default TBLite Parameters
    default_tblite_params = { "method": "GFN2-xTB", "verbosity": 0 }
    tblite_params = {**default_tblite_params, **tblite_params}

    # Default MinimaHopping Parameters
    default_mh_params = {
        "T0": 1000.0, "beta1": 1.1, "beta2": 1.1, "beta3": 1.0 / 1.1, "Ediff0": 0.5,
        "alpha1": 0.98, "alpha2": 1.0 / 0.98, "mdmin": 2, "logfile": "hop.log",
        "minima_threshold": 0.5, "timestep": 1.0, "minima_traj": "minima.traj",
        "fmax": 0.05,
    }
    mh_params = {**default_mh_params, **mh_params}

    # --- Determine Starting Atoms based on restart_calc ---
    atoms_to_run = None
    initial_atoms_source = "Unknown"

    if restart_calc:
        print("\n--- Relaxation Restart requested (restart_calc=True) ---")
        restart_atoms = None
        minima_file_to_check = mh_params.get('minima_traj', 'minima.traj')
        qn_pattern = 'qn*.traj' # Standard MH quench file pattern

        # 1. Check target minima.traj
        print(f"Checking for non-empty '{minima_file_to_check}'...")
        if os.path.exists(minima_file_to_check) and os.path.getsize(minima_file_to_check) > 0:
            try:
                print(f"Attempting restart from last frame of '{minima_file_to_check}'.")
                restart_atoms = read(minima_file_to_check, index=-1)
                initial_atoms_source = f"'{minima_file_to_check}' (last frame)"
                print(f"Successfully loaded restart configuration for relaxation.")
            except Exception as e:
                print(f"Warning: Could not read from '{minima_file_to_check}' for restart. Error: {e}")
                restart_atoms = None

        # 2. If no restart from minima.traj, check qn*.traj files
        if restart_atoms is None:
            print(f"\n'{minima_file_to_check}' not used for restart. Checking '{qn_pattern}'...")
            qn_files = glob.glob(qn_pattern)
            if qn_files:
                qn_files.sort()
                print(f"Found potential quench files: {qn_files}")
                for qn_file in reversed(qn_files):
                    print(f"Checking '{qn_file}'...")
                    if os.path.exists(qn_file) and os.path.getsize(qn_file) > 0:
                        try:
                            print(f"Attempting restart from last frame of '{qn_file}'.")
                            restart_atoms = read(qn_file, index=-1)
                            initial_atoms_source = f"'{qn_file}' (last frame)"
                            print(f"Successfully loaded restart configuration for relaxation.")
                            break
                        except Exception as e:
                            print(f"Warning: Could not read from '{qn_file}' for restart. Error: {e}")
                            restart_atoms = None
                    else:
                         print(f"'{qn_file}' is empty or does not exist.")
            else:
                print(f"No files found matching '{qn_pattern}'.")

        # 3. Assign atoms_to_run or Raise Error if restart failed
        if restart_atoms is not None and isinstance(restart_atoms, Atoms):
            print(f"\n--- Using configuration from {initial_atoms_source} for relaxation restart ---")
            atoms_to_run = restart_atoms
        else:
            raise FileNotFoundError(
                f"Restart requested for relaxation (restart_calc=True) but failed. "
                f"Could not find or read a valid configuration from "
                f"'{minima_file_to_check}' or any '{qn_pattern}' files."
            )
    else:
        # restart_calc is False (default)
        print("\n--- Starting relaxation from provided configuration (restart_calc=False) ---")
        atoms_to_run = atoms # Use the input atoms directly
        initial_atoms_source = "Provided to relax function"

    # --- Sanity Check ---
    if atoms_to_run is None or not isinstance(atoms_to_run, Atoms):
         raise ValueError("Critical Error: Could not determine a valid Atoms object for relaxation.")

    # --- Relaxation Logic ---
    print(f"\nInitiating relaxation using atoms from: {initial_atoms_source}")
    # ... (rest of the relaxation logic: setting calculator, running hop, post-processing) ...
    print(f"Using TBLite parameters: {tblite_params}")
    print(f"Using MinimaHopping parameters: {mh_params}")
    print(f"Total steps for this MinimaHopping run: {totalsteps}")

    calculator = TBLite(**tblite_params)
    atoms_to_run.set_calculator(calculator)

    hop = MinimaHopping(atoms=atoms_to_run, **mh_params)
    try:
        hop(totalsteps=totalsteps)
        print("\nMinimaHopping run completed.")
    except Exception as e:
        print(f"\nError during MinimaHopping execution: {e}")
        raise e # Re-raise the error to signal failure

    # --- Post-processing ---
    print("\n--- Relaxation Post-processing ---")
    if MHPlot:
        try:
            mhplot = MHPlot()
            plot_filename = f"{output_filename_prefix}_summary.png"
            mhplot.save_figure(plot_filename)
            print(f"Saved MinimaHopping summary plot to '{plot_filename}'")
        except FileNotFoundError:
             print(f"Warning: Log file '{mh_params.get('logfile', 'hop.log')}' not found. Cannot generate plot.")
        except Exception as e:
            print(f"Warning: Could not generate or save MHPlot. Error: {e}")
    else:
        print("MHPlot not available, skipping plot generation.")

    final_minima_traj_file = mh_params.get('minima_traj', 'minima.traj')
    final_output_xyz = f"{output_filename_prefix}_minima.xyz"
    try:
        if os.path.exists(final_minima_traj_file) and os.path.getsize(final_minima_traj_file) > 0:
            traj = Trajectory(final_minima_traj_file, 'r')
            if len(traj) > 0:
                final_minimum_atoms = traj[-1]
                write(final_output_xyz, final_minimum_atoms)
                print(f"Saved final minimum structure from '{final_minima_traj_file}' (this run) to '{final_output_xyz}'")
                # Return the specific atoms object found as the minimum
                atoms_to_run = final_minimum_atoms
            else:
                 print(f"Warning: Trajectory file '{final_minima_traj_file}' is empty after the run.")
        else:
             print(f"Warning: Final minima trajectory file '{final_minima_traj_file}' not found or empty after the run.")
    except Exception as e:
        print(f"Warning: Could not read final trajectory '{final_minima_traj_file}' or write XYZ file '{final_output_xyz}'. Error: {e}")

    # Return the atoms object corresponding to the final minimum found
    return atoms_to_run
