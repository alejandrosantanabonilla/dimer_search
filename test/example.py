from dimer_search import *  

# Example usage:
processor = DimerProcessor(pdb_filename='dimer.pdb', atom_indices=[0, 1, 2])
joined_molecule = processor.rotate_and_join(z_angle=35)

relaxed_molecule = processor.relax(method="GFN1-xTB", Ediff0=0.5, T0=3000, totalsteps=2)
print("Relaxation complete. Files written.")

relaxed_molecule = processor.relax()  # Using default values
print("Relaxation complete. Files written.")
