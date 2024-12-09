#!/bin/bash
# Molecular Dynamics Simulation Workflow Using GROMACS

# Step 1: Remove water molecules (HOH) from the PDB file to clean it.
grep -v HOH 1L2Y.pdb > 1L2Y_clean.pdb
# This removes all water molecules (HOH) from the input PDB file.

# Step 2: Convert the cleaned PDB file to GROMACS format and generate topology.
gmx pdb2gmx -f 1L2Y_clean.pdb -o 1L2Y_processed.gro -water spce -ignh
# -water spce: Use SPC/E water model for solvation.
# -ignh: Ignore hydrogens in the input PDB file and regenerate them.

# Step 3: Define a cubic simulation box with a 0.8 nm distance to the molecule.
gmx editconf -f 1L2Y_processed.gro -o 1L2Y_newbox.gro -c -d 0.8 -bt cubic
# -c: Center the molecule in the box.
# -d 0.8: Minimum distance of 0.8 nm between the molecule and the box edge.
# -bt cubic: Use a cubic simulation box.

# Step 4: Solvate the box with SPC216 water molecules.
gmx solvate -cp 1L2Y_newbox.gro -cs spc216.gro -o 1L2Y_solvated.gro -p topol.top
# -cp: Input the box containing the molecule.
# -cs: Use the SPC216 water configuration file.
# -p: Update the topology file to include water molecules.

# Step 5: Prepare a tpr file for ion addition to neutralize the system.
gmx grompp -f ions.mdp -c 1L2Y_solvated.gro -p topol.top -o ions.tpr
# -f ions.mdp: Parameter file for ion addition.

# Step 6: Add ions (Na+ and Cl-) to neutralize the system.
echo "SOL" | gmx genion -s ions.tpr -o 1L2Y_ions.gro -p topol.top -pname NA -nname CL -neutral
# Replace solvent molecules (SOL) with ions.
# -pname NA: Specify Na+ as the positive ion.
# -nname CL: Specify Cl- as the negative ion.
# -neutral: Ensure the system is electrically neutral.

# Step 7: Perform energy minimization to remove bad contacts in the system.
gmx grompp -f minim.mdp -c 1L2Y_ions.gro -p topol.top -o em.tpr
gmx mdrun -v -deffnm em
# -f minim.mdp: Parameter file for energy minimization.
# -deffnm em: Set the output file prefix to "em."

# Step 8: Equilibrate the system in the NVT ensemble (constant Number, Volume, Temperature).
gmx grompp -f nvt.mdp -c em.gro -r em.gro -p topol.top -o nvt.tpr
gmx mdrun -v -deffnm nvt -nb gpu -pme gpu -bonded gpu -update gpu -tunepme yes
# -f nvt.mdp: Parameter file for NVT equilibration.
# Use GPU acceleration for faster computations.

# Step 9: Equilibrate the system in the NPT ensemble (constant Number, Pressure, Temperature).
gmx grompp -f npt.mdp -c nvt.gro -r nvt.gro -t nvt.cpt -p topol.top -o npt.tpr
gmx mdrun -v -deffnm npt -nb gpu -pme gpu -bonded gpu -update gpu -tunepme yes
# -t nvt.cpt: Use the checkpoint file from the NVT run.

# Step 10: Perform the production molecular dynamics simulation.
gmx grompp -f md.mdp -c npt.gro -t npt.cpt -p topol.top -o md_0_1.tpr
gmx mdrun -v -deffnm md_0_1 -nb gpu -pme gpu -bonded gpu -update gpu -tunepme yes

# Post-Processing Steps

# Step 11: Center the trajectory and remove periodic boundary effects.
gmx trjconv -s md_0_1.tpr -f md_0_1.xtc -o md_0_1_center.xtc -pbc mol -center
# -pbc mol: Remove periodic boundary conditions for molecules.
# -center: Center the molecule in the simulation box.

# Step 12: Calculate the RMSD (Root Mean Square Deviation) over time.
gmx rms -s md_0_1.tpr -f md_0_1_center.xtc -o rmsd.xvg -tu ns
# -tu ns: Set time units to nanoseconds.

# Step 13: Calculate RMSF (Root Mean Square Fluctuation) for each residue.
gmx rmsf -s md_0_1.tpr -f md_0_1_center.xtc -o rmsf.xvg -res
# -res: Calculate per-residue fluctuations.

# Step 14: Calculate the Radius of Gyration over time.
gmx gyrate -s md_0_1.tpr -f md_0_1_center.xtc -o gyrate.xvg
# Measure the compactness of the molecule.

# Step 15: Calculate the number of hydrogen bonds over time.
gmx hbond -s md_0_1.tpr -f md_0_1_center.xtc -num hydrogen_bonds.xvg -tu ns
# -num: Output the number of hydrogen bonds.

# Step 16: Calculate the Solvent Accessible Surface Area (SASA) over time.
gmx sasa -s md_0_1.tpr -f md_0_1_center.xtc -o area.xvg -tu ns
# -tu ns: Output time units in nanoseconds.
# The the md_simulation file
chmod +x md_simulation.sh
# Run the simulation workflow script.
./md_simulation.sh