#!/bin/bash
# This script sets up and runs a Martini coarse-grained molecular dynamics simulation using GROMACS.

# === Step 1: Download the Protein Structure File ===
echo "Downloading the protein structure file..."
wget http://www.rcsb.org/pdb/files/181L.pdb
# Purpose:
# Downloads the protein structure file (PDB format) from the Protein Data Bank (PDB).

# === Step 2: Clean the Protein Structure File ===
echo "Cleaning the protein structure file..."
grep "^ATOM" 181L.pdb > 181L_clean.pdb
# Explanation:
# - Removes non-atom entries (e.g., HETATM or REMARK) to retain only atomic coordinates for the protein.
# - Outputs a cleaned PDB file named `181L_clean.pdb`.

# === Step 3: Coarse-Grain the Protein Structure ===
echo "Coarse-graining the protein structure..."
martinize2 -f 181L_clean.pdb -o t4l_only.top -x t4l_cg.pdb -p backbone -ff martini3001 \
-ss CCHHHHHHHHHCCEEEEEECTTSCEEEETTEEEESSSCHHHHHHHHHHHHTSCCTTBCCHHHHHHHHHHHHHHHHHHHHHCTTTHHHHHHSCHHHHHHHHHHHHHHHHHHHHTCHHHHHHHHTTCHHHHHHHHHSSHHHHHSHHHHHHHHHHHHHSSSGGGC
# Explanation:
# - `-f 181L_clean.pdb`: Input cleaned PDB file.
# - `-o t4l_only.top`: Output topology file for the coarse-grained protein.
# - `-x t4l_cg.pdb`: Output coarse-grained structure file in PDB format.
# - `-p backbone`: Specifies the treatment of the protein backbone with Martini parameters.
# - `-ff martini3001`: Uses Martini 3.0.0 force field.
# - `-ss`: Specifies secondary structure using DSSP format.

# === Step 4: Generate the Simulation Box ===
echo "Generating the simulation box..."
gmx editconf -f t4l_cg.pdb -d 2.0 -bt cubic -o t4l_cg.gro
# Explanation:
# - `gmx editconf`: Creates or modifies a simulation box.
# - `-f t4l_cg.pdb`: Input coarse-grained structure file.
# - `-d 2.0`: Adds 2.0 nm padding between the protein and box edges.
# - `-bt cubic`: Generates a cubic simulation box.
# - `-o t4l_cg.gro`: Outputs the structure with the defined box in `.gro` format.

# === Step 5: Preprocess for Energy Minimization ===
echo "Preprocessing for energy minimization..."
gmx grompp -p t4l_only.top -f minimization.mdp -c t4l_cg.gro -o minimization-vac.tpr -r t4l_cg.gro
# Explanation:
# - `gmx grompp`: Prepares input for molecular dynamics.
# - `-p t4l_only.top`: Topology file defining the system components.
# - `-f minimization.mdp`: Parameter file for minimization.
# - `-c t4l_cg.gro`: Input structure file from the previous step.
# - `-o minimization-vac.tpr`: Output file for minimization.
# - `-r t4l_cg.gro`: Reference structure for position restraints.

# === Step 6: Run Energy Minimization ===
echo "Running energy minimization..."
gmx mdrun -deffnm minimization-vac -v
# Explanation:
# - `gmx mdrun`: Runs the molecular dynamics simulation.
# - `-deffnm minimization-vac`: Sets output file prefix.
# - `-v`: Enables verbose output.

# === Step 7: Solvate the System ===
echo "Adding solvent to the system..."
gmx solvate -cp minimization-vac.gro -cs water.gro -radius 0.21 -o solvated.gro
# Explanation:
# - `gmx solvate`: Adds water molecules to the simulation box.
# - `-cp minimization-vac.gro`: Input minimized structure.
# - `-cs water.gro`: Water box compatible with Martini force field.
# - `-radius 0.21`: Minimum distance between solute and solvent beads.
# - `-o solvated.gro`: Output solvated structure.

# === Step 8: Preprocess for Equilibration ===
echo "Preparing equilibration input..."
gmx grompp -p system.top -c solvated.gro -f equilibration.mdp -o equilibration.tpr -r solvated.gro
# Explanation:
# - `-p system.top`: Topology file updated with solvent.
# - `-c solvated.gro`: Input solvated structure file.
# - `-f equilibration.mdp`: Parameter file for equilibration.
# - `-o equilibration.tpr`: Output file for equilibration.

# === Step 9: Run Equilibration ===
echo "Running equilibration..."
gmx mdrun -deffnm equilibration -v
# Explanation:
# - Runs the molecular dynamics simulation for equilibration.

# === Step 10: Production MD Simulation ===
echo "Running production MD simulation..."
gmx grompp -p system.top -c equilibration.gro -f dynamic.mdp -o dynamic.tpr -r equilibration.gro
gmx mdrun -deffnm dynamic -v
# Explanation:
# - `gmx grompp`: Prepares input for production MD.
# - `gmx mdrun`: Executes the production simulation.

# === Step 11: Post-Processing and Visualization ===
echo "Processing trajectories for visualization..."
gmx trjconv -f dynamic.xtc -s dynamic.tpr -fit rot+trans -o viz.xtc
vmd equilibration.gro viz.xtc
# Explanation:
# - `gmx trjconv`: Processes trajectory files (e.g., removes periodic boundary artifacts).
# - `-fit rot+trans`: Aligns the system to the reference structure.
# - `vmd`: Visualizes the trajectory using VMD.

# === Step 12: Convert the .gro File to .pdb with Bond Information ===
echo "Converting the equilibration structure to PDB with bond information..."
gmx trjconv -f equilibration.gro -s dynamic.tpr -pbc whole -o equilibration.pdb -conect
# Explanation:
# - `gmx trjconv`: Converts the structure or trajectory file into different formats.
# - `-f equilibration.gro`: Input structure file from equilibration.
# - `-s dynamic.tpr`: TPR file with topology and bonding information.
# - `-pbc whole`: Reconstructs molecules across periodic boundaries.
# - `-o equilibration.pdb`: Output PDB file.
# - `-conect`: Adds bond information to the PDB file.

echo "Removing the ENDMDL line from the PDB file..."
sed "/ENDMDL/d" -i equilibration.pdb
# Explanation:
# - `sed`: A stream editor used to process and transform text.
# - `"/ENDMDL/d"`: Deletes lines containing "ENDMDL" in the file.
# - `-i`: Edits the file in place (no need to create a temporary file).

# === Step 13: Visualize the Trajectory Using VMD ===
echo "Visualizing the trajectory in VMD..."
vmd equilibration.pdb viz.xtc
# Explanation:
# - `vmd`: Launches the Visual Molecular Dynamics (VMD) software for visualization.
# - `equilibration.pdb`: Structure file with bonds included.
# - `viz.xtc`: Processed trajectory file.

# === Step 14: CConvert the Trajectory to PDB for PyMOL Visualization ===

echo "Converting trajectory to PDB format for PyMOL..."
gmx trjconv -f dynamic.xtc -s dynamic.tpr -fit rot+trans -o viz.pdb -conect
# Explanation:
# - `-f dynamic.xtc`: Input trajectory file.
# - `-s dynamic.tpr`: Input TPR file with topology and reference structure.
# - `-fit rot+trans`: Aligns the trajectory by removing rotation and translation.
# - `-o viz.pdb`: Output trajectory in PDB format.
# - `-conect`: Adds bond information for better visualization.

echo "Opening the PDB trajectory in PyMOL..."
pymol viz.pdb
# Explanation:
# - `pymol`: Launches the PyMOL visualization tool with the PDB trajectory.

# === Step 15 : Index File 1: For RMSD of Protein and Backbone Beads ===
echo "Creating index file for RMSD of Protein and Backbone Beads..."
gmx make_ndx -f equilibration.gro -o protein_bb.ndx
# Explanation:
# - `gmx make_ndx`: GROMACS tool for creating or modifying index files.
# - `-f equilibration.gro`: Input structure file from equilibration.
# - `-o protein_bb.ndx`: Output index file for analysis.

# Interactive Selections:
# 1. Select `Protein` group:
#    > 1 (Protein)
# 2. Add Backbone Beads (Martini backbone typically uses BB bead type):
#    > a BB
# 3. Combine the Protein and Backbone Beads groups:
#    > 1 | a BB
# 4. Name the combined group:
#    > name 15 Protein_Backbone_Beads
# 5. Save and quit:
#    > q

# === Step 16 : Index File 2: For RMSF by Residues ===
echo "Creating index file for RMSF by residues..."
gmx make_ndx -f equilibration.gro -o rmsf_residues.ndx
# Explanation:
# - `-f equilibration.gro`: Input structure file.
# - `-o rmsf_residues.ndx`: Index file for residue-based RMSF analysis.

# Interactive Selections:
# 1. Select all residues of interest (e.g., residues 1–100):
#    > r 1-100
# 2. Name the group:
#    > name 16 Residues_1_100
# 3. Save and quit:
#    > q

# === Step 17 : RMSD for Protein and Backbone Beads ===
echo "Calculating RMSD for Protein and Backbone Beads..."
gmx rms -s dynamic.tpr -f dynamic.xtc -n protein_bb.ndx -o rmsd_protein_bb.xvg -b 20000 -e 60000
# Explanation:
# - `gmx rms`: Calculates RMSD over the trajectory.
# - `-s dynamic.tpr`: Reference structure for alignment.
# - `-f dynamic.xtc`: Trajectory file.
# - `-n protein_bb.ndx`: Index file with the combined group for protein and backbone beads.
# - `-o rmsd_protein_bb.xvg`: Output file for RMSD values.
# - `-b 20000 -e 60000`: Restricts analysis to 20–60 ns.

# === Step 18 : RMSF for Residues ===
echo "Calculating RMSF for residues..."
gmx rmsf -s dynamic.tpr -f dynamic.xtc -n rmsf_residues.ndx -o rmsf_residues.xvg -b 20000 -e 60000
# Explanation:
# - `gmx rmsf`: Calculates RMSF per residue.
# - `-s dynamic.tpr`: Reference structure for alignment.
# - `-f dynamic.xtc`: Trajectory file.
# - `-n rmsf_residues.ndx`: Index file with the residue group.
# - `-o rmsf_residues.xvg`: Output file for RMSF values.
# - `-b 20000 -e 60000`: Restricts analysis to 20–60 ns.

# === Step 19 : Identify H-bonds Between Groups ===
echo "Analyzing hydrogen bonds between selected groups..."
gmx hbond -s dynamic.tpr -f dynamic.xtc -n protein_bb.ndx -num hbonds.xvg -b 20000 -e 60000
# Explanation:
# - `gmx hbond`: Identifies hydrogen bonds between groups.
# - `-s dynamic.tpr`: Reference structure.
# - `-f dynamic.xtc`: Trajectory file.
# - `-n protein_bb.ndx`: Index file containing the groups of interest (e.g., Protein and Backbone Beads).
# - `-num hbonds.xvg`: Output file with the number of H-bonds over time.
# - `-b 20000 -e 60000`: Restricts analysis to 20–60 ns.

# === Step 20 : ===


# === Step 21 : ===


# === Step 22 : ===