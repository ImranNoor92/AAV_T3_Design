#!/bin/bash
echo "$RANDOM"  # Supported in bash. No warnings.
#This file will be used to write and store codes for removing periodic boundary effects in the molecular dynamics simulation trajectory.
#PBC: defines the boundary conditions of the simulation box, which can introduce artifacts in the trajectory due to molecules crossing box boundaries.
#The gmx trjconv command in GROMACS can be used to remove these artifacts and center the trajectory on the molecule of interest.
#The -pbc mol flag specifies that the periodic boundary conditions should be removed for molecules, and the -center flag centers the molecule in the simulation box.

# === Step 20 : Center the Protein and Remove PBC Artifacts ===
echo "Centering the protein and removing PBC artifacts..."
gmx trjconv -f dynamic.xtc -s dynamic.tpr -o centered.xtc -center -pbc mol
# Explanation:
# - `gmx trjconv`: Processes trajectory files to remove artifacts or reformat data.
# - `-f dynamic.xtc`: Input trajectory file.
# - `-s dynamic.tpr`: Reference structure and topology file.
# - `-o centered.xtc`: Output trajectory file with the protein centered and PBC artifacts removed.
# - `-center`: Centers the protein in the simulation box.
# - `-pbc mol`: Ensures each molecule (e.g., the protein) is kept whole across PBC boundaries.
# Group selection:
# - Select "Protein" (1) to center the protein.
# - Select "System" (0) for output.

# Interactive Selection:
# 1. Select the group for centering:
#    - Choose `Protein` (or your custom group, e.g., `Protein_Backbone_Beads`) when prompted.
# 2. Select the group to output:
#    - Choose `System` to include the protein and solvent in the output.

# === Step 21 : Process for Visualization ===
echo "Processing trajectory for visualization..."
gmx trjconv -f centered.xtc -s dynamic.tpr -o viz.xtc -fit rot+trans
# Explanation:
# - `-fit rot+trans`: Aligns the trajectory by removing overall rotation and translation of the system.
# - `-o viz.xtc`: Output trajectory file optimized for visualization.
# - Select "Protein" (1) to center the protein.
# - Select "System" (0) for output.

# === Step 22 : Visualize the Corrected Trajectory ===
echo "Visualizing centered trajectory in VMD..."
vmd equilibration.gro viz.xtc
# Explanation:
# - Opens the trajectory (`viz.xtc`) in VMD with the initial structure (`equilibration.gro`) as reference.


# === Step 23: Use PyMol ===
echo "Converting trajectory to PDB format for PyMOL..."
gmx trjconv -f centered.xtc -s dynamic.tpr -o viz.pdb -conect
pymol viz.pdb
# Explanation:
# - Converts the trajectory to PDB format for compatibility with PyMOL.
# - Adds bond information using the `-conect` flag.#This step is important for analyzing and visualizing the trajectory data without artifacts from periodic boundary effectsecho "Creating index file for RMSD of Protein and Backbone Beads..."
gmx make_ndx -f equilibration.gro -o protein_bb.ndx
#alternatively
gmx rms -s dynamic.tpr -f dynamic.xtc -n protein_bb.ndx -o rmsd_protein_bb.xvg -b 20000 -e 60000

# Explanation:
# - `gmx make_ndx`: Generates or modifies index files.
# - `-f equilibration.gro`: Input structure file from the equilibration step.
# - `-o protein_bb.ndx`: Outputs an index file for protein and backbone beads.

# Interactive Selections:
# 1. Select the `Protein` group:
#    > 1 (Protein)
# 2. Add Backbone Beads (Martini backbone typically uses BB bead type):
#    > a BB
# 3. Combine Protein and Backbone Beads groups:
#    > 1 | a BB
# 4. Name the combined group:
#    > name 15 Protein_Backbone_Beads
# 5. Save and quit:
#    > qecho "Creating index file for RMSF by residues..."
gmx make_ndx -f equilibration.gro -o rmsf_residues.ndx
# Explanation:
# - `gmx make_ndx`: Creates or modifies index files.
# - `-f equilibration.gro`: Input structure file from the equilibration step.
# - `-o rmsf_residues.ndx`: Outputs an index file for residue-specific RMSF analysis.

# Interactive Selections:
# 1. Select all residues of interest (e.g., residues 1–100):
#    > r 1-100
# 2. Name the group:
#    > name 16 Residues_1_100
# 3. Save and quit:
#    > qecho "Calculating RMSD for Protein and Backbone Beads..."
gmx rms -s dynamic.tpr -f centered.xtc -n protein_bb.ndx -o rmsd_protein_bb.xvg -b 20000 -e 60000
# Explanation:
# - `gmx rms`: Computes Root Mean Square Deviation (RMSD) for the specified group.
# - `-s dynamic.tpr`: Reference structure file for alignment.
# - `-f centered.xtc`: PBC-corrected and centered trajectory file.
# - `-n protein_bb.ndx`: Index file with the combined Protein and Backbone Beads group.
# - `-o rmsd_protein_bb.xvg`: Output RMSD values file.
# - `-b 20000 -e 60000`: Restricts RMSD calculation to 20–60 ns.echo "Calculating RMSF for residues..."
gmx rmsf -s dynamic.tpr -f centered.xtc -n rmsf_residues.ndx -o rmsf_residues.xvg -b 20000 -e 60000
# Explanation:
# - `gmx rmsf`: Calculates Root Mean Square Fluctuation (RMSF) per residue.
# - `-s dynamic.tpr`: Reference structure.
# - `-f centered.xtc`: PBC-corrected trajectory file.
# - `-n rmsf_residues.ndx`: Index file for the residue group.
# - `-o rmsf_residues.xvg`: Output file for residue-specific RMSF values.
# - `-b 20000 -e 60000`: Restricts RMSF calculation to 20–60 ns.echo "Analyzing hydrogen bonds between selected groups..."
gmx hbond -s dynamic.tpr -f centered.xtc -n protein_bb.ndx -num hbonds.xvg -b 20000 -e 60000
# Explanation:
# - `gmx hbond`: Identifies hydrogen bonds between specified groups.
# - `-s dynamic.tpr`: Reference topology and structure file.
# - `-f centered.xtc`: PBC-corrected trajectory file.
# - `-n protein_bb.ndx`: Index file defining the groups of interest (e.g., Protein and Backbone Beads).
# - `-num hbonds.xvg`: Outputs the number of hydrogen bonds over time.
# - `-b 20000 -e 60000`: Restricts H-bond analysis to the 20–60 ns time frame.echo "Processing trajectory for visualization..."
gmx trjconv -f centered.xtc -s dynamic.tpr -o viz.xtc -fit rot+trans
# Explanation:
# - `gmx trjconv`: Aligns and processes trajectory for visualization.
# - `-fit rot+trans`: Removes global rotation and translation.
# - `-o viz.xtc`: Outputs the aligned trajectory file for visualization.echo "Visualizing trajectory with VMD..."
vmd equilibration.gro viz.xtc
# Explanation:
# - Loads the aligned trajectory (`viz.xtc`) into VMD using the initial structure file (`equilibration.gro`).echo "Visualizing trajectory with VMD..."
vmd equilibration.gro viz.xtc
# Explanation:
# - Loads the aligned trajectory (`viz.xtc`) into VMD using the initial structure file (`equilibration.gro`).echo "Converting trajectory to PDB for PyMOL..."
gmx trjconv -f centered.xtc -s dynamic.tpr -o viz.pdb -conect
pymol viz.pdb
# Explanation:
# - Converts the trajectory to PDB format for compatibility with PyMOL.
# - `-conect`: Adds bond information to the PDB file..