#!/bin/bash
# AAV2 (1LP3) Full Capsid Assembly Simulation with MARTINI in GROMACS

# === Step 1: Download the Full AAV2 Capsid PDB ===
echo "Downloading AAV2 capsid PDB..."
wget https://files.rcsb.org/download/1LP3.pdb -O aav2.pdb

# === Step 2: Clean the PDB File ===
echo "Cleaning the PDB file..."
grep "^ATOM" aav2.pdb > aav2_clean.pdb

# === Step 3: Coarse-Grain the Full Capsid ===
echo "Converting AAV2 capsid to coarse-grained MARTINI model..."
martinize2 -f aav2_clean.pdb -o aav2_cg.top -x aav2_cg.pdb -ff martini3001 -p backbone

# === Step 4: Define a Simulation Box ===
echo "Defining a large simulation box for the capsid..."
gmx editconf -f aav2_cg.pdb -d 5.0 -bt dodecahedron -o aav2_capsid_box.gro

# === Step 5: Add MARTINI-Compatible Solvent ===
echo "Adding MARTINI water..."
gmx solvate -cp aav2_capsid_box.gro -cs water.gro -o solvated.gro

# === Step 6: Energy Minimization ===
echo "Preparing for energy minimization..."
gmx grompp -f minimization.mdp -c solvated.gro -p aav2_cg.top -o min.tpr
gmx mdrun -deffnm min -v

# === Step 7: Equilibration ===
echo "Running equilibration..."
gmx grompp -f equilibration.mdp -c min.gro -p aav2_cg.top -o equilibration.tpr
gmx mdrun -deffnm equilibration -v

# === Step 8: Production MD Simulation ===
echo "Running long production MD simulation..."
gmx grompp -f production.mdp -c equilibration.gro -p aav2_cg.top -o production.tpr
gmx mdrun -deffnm production -v

# === Step 9: Center the Protein and Remove PBC Artifacts ===
echo "Centering the protein and removing PBC artifacts..."
gmx trjconv -f production.xtc -s production.tpr -o centered.xtc -center -pbc mol

# === Step 10: Process for Visualization ===
echo "Processing trajectory for visualization..."
gmx trjconv -f centered.xtc -s production.tpr -o viz.xtc -fit rot+trans
vmd equilibration.gro viz.xtc

# === Step 11: Create Index Files for Analysis ===
echo "Creating index file for RMSD of Protein and Backbone Beads..."
gmx make_ndx -f equilibration.gro -o protein_bb.ndx
echo "Creating index file for RMSF by residues..."
gmx make_ndx -f equilibration.gro -o rmsf_residues.ndx

# === Step 12: Compute RMSD for Protein and Backbone Beads ===
echo "Calculating RMSD for Protein and Backbone Beads..."
gmx rms -s production.tpr -f centered.xtc -n protein_bb.ndx -o rmsd_protein_bb.xvg -b 20000 -e 60000

# === Step 13: Compute RMSF for Residues ===
echo "Calculating RMSF for residues..."
gmx rmsf -s production.tpr -f centered.xtc -n rmsf_residues.ndx -o rmsf_residues.xvg -b 20000 -e 60000

# === Step 14: Cluster Analysis for Capsid Assembly ===
echo "Analyzing capsid self-assembly..."
gmx cluster -f production.xtc -s production.tpr -o clusters.xvg

echo "Simulation complete. Check visualization and clustering results."
