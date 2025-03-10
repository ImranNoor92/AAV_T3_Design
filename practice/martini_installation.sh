#!/bin/bash
echo "$RANDOM"  # Supported in bash. No warnings.
######################
#Install martini 3.0.0 force field and vermouth package for parameterization.
#Step 1: Clone the MARTINI Repository
git clone https://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@github.com/ImranNoor92/martini_ff.git
#Purpose: Clones the MARTINI repository from GitHub using your personal access token (PAT) for authentication.
#Result: A new folder named `martini_ff` is created, containing the MARTINI files.

#Step 2: Create and Navigate to the Working Directory
mkdir martini_tutorial
cd martini_tutorial/
#Purpose: Creates a new directory named `martini_tutorial` and navigates into it.
#Purpose: Downloads the protein structure file `181L.pdb` from the Protein Data Bank (PDB).

#Step 3: Clean the Protein Structure File
grep "^ATOM" 181L.pdb > 181L_clean.pdb
#Purpose: Extracts only the atomic coordinates (`ATOM` lines) from the `181L.pdb` file, removing any irrelevant information (e.g., HETATM, CONECT lines).

#Step 4: 
sudo apt update
sudo apt install python3-pip
pip3 --version
pip install vermouth
#Purpose: Installs the Python package manager `pip` and the `vermouth` package for MARTINI force field parameterization.

#Step 5: Install Martinize2 Using Vermouth.
pip install git+https://github.com/marrink-lab/vermouth-martinize.git#vermouth
#Purpose: Installs the MARTINI2 tool `vermouth-martinize` using the `vermouth` package.

#Step 6: Add `martinize2` to PATH
nano ~/.bashrc
export PATH=$PATH:~/.local/bin
source ~/.bashrc
#Purpose: Adds the `martinize2` executable to the system PATH for easy access. 
#Step 7: Test `martinize2`
martinize2 -h
#Purpose: Checks if the `martinize2` tool is successfully installed and accessible.

#Step 8: DSSP Installation
git clone https://github.com/cmbi/dssp.git
cd dssp
./autogen.sh
./configure
make
sudo apt install dssp
dssp --version
# Purpose: Compiles DSSP from source using the Autotools build system.
#   ./autogen.sh: Prepares the build environment.
#   ./configure: Configures the build system.
#   make: Compiles the DSSP binary.
# Result: Builds the dssp executable.
#Verifies the DSSP installation and displays the version.

