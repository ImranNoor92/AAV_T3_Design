# **🚀 Toni Lab Workstation Guide**

## **🔗 Connecting to the UM Workstation**
To access the **Toni Lab workstation**, follow these steps:

1. **Open a terminal** on your local machine.  
2. **Connect to the UM VPN** (required for remote access).  
3. **Login to the workstation** by entering:  
   ```bash
   ssh mxn833@asbio01.as.miami.edu
   ```
4. Enter your **password** when prompted.

### **📂 Navigating to Your Working Directory**
Once logged in, switch to your working directory:
```bash
cd /data/data/noor
```

---

## **📤 Transferring Files Between Local Machine and Workstation**
To **copy files from your laptop to asbio01**, use:
```bash
scp -r aav2-monomer-charm mxn833@asbio01.as.miami.edu:/data/data/noor
```

To **copy files from asbio01 to your local machine**, use:
```bash
scp -r mxn833@asbio01.as.miami.edu:/data/data/noor/aav2-monomer-charm .
```
💡 *Note: Ensure there is a space before the period (.) at the end!*

---

## **🛠 Running GROMACS Simulations**
### **🔧 Example `grompp` Command**
To preprocess your simulation files:
```bash
gmx_mpi grompp
```
⚠️ *If you encounter warnings, you may need to use:*
```bash
gmx_mpi grompp -maxwarn -1
```

### **🏃 Running the `mdrun` Command**
Execute a simulation with:
```bash
nohup gmx_mpi mdrun -rdd 3.0 -deffnm step4.0_minimization &
```
- `nohup` ensures the process continues running even if you disconnect.
- The ampersand (`&`) runs the job in the background.

---

## **📜 Checking Simulation Progress & Issues**
To **check which step your simulation is on**, run:
```bash
cat step4.0_minimization.log
```

To **check for issues and monitor max force** (useful during minimization):
```bash
cat nohup.out
```

---

## **⚙️ Modifying MDP Files**
The **MDP (Molecular Dynamics Parameter) file** controls your simulation settings.  
- Adjust **`nsteps`** (total simulation steps) or **`emtol`** (energy minimization tolerance).
- You can only have **one** active at a time.

---

## **🎬 Production Simulation Guidelines**
### **✅ Key Parameters**
- **Timestep**: Keep at `0.02 ps` *(20 fs)*.
- **Total Simulation Time**: Run the **full capsid system** for **2 microseconds**.
- **Steps (`nsteps`)**: Adjust accordingly to reach 2 μs.
- **Frame Output Frequency**:  
  - Ensure **2000 frames over 2 μs**.
  - Modify `nstxout-compressed` to control output frames.
  - Keep other `nst...` parameters the same as `nstxout-compressed`.
  - Maintain `compressed-x-precision` as is.

🚨 **IMPORTANT:** When the **production run** starts, inform the supervisor for verification.

---

## **📊 Post-Simulation Analysis**
Once the simulation **completes**, calculate:

1. **RMSD** (*Root Mean Square Deviation*) of **capsid backbone beads** over **2 μs**.  
2. **Radius of Gyration** for **all beads** over **2 μs**.

---

## **📌 Next Steps (Optional, Not Urgent)**
After production, for additional analysis:

✔ **Pick two monomers** and create sections in the index file.  
✔ Compute **average distance** between their **center of mass (COM)** using **backbone beads** over the **last 1 μs** of simulation.

---

## **📌 Quick Reference Commands**
| Task | Command |
|------|---------|
| Login to workstation | `ssh mxn833@asbio01.as.miami.edu` |
| Copy files (local → workstation) | `scp -r <folder> mxn833@asbio01.as.miami.edu:/data/data/noor` |
| Copy files (workstation → local) | `scp -r mxn833@asbio01.as.miami.edu:/data/data/noor/<folder> .` |
| Run GROMACS pre-processing (`grompp`) | `gmx_mpi grompp` |
| Run GROMACS with warnings allowed | `gmx_mpi grompp -maxwarn -1` |
| Run simulation (`mdrun`) | `nohup gmx_mpi mdrun -rdd 3.0 -deffnm step4.0_minimization &` |
| Check simulation step | `cat step4.0_minimization.log` |
| Check issues and max force | `cat nohup.out` |

---

## **🎯 Final Notes**
- Keep **simulation parameters** accurate for **consistency**.
- **Always verify settings** before starting a long run.
- Contact a **supervisor** once the production run **begins**.
- **Use nohup** for long-running jobs to prevent interruptions.

🛠 **Happy Simulating! 🚀**

