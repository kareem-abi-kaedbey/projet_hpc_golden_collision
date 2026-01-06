# PPAR Project: Parallel MITM Attack

**Authors:** [Victor] [ZHOU] & [Kareem] [ABI KAEDBEY]
**Date:** January 2026

---

## 1. Project Description
Our primary objective was to find the golden collision for the highest possible value of N. We developed a distributed hybrid architecture combining MPI and OpenMP. Our implementation strategy focuses on two critical performance factors: minimizing the memory footprint per node and controlling network congestion.

The project has been optimized and validated to solve instances up to **N=38** using 40 nodes on the cluster.

### Included Files
* `mitm.c`: Main source code (C / MPI / OpenMP).
* `Makefile`: Automation script for compilation.
* `lanceur.sh`: Batch submission script for OAR (configured for N=38).
* `rapport.pdf`: Detailed technical report.

---

## 2. Compilation

The project provides a standard `Makefile`.

```bash
make
```

This will generate the binary executable **`mitm`**.

> **Cleaning:** To remove the executable and object files:
> ```bash
> make clean
> ```

---

## 3. Execution: Interactive Mode (Quick Tests)

This mode is recommended for verifying the program on small instances (e.g., N=35 or N=36) without waiting for a long reservation queue.

### Step 1: Resource Reservation
Reserve nodes (e.g., 4 nodes for 30 minutes) in interactive mode:

```bash
oarsub -I -l nodes=4,walltime=00:30:00
```

### Step 2: Launch Calculation
Once connected to the master node, execute the following commands.
**Important:** The MPI mapping (`--map-by ppr:1:node:PE=18`) is crucial to ensure each MPI process utilizes all 18 cores via OpenMP.

```bash

#1 Set the number of OpenMP threads (18 cores on 'gros' cluster)
export OMP_NUM_THREADS=18

#2 Recompile with Makefile Don't forget to adapt the number of Chunks : NB_CHUNKS inside mitm.c
make

#3 Run the program (Example for N=35)
mpiexec --hostfile $OAR_NODE_FILE --map-by ppr:1:node:PE=18 --bind-to core \
    ./mitm --n 35 --C0 539d6089e2b1db91 --C1 2ff8b5c9b980f116
```


## 4. Execution: Batch Mode (Reproducing N=38)

This mode is used for intensive calculations requiring the allocation of many nodes (40+) and fine-grained memory management.

We provide the **`lanceur.sh`** script configured for N=38.

### Step 1: Check or Modify the Script
If you wish to change the node count or difficulty, edit the script:
```bash
nano lanceur.sh
```
* To modify resources: `#OAR -l /nodes=40...`
* To modify program parameters: the `mpiexec ...` line.

### Step 2: Submit Job
Ensure the script is executable, then submit it to the job scheduler (OAR).
**Note:** `chmod +x` is important to avoid "command not found" errors.

```bash
chmod +x lanceur.sh
oarsub -S ./lanceur.sh
```
The results (standard output) will be written to the file `resultat_n38_2_final.<JOB_ID>.log`.

## 5. Technical Notes & Configuration

### Network Optimization (Chunking)
To avoid saturating the Infiniband network during global exchanges (`MPI_Alltoall`), the code divides the search space into subsets (Chunks) processed sequentially.
* This value (`NUM_CHUNKS`) is defined in `mitm.c` and must be adapt with the value of n.
* It is calibrated to ensure packets remain within the optimal network performance zone (< 200 MB per exchange).


## 6. IMPORTANT NOTE:

We successfully validated N=38 on 40 nodes during our first test run.
However, a subsequent attempt with the exact same configuration failed.
To ensure a stable and successful execution for grading, we highly recommend
using 50 nodes instead of 40.