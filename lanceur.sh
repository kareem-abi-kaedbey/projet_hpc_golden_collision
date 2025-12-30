#!/bin/bash
#OAR -l /nodes=40,walltime=02:30:00
#OAR -O resultat_n38_final.%jobid%.log
#OAR -E erreur_n38_final.%jobid%.log

cd $OAR_WORKDIR

mpicc -Wall -O3 -fopenmp -o mitm mitm.c

export OMP_NUM_THREADS=18

echo "=== Lancement FINAL N=38 sur 40 nœuds ==="
echo "Date début : $(date)"

mpiexec --hostfile $OAR_NODE_FILE --map-by ppr:1:node:PE=18 --bind-to core ./mitm --n 38 --C0 1f042c951edfef75 --C1 187db4f5ebf4bc01

echo "Date fin : $(date)"
