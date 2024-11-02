#!/bin/bash
#SBATCH --nodes=7                  # Total number of nodes
#SBATCH --ntasks-per-node=7        # Number of tasks per node
#SBATCH --cpus-per-task=4          # Number of OpenMP threads per MPI task
#SBATCH -p debug
#SBATCH --exclusive
# # #SBATCH --time=00:05:00
# # #SBATCH --account=tra24_epicure

RUN=srun

N_OF_BODIES=(4 16 64 256 1024 4096 16384)
GRID_SIZE=(200 400 800 1600 3200 6400 12800)

TIMESTAMP=$(date +"%Y_%m_%d___%H:%M:%S")
filename="./$TIMESTAMP.csv"
gifname="./$TIMESTAMP.gif"

index=1
iter=20000

export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK 
$RUN -n $SLURM_NTASKS ~/ThreeBodyProblemSimulation/simulation.out ${N_OF_BODIES[$index]} $iter $filename

source ~/venv1/bin/activate
python ~/ThreeBodyProblemSimulation/plot.py $filename $gifname ${GRID_SIZE[$index]}
