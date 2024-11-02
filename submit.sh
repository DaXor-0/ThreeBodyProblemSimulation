#!/bin/bash
# #SBATCH -N 7
# #SBATCH -n 28
# #SBATCH -p debug
# # #SBATCH --gres=gpu:1
# # #SBATCH --exclusive
# # #SBATCH --time=00:05:00
# # #SBATCH --account=tra24_epicure


N_OF_BODIES=(4 16 64 256 1024 4096 16384)
GRID_SIZE=(200 400 800 1600 3200 6400 12800)

RUN=mpicc
N_PROC=2
#RUN=srun

index=1

iter=20000
TIMESTAMP=$(date +"%Y_%m_%d___%H:%M:%S")
filename="./$TIMESTAMP.csv"
gifname="./$TIMESTAMP.gif"

$RUN -n $N_PROC ~/ThreeBodyProblemSimulation/simulation.out ${N_OF_BODIES[$index]} $iter $filename

source ~/venv1/bin/activate
python ~/ThreeBodyProblemSimulation/plot.py $filename $gifname ${N_OF_BODIES[$index]}
