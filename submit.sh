#!/bin/bash
#SBATCH -N 1
#SBATCH -n 32
#SBATCH -p boost_usr_prod
#SBATCH --gres=gpu:1
# #SBATCH --exclusive
#SBATCH --time=00:05:00
#SBATCH --account=tra24_epicure

module purge
module load openmpi

n_of_bodies=128
iter=20000
TIMESTAMP=$(date +"%Y_%m_%d___%H:%M:%S")
filename="./$TIMESTAMP.csv"
gifname="./$TIMESTAMP.gif"

srun -n 32 ~/ThreeBodyProblemSimulation/simulation.out $n_of_bodies $iter $filename

module load python
source ~/myvenv/bin/activate
python ~/ThreeBodyProblemSimulation/plot.py $filename $gifname
deactivate
