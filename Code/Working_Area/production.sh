#!/bin/bash
#SBATCH -t 02:30:00 # execution time. Ex: 1 hour
#SBATCH -n 1 -c 32 # number of tasks, number of cores
#SBATCH --ntasks-per-node=1
#SBATCH --gres=gpu:a100:1
#SBATCH --mem=4G

module load cesga/2022 gcc/system openmpi/4.1.4 gromacs/2024.2-PLUMED-2.9.2

#srun gmx_mpi mdrun -pin on -cpi -noappend -s ../Output/prod.tpr
srun gmx_mpi mdrun -pin on -cpi ../Output/state.cpt -s ../Output/prod.tpr -o ../Output/traj.xtc -g ../Output/md.log -c ../Output/confout.gro #-cpi ../Output/state.cpt