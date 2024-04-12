#!/bin/bash
#SBATCH --job-name=PP_ALL_224
#SBATCH --partition=mapu
#SBATCH -N 1
#SBATCH --ntasks-per-node=32

module load mkl
module load compiler
module load hdf5/1.14.1-2_intel
module load Python/3.11.4
module load openmpi/4.1.5

export OMP_NUM_THREADS=32

cd /home/saguilar/Network/krome/build_dust/
python3 post_process.py /home/saguilar/MAGISTER_DATA/snapshots/ > /home/saguilar/Network/krome/build_dust/logfile
