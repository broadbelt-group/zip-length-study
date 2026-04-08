#!/bin/bash
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -t 48:00:00
#SBATCH --mem=50G
#SBATCH -A TG-CTS120055
#SBATCH -p shared 
#SBATCH --job-name="DPn_100"
#SBATCH --mail-user=shivanikozarekar2026@u.northwestern.edu
#SBATCH --mail-type=FAIL
#SBATCH --mail-type=BEGIN
#SBATCH --mail-type=END

# Load gcc version that has c++17
module reset
module load gcc/10.2.0
module load hdf5/1.10.7

# Run program
./program
