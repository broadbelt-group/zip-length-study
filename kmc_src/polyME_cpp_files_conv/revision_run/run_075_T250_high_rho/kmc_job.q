#!/bin/bash
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -t 48:00:00
#SBATCH --mem=5G
#SBATCH -A TG-CTS120055
#SBATCH -p shared 
#SBATCH --job-name="param_fit"
#SBATCH --mail-user=shivanikozarekar2026@u.northwestern.edu
#SBATCH --mail-type=FAIL

# Load gcc version that has c++17
module reset
module load gcc/10.2.0

# Run program
./program
