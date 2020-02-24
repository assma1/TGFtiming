#!/bin/bash

#SBATCH -p cloud,physical

#SBATCH --time=7-00:00:00

#SBATCH --nodes=1

#SBATCH --ntasks=1

#SBATCH --cpus-per-task=1

#SBATCH --mem=10240

#SBATCH --mail-user=a.smart1@student.unimelb.edu.au
#SBATCH --mail-type=TIME_LIMIT_50,TIME_LIMIT_90

i=$1

module purge

module load R/3.4.1-GCC-4.9.2

cd /data/cephfs/punim0864/TGFtiming/

Rscript --vanilla src/spartan_TGF_parameters.R $i