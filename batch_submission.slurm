#!/bin/bash

for i in {1..384}
do
JOB_NAME="TGF_${i}"
output="/data/gpfs/projects/punim0864/TGFtiming/out/$JOB_NAME"
sbatch --job-name=$JOB_NAME --output=$output /data/gpfs/projects/punim0864/TGFtiming/src/slurm/job_submission.slurm $i
done
