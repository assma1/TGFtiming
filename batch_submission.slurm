#!/bin/bash

for i in {1..972}
do
JOB_NAME="TGF_${i}"
output="/data/cephfs/punim0864/TGFtiming/out/$JOB_NAME"
sbatch --job-name=$JOB_NAME --output=$output /data/cephfs/punim0864/TGFtiming/src/slurm/job_submission.slurm $i
done
