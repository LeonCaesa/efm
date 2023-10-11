#!/bin/bash -l

#$ -P dmfgrp
#$ -l h_rt=72:00:00
#$ -j y


echo "=========================================================="
echo "Starting on       : $(date)"
echo "Running on node   : $(hostname)"
echo "Current directory : $(pwd)"
echo "Current job ID    : $JOB_ID"
echo "Current job name  : $JOB_NAME"
echo "Task index number : $TASK_ID"
echo "=========================================================="



module load R/4.1.2
exp_idx=4
n_repeats=$SGE_TASK_ID
#n_repeats=1

for (( d=16; d<=500; d+= 50)); do
    echo "Running experiments for:  (exp_idx, repeat_idex, d) = (${exp_idx}, ${n_repeats}, ${d})"
    R -q --slave --vanilla --args $exp_idx $n_repeats $d < /projectnb/dmfgrp/efm/R/covexp.R
done


#qsub -t 1-5 submit_cov.sh


echo "=========================================================="
echo "Finished on       : $(date)"
echo "=========================================================="
