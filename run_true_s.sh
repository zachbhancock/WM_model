#!/bin/bash

#this is an array job - runs one simulation per job

#--------------- EXECUTABLE ---------------

#define unique slurm jobid
UNIQUEJOBID=${SLURM_ARRAY_JOB_ID}_"slimIter_"${SLURM_ARRAY_TASK_ID}
#define working dir on execute node
WORKINGDIR=$SLURM_JOB_ID"_"$UNIQUEJOBID
#define prefix name for outputs
TREEFILE="/wmModel_"${UNIQUEJOBID}"_sigma_"$SIGMA"_K_"$K

#check if these directories have been created; if not, make them
if [ ! -d $WORKINGDIR ]; then mkdir $WORKINGDIR; fi
if [ ! -d $OUTDIR ]; then mkdir $OUTDIR; fi

### PARSE TREE FILES IN PYTHON ###

#module load GCCcore/11.1.0
#module load Python/3.9.6

d=`date +%m,%d,%Y,%H,%M`
echo "TIME,START,PYTHON,$d"

export PATH=$PATH:$HOME/miniconda3/bin

source activate myconda

python3 true_s.py $WORKINGDIR $TREEFILE

conda deactivate

echo "DONE RUNNING PYTHON SCRIPT"
d=`date +%m,%d,%Y,%H,%M`
echo "TIME,END,PYTHON,$d"

#generates:
#pi and geographic distance matrices

#copy all output files back to storage node
cp wmModel_* $OUTDIR
cp wmModel_* $WORKINGDIR
wait
#rm -rf $WORKINGDIR

echo "ALL DONE WITH THE WHOLE JOB"


#print some environment variables to stdout for records
echo ----------------------------------------------------------------------------------------
echo PRINTING SUBSET OF ENVIRONMENT VARIABLES:
(set -o posix ; set | grep -v ^_ | grep -v ^EB | grep -v ^BASH | grep -v PATH | grep -v LS_COLORS)