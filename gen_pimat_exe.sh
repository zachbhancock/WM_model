#!/bin/bash

#this is an array job - runs one simulation per job

#--------------- EXECUTABLE ---------------

#define unique slurm jobid
UNIQUEJOBID=${SLURM_ARRAY_JOB_ID}_"slimIter_"${SLURM_ARRAY_TASK_ID}
#define working dir on execute node
WORKINGDIR=/home/hancockz/$SLURM_JOB_ID/$UNIQUEJOBID
#define prefix name for outputs
TREEFILE="/wmModel_"${UNIQUEJOBID}"_sigma_"$SIGMA"_K_"$K

#check if these directories have been created; if not, make them
if [ ! -d $WORKINGDIR ]; then mkdir $WORKINGDIR; fi
if [ ! -d $OUTDIR ]; then mkdir $OUTDIR; fi

#copy scripts to execute node
cp $SLURM_SUBMIT_DIR/test.slim $WORKINGDIR
cp $SLURM_SUBMIT_DIR/processing_slim.py $WORKINGDIR

#move to execute node
cd $WORKINGDIR

### RUN SLIM ###

# load all of the programs that we want to use
# note - if defining variables with -d for slim, don't include a defineConstant() line in slim script
# text variables defined with -d should use syntax - "MYVAR='$MYVAR'"
#module purge
#module load GCC/6.4.0-2.28  OpenMPI/2.1.2
#module load GCC/6.4.0-2.28  OpenMPI/2.1.2-CUDA
#module load GNU/6.4.0-2.28  OpenMPI/2.1.2
#module load GNU/6.4.0-2.28  OpenMPI/2.1.2-CUDA
#module load SLiM/3.7.1

#format: month, day, year, hour (24 hr format), minutes
d=`date +%m,%d,%Y,%H,%M`
echo "TIME,START,SLIM,$d"

export PATH=$PATH:$HOME/miniconda3/bin

source activate myconda

slim -d K=$K -d SIGMA=$SIGMA -d "WORKINGDIR='$WORKINGDIR'" -d "TREEFILE='$TREEFILE'" -d "n='$SLURM_ARRAY_TASK_ID'" test.slim

echo "DONE RUNNING SLIM"
d=`date +%m,%d,%Y,%H,%M`
echo "TIME,END,SLIM,$d"
#generates tree sequence for each sigma and K combo

echo these are the files that I currently have
du -a

### PARSE TREE FILES IN PYTHON ###

#module load GCCcore/11.1.0
#module load Python/3.9.6

d=`date +%m,%d,%Y,%H,%M`
echo "TIME,START,PYTHON,$d"

python3 processing_slim.py $WORKINGDIR $TREEFILE

conda deactivate

echo "DONE RUNNING PYTHON SCRIPT"
d=`date +%m,%d,%Y,%H,%M`
echo "TIME,END,PYTHON,$d"

#generates:
#pi and geographic distance matrices

#copy all output files back to storage node
cp wmModel_* $OUTDIR
wait
rm -rf $WORKINGDIR

echo "ALL DONE WITH THE WHOLE JOB"


#print some environment variables to stdout for records
echo ----------------------------------------------------------------------------------------
echo PRINTING SUBSET OF ENVIRONMENT VARIABLES:
(set -o posix ; set | grep -v ^_ | grep -v ^EB | grep -v ^BASH | grep -v PATH | grep -v LS_COLORS)