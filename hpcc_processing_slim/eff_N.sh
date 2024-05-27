#!/bin/bash

#--------------- EXECUTABLE ---------------

#define unique slurm jobid
UNIQUEJOBID=${SLURM_JOB_ID}
#define working dir on execute node
WORKINGDIR=compiled_output
#define prefix name for outputs
#TREEFILE="/wmModel_sigma_"$SIGMA"_K_"$K

#check if these directories have been created; if not, make them
if [ ! -d $WORKINGDIR ]; then mkdir $WORKINGDIR; fi
if [ ! -d $OUTDIR ]; then mkdir $OUTDIR; fi

#copy scripts to execute node
cp $SLURM_SUBMIT_DIR/effect_N.py $WORKINGDIR
cp /nfs/turbo/lsa-bradburd/Zach/WM_output/compiled_output/wmModel_*_sigma_"$SIGMA"_K_"$K" $WORKINGDIR
cp /nfs/turbo/lsa-bradburd/Zach/WM_output/compiled_output/wmModel_*_sigma_"$SIGMA"_K_"$K"_popsize $WORKINGDIR


#move to execute node
cd $WORKINGDIR

### RUN STAN MODEL IN R ###

#tell me what files are in WORKINGDIR
echo "printing working dir stuff"

pwd
ls

### PARSE TREE FILES IN PYTHON ###

#module load GCCcore/11.1.0
#module load Python/3.9.6

export PATH=$PATH:$HOME/miniconda3/bin

source activate myconda

d=`date +%m,%d,%Y,%H,%M`
echo "TIME,START,PYTHON,$d"

python3 effect_N.py

conda deactivate

echo "DONE RUNNING PYTHON SCRIPT"
d=`date +%m,%d,%Y,%H,%M`
echo "TIME,END,PYTHON,$d"

#copy all output files back to storage node
cp wmModel_* $OUTDIR
wait
rm -rf $WORKINGDIR

echo "ALL DONE WITH THE WHOLE JOB"


#print some environment variables to stdout for records
echo ----------------------------------------------------------------------------------------
echo PRINTING SUBSET OF ENVIRONMENT VARIABLES:
(set -o posix ; set | grep -v ^_ | grep -v ^EB | grep -v ^BASH | grep -v PATH | grep -v LS_COLORS)