#!/bin/bash

#this is an array job - runs one simulation per job

#--------------- EXECUTABLE ---------------

#define unique slurm jobid
UNIQUEJOBID=${SLURM_ARRAY_JOB_ID}_"slimIter_"${SLURM_ARRAY_TASK_ID}
#define working dir on execute node
WORKINGDIR=/tmp/local/$SLURM_JOB_ID/$UNIQUEJOBID
#define prefix name for outputs
TREEFILE="/wmModel_"${UNIQUEJOBID}"_sigma_"$SIGMA"_K_"$K

#check if these directories have been created; if not, make them
if [ ! -d $WORKINGDIR ]; then mkdir $WORKINGDIR; fi
if [ ! -d $OUTDIR ]; then mkdir $OUTDIR; fi

#copy scripts to execute node
cp $SLURM_SUBMIT_DIR/test.slim $WORKINGDIR
cp $SLURM_SUBMIT_DIR/processing_slim.py $WORKINGDIR
cp $SLURM_SUBMIT_DIR/exe_WM.R $WORKINGDIR
cp $SLURM_SUBMIT_DIR/wmModel_plots.R $WORKINGDIR
cp $SLURM_SUBMIT_DIR/wm_lib.R $WORKINGDIR
cp $SLURM_SUBMIT_DIR/wm_hom_cmpPar_mod_block_scaled.R $WORKINGDIR
cp $SLURM_SUBMIT_DIR/wm_hom_cmpPar_cmpLnL_mod_block_scaled.R $WORKINGDIR

#move to execute node
cd $WORKINGDIR

### RUN SLIM ###

# load all of the programs that we want to use
# note - if defining variables with -d for slim, don't include a defineConstant() line in slim script
# text variables defined with -d should use syntax - "MYVAR='$MYVAR'"
module purge
module load GCC/6.4.0-2.28
module load OpenMPI/2.1.2
module load SLiM/2019dev

#format: month, day, year, hour (24 hr format), minutes
d=`date +%m,%d,%Y,%H,%M`
echo "TIME,START,SLIM,$d"

slim -d K=$K -d SIGMA=$SIGMA -d "WORKINGDIR='$WORKINGDIR'" -d "TREEFILE='$TREEFILE'" -d "n='$SLURM_ARRAY_TASK_ID'" test.slim

echo "DONE RUNNING SLIM"
d=`date +%m,%d,%Y,%H,%M`
echo "TIME,END,SLIM,$d"
#generates tree sequence for each sigma and K combo

### PARSE TREE FILES IN PYTHON ###

#module load GCCcore/11.1.0
#module load Python/3.9.6

export PATH=$PATH:$HOME/anaconda3/bin

source activate myconda

d=`date +%m,%d,%Y,%H,%M`
echo "TIME,START,PYTHON,$d"

python3 processing_slim.py $WORKINGDIR $TREEFILE

conda deactivate

echo "DONE RUNNING PYTHON SCRIPT"
d=`date +%m,%d,%Y,%H,%M`
echo "TIME,END,PYTHON,$d"

#generates:
#pi and geographic distance matrices

### RUN STAN MODEL IN R ###

# load all of the programs that we want to use
module purge
module load GCC/10.2.0  OpenMPI/4.0.5
module load foss/2020b
module load R/4.0.3

d=`date +%m,%d,%Y,%H,%M`
echo "TIME,START,EXE_WM,$d"

Rscript exe_WM.R $WORKINGDIR $TREEFILE
echo "DONE RUNNING exe_WM.R SCRIPT"
d=`date +%m,%d,%Y,%H,%M`
echo "TIME,END,EXE_WM,$d"
echo ""
#generates:
#-out.Robj
#-estplots.pdf
#for each sigma and K combo

d=`date +%m,%d,%Y,%H,%M`
echo "TIME,START,COLLECT_PLOTS,$d"

Rscript wmModel_plots.R $WORKINGDIR $TREEFILE
#echo "DONE RUNNING collectpi_plots.R SCRIPT"
#-est_values.txt
#-est_plot.pdf
#for each sigma / K combo, slimIter

d=`date +%m,%d,%Y,%H,%M`
echo "TIME,END,COLLECT_PLOTS,$d"

#copy all output files back to storage node
cp wmModel_* $OUTDIR
wait
rm -rf $WORKINGDIR

echo "ALL DONE WITH THE WHOLE JOB"


#print some environment variables to stdout for records
echo ----------------------------------------------------------------------------------------
echo PRINTING SUBSET OF ENVIRONMENT VARIABLES:
(set -o posix ; set | grep -v ^_ | grep -v ^EB | grep -v ^BASH | grep -v PATH | grep -v LS_COLORS)