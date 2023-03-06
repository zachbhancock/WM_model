#!/bin/bash

#--------------- EXECUTABLE ---------------

#define unique slurm jobid
UNIQUEJOBID=${SLURM_JOB_ID}
#define working dir on execute node
WORKINGDIR=/tmp/local/$SLURM_JOB_ID
#define prefix name for outputs
#TREEFILE="/wmModel_sigma_"$SIGMA"_K_"$K

#check if these directories have been created; if not, make them
if [ ! -d $WORKINGDIR ]; then mkdir $WORKINGDIR; fi
if [ ! -d $OUTDIR ]; then mkdir $OUTDIR; fi

#copy scripts to execute node
cp $SLURM_SUBMIT_DIR/est_col_edges.R $WORKINGDIR
cp "$INDIR"/wmModel_*_sigma_"$SIGMA"_K_"$K"-pi.csv $WORKINGDIR
cp "$INDIR"/wmModel_*_sigma_"$SIGMA"_K_"$K"-pi_locs.txt $WORKINGDIR
cp "$INDIR"/wmModel*_sigma_"SIGMA"_K_"$K"_out.Robj $WORKINGDIR


#move to execute node
cd $WORKINGDIR

### RUN STAN MODEL IN R ###

# load all of the programs that we want to use
module purge
module load GCC/10.2.0  OpenMPI/4.0.5
module load foss/2020b
module load R/4.0.3

d=`date +%m,%d,%Y,%H,%M`
echo "TIME,START,EXE_WM,$d"

#Rscript exe_WM.R $WORKINGDIR $TREEFILE $MODEL_FLAVOR
Rscript est.col.edges.R $WORKINGDIR $MODEL_FLAVOR
echo "DONE RUNNING exe_WM.R SCRIPT"
d=`date +%m,%d,%Y,%H,%M`
echo "TIME,END,EXE_WM,$d"
echo ""
#generates:
#-out.Robj
#-estplots.pdf
#for each sigma and K combo

cp wmModel_* $OUTDIR

d=`date +%m,%d,%Y,%H,%M`
echo "TIME,START,COLLECT_PLOTS,$d"

#copy all output files back to storage node
cp wmModel_* $OUTDIR
wait
rm -rf $WORKINGDIR

echo "ALL DONE WITH THE WHOLE JOB"


#print some environment variables to stdout for records
echo ----------------------------------------------------------------------------------------
echo PRINTING SUBSET OF ENVIRONMENT VARIABLES:
(set -o posix ; set | grep -v ^_ | grep -v ^EB | grep -v ^BASH | grep -v PATH | grep -v LS_COLORS)