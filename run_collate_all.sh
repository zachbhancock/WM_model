#!/bin/bash

#--------------- EXECUTABLE ---------------

#define unique slurm jobid
UNIQUEJOBID=${SLURM_JOB_ID}
#define working dir on execute node
WORKINGDIR=/slim_output
#copy scripts to execute node
cp $SLURM_SUBMIT_DIR/collate_all_txt_outputs.R $WORKINGDIR

#check if these directories have been created; if not, make them
if [ ! -d $WORKINGDIR ]; then mkdir $WORKINGDIR; fi
if [ ! -d $OUTDIR ]; then mkdir $OUTDIR; fi
if [ ! -d $FINALDIR ]; then mkdir $FINALDIR; fi

cp "$OUTDIR"/wmModel_*-est.txt $WORKINGDIR

#move to execute node
cd $WORKINGDIR

### RUN STAN MODEL IN R ###

# load all of the programs that we want to use
#module purge
#module load GCC/10.2.0  OpenMPI/4.0.5
#module load foss/2020b
module load R/4.0.3

#rscript stuff
Rscript collate_all_txt_outputs.R $WORKINGDIR
d=`date +%m,%d,%Y,%H,%M`
echo "TIME,END,COLLECT_PLOTS,$d"

cp alloutput* $FINALDIR
wait
rm -rf $WORKINGDIR

echo "ALL DONE WITH THE WHOLE JOB"


#print some environment variables to stdout for records
echo ----------------------------------------------------------------------------------------
echo PRINTING SUBSET OF ENVIRONMENT VARIABLES:
(set -o posix ; set | grep -v ^_ | grep -v ^EB | grep -v ^BASH | grep -v PATH | grep -v LS_COLORS)