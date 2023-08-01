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
#cp $SLURM_SUBMIT_DIR/run_final_hpcc.r $WORKINGDIR
#cp $SLURM_SUBMIT_DIR/wmModel_plots.R $WORKINGDIR
#cp $SLURM_SUBMIT_DIR/slim_output/wmModel_*_sigma_"$SIGMA"_K_"$K"-pi.csv $WORKINGDIR
#cp $SLURM_SUBMIT_DIR/slim_output/wmModel_*_sigma_"$SIGMA"_K_"$K"-pi_locs.txt $WORKINGDIR


#move to execute node
cd $WORKINGDIR

### RUN STAN MODEL IN R ###

#tell me what files are in WORKINGDIR
echo "printing working dir stuff"

pwd
ls

# load all of the programs that we want to use
#module purge
#module load GCC/10.2.0  OpenMPI/4.0.5
#module load foss/2020b
module load R/4.2.0

cp wmModel_* $OUTDIR

d=`date +%m,%d,%Y,%H,%M`
echo "TIME,START,COLLECT_PLOTS,$d"

Rscript run_final_hpcc.r $WORKINGDIR
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