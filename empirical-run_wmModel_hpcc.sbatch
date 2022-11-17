#!/bin/bash

# Shell script to run R script for WM gendiv extrapolation model

#load programs
module purge
module load foss/2020b
module load R/4.0.3

#set up work environment for job depending on if you want to load files from where they live on cluster or from tmp dir on execute node
if [ $COPY_FILES_TO_EXECUTE_NODE = yes ]
then
	
	#make unique working directory on execute node
	WORKDIR=/tmp/local/$SLURM_JOB_ID #working dir on (remote) execute node
	if [ ! -d $WORKDIR ]; then mkdir $WORKDIR; fi
	
	#copy input files to tmpdir to speed up run time
	cp $INDIR/samplenamekey.txt $WORKDIR
	cp $INDIR/"max_and_pw_dists"* $WORKDIR
	cp $INDIR/"popgenstats."$MINPROPINDIVSSCOREDIN* $WORKDIR
	cp $INDIR/"bpstats."$MINPROPINDIVSSCOREDIN* $WORKDIR
	cp $SLURM_SUBMIT_DIR/exe_WM-emperical.R $WORKDIR
	cp $SLURM_SUBMIT_DIR/wm_lib.R $WORKDIR
	cp $SLURM_SUBMIT_DIR/models/wm_hom_cmpPar_mod_block_scaled.R $WORKDIR
	cp $SLURM_SUBMIT_DIR/models/wm_hom_cmpPar_cmpLnL_mod_block.R $WORKDIR
	
	#move to execute node
	cd $WORKDIR
	
	#define input directory for R for this job
	INDIR=$WORKDIR
	
	#define output directory for R for this job
	OUTDIR_FINAL=$OUTDIR
	OUTDIR_SUFFIX=$(echo $OUTDIR | rev | cut -d "/" -f1 | rev)
	OUTDIR=$WORKDIR
		
fi

#print files currently on execute node
echo "files currently where I am working from are"
du -a

#run R, with the name of my R script
Rscript $WORKDIR/exe_WM-emperical.R $RUN_NAME $STORAGENODE $INDIR $OUTDIR $MINPROPINDIVSSCOREDIN $MODEL_FLAVOR
wait
echo "done with exe_WM-emperical.R"

#print files currently on execute node
echo "files currently where I am working from are"
du -a

#copy all output files back to the output dir on storage node for this dataset,
#and also put a copy of just the figure .pdfs into one master output dir on storage node
if [ $COPY_FILES_TO_EXECUTE_NODE = yes ]
then
	
	cp -r $OUTDIR/WMfit* $OUTDIR_FINAL
	#cp -r $OUTDIR/plots* $FIGDIR
	
	wait
	rm -rf $WORKDIR
fi

# END

#print some environment variables to stdout for records
echo ----------------------------------------------------------------------------------------
echo PRINTING SUBSET OF ENVIRONMENT VARIABLES:
(set -o posix ; set | grep -v ^_ | grep -v ^EB | grep -v ^BASH | grep -v PATH | grep -v LS_COLORS)
