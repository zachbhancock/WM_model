#!/bin/bash

#define variables:
storagenode=compiled_output #path to main node where input files live

logfilesdir=log_true_s #name of directory to create and then write log files to

datesuffix=$(date +%m-%d-%Y.%H)
outdir=true_s_out #name of directory to create and write all outputs to
indir=$storagenode
finaldir=$storagenode/summary_files #where the models go to die

model_flavor=wishart #value wishart or cmplnl

#define some values to pass into slim
vector_of_K_values=( 2.0 5.0 10.0 25.0 )
vector_of_sigma_values=( 0.5 0.75 1.0 1.25 1.5 2.0 )

cpus=2 #number of CPUs to request/use per dataset
ram_per_cpu=16G #amount of RAM to request/use per CPU
time=5:00:00

#slurm variable key:
# %A = SLURM_ARRAY_JOB_ID
# %a = SLURM_ARRAY_TASK_ID

#---------------------------------------------------------

jobname=wmModel #label for SLURM book-keeping
executable=run_true_s.sh #script to run

#check if logfiles directory has been created in submit dir yet; if not, make one
if [ ! -d ./$logfilesdir ]; then mkdir ./$logfilesdir; fi

#submit job to cluster
for sigma in "${vector_of_sigma_values[@]}"
do
	for K in "${vector_of_K_values[@]}"
	do
	jid_pi=$(sbatch --job-name=$jobname \
					--export=CPUS=$cpus,STORAGENODE=$storagenode,OUTDIR=$outdir,INDIR=$indir,LOGFILESDIR=$logfilesdir,K=$K,SIGMA=$sigma,MODEL_FLAVOR=$model_flavor \
					--cpus-per-task=$cpus \
					--mem-per-cpu=$ram_per_cpu \
					--output=./$logfilesdir/${jobname}_%A_sigma_${sigma}_K_${K}.out \
					--error=./$logfilesdir/${jobname}_%A_sigma_${sigma}_K_${K}.err \
					--time=$time \
					$executable \
					|cut -d " " -f 4)
	jid_pi=:${jid_pi}
	echo "submitted job has sigma $sigma and K $K and jid_pi $jid_pi"
	echo ""
	declare "all_pis_jids=${jid_pi}${all_pis_jids}"
	done
	
done

echo "DONE SUBMITTING ALL $jobname jobs"
echo "all_pis_jids is $all_pis_jids"


#---------------------------------------------------------
#submit job to cluster after all above jobs finish that compiles all final outputs into one summary .txt file

#note - dependency afterany says wait for all jobs to finish and then run (doesn't matter if jobs have exit status of 0 or not)

#DONE

