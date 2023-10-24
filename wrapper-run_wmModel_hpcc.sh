#!/bin/bash

#define variables:
storagenode=exp_1_out #path to main node where input files live

logfilesdir=logfiles_exp_10_wmModel #name of directory to create and then write log files to

datesuffix=$(date +%m-%d-%Y.%H)
outdir=wm_exp_1 #name of directory to create and write all outputs to
indir=$storagenode
finaldir=$storagenode/summary_files #where the models go to die

model_flavor=wishart #value wishart or cmplnl

#define some values to pass into slim
vector_of_K_values=( 5 )
vector_of_sigma_values=( 1 )

cpus=2 #number of CPUs to request/use per dataset
ram_per_cpu=5G #amount of RAM to request/use per CPU
time=60:00:00

#slurm variable key:
# %A = SLURM_ARRAY_JOB_ID
# %a = SLURM_ARRAY_TASK_ID

#---------------------------------------------------------

jobname=wm_Model #label for SLURM book-keeping
executable=run_wmModel_hpcc.sh #script to run

#check if logfiles directory has been created in submit dir yet; if not, make one
if [ ! -d ./$logfilesdir ]; then mkdir ./$logfilesdir; fi

#submit job to cluster
for sigma in "${vector_of_sigma_values[@]}"
do
	for K in "${vector_of_K_values[@]}"
	do
	jid_pi=$(sbatch --job-name=$jobname \
					--account=bradburd1 \
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

jobname=run-compileoutputs #label for SLURM book-keeping
executable=run_collate_all.sh #script to run

sbatch --job-name=$jobname \
					--export=CPUS=$cpus,STORAGENODE=$storagenode,OUTDIR=$outdir,LOGFILESDIR=$logfilesdir,FINALDIR=$finaldir \
					--cpus-per-task=$cpus \
					--mem-per-cpu=$ram_per_cpu \
					--output=./$logfilesdir/${jobname}_%A.out \
					--error=./$logfilesdir/${jobname}_%A.err \
					--dependency=afterany$(eval echo \$all_pis_jids) \
					--kill-on-invalid-dep=yes \
					--time=$time \
					$executable

#DONE


