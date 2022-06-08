#!/bin/bash

#define variables:
storagenode=/mnt/scratch/hancockz #path to main node where input files live

logfilesdir=logfiles_wmModel #name of directory to create and then write log files to

datesuffix=$(date +%m-%d-%Y.%H)
outdir=$storagenode/allpimat_wmModel_outputs-$datesuffix #name of directory to create and write all outputs to

n_iterations=10 #number of iterations to run of each identical simulation

#define some values to pass into slim
vector_of_K_values=( 2.0 5.0 10.0 25.0 50.0 100.0)
vector_of_sigma_values=( 0.5 0.75 1.0 1.25 1.5 2.0)

cpus=2 #number of CPUs to request/use per dataset
ram_per_cpu=16G #amount of RAM to request/use per CPU
time=168:00:00

#slurm variable key:
# %A = SLURM_ARRAY_JOB_ID
# %a = SLURM_ARRAY_TASK_ID

#---------------------------------------------------------

jobname=run-collectpi #label for SLURM book-keeping
executable=gen_pimat_exe.sh #script to run

#check if logfiles directory has been created in submit dir yet; if not, make one
if [ ! -d ./$logfilesdir ]; then mkdir ./$logfilesdir; fi

#submit job to cluster
for sigma in "${vector_of_sigma_values[@]}"
do
	for K in "${vector_of_K_values[@]}"
	do
	jid_pi=$(sbatch --job-name=$jobname \
					--export=CPUS=$cpus,STORAGENODE=$storagenode,OUTDIR=$outdir,LOGFILESDIR=$logfilesdir,K=$K,SIGMA=$sigma \
					--cpus-per-task=$cpus \
					--mem-per-cpu=$ram_per_cpu \
					--array=[0-$((n_iterations-1))]%n_iterations \
					--output=./$logfilesdir/${jobname}_%A_slimIter_%a_sigma_${sigma}_K_${K}.out \
					--error=./$logfilesdir/${jobname}_%A_slimIter_%a_sigma_${sigma}_K_${K}.err \
					--time=$time \
					$executable \
					|awk -v "nID=$n_iterations" '{for (i=0; i<nID; i++) printf(":%s",$4"_"i)}')
	declare "all_pis_jids=${jid_pi}${all_pis_jids}"
	
	echo "submitted job has sigma value $sigma and K value $K; running $n_iterations slim iterations"
	echo ""
	done
done

echo "all_pis_jids is $all_pis_jids"
echo ""

#DONE


