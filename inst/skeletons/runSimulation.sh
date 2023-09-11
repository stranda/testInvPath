#!/bin/bash
###
### This bash script runs testInvPath simulations in a slurm-based HPC environment
###


# Check if a parameter is provided
if [ "$#" -ne 1 ]; then
    echo "Usage: $0 <number_of_jobs>"
    exit 1
fi

# Get the current directory name to use as the job name
d=`pwd`
job_name=`basename ${d}`

#get the container name
container=${HOME}/bin/testinvpath.simg

# Generate the Slurm script and pipe it directly to sbatch
cat <<EOL | sbatch
#!/bin/bash

#SBATCH --partition=stdmemq
#SBATCH --job-name=$job_name
#SBATCH --array=1-$1
#SBATCH --time=48:00:00
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=7500	# 5.5gb memory per core
#SBATCH --export=ALL
#SBATCH --output="${d}/logs/slurm-%A_%a.out"

scratchroot=/scratch/stranda

echo $d

cd $d/src

module load singularity
singularity run ${container} Rscript runReps.R 1  ${SLURM_JOB_ID}${SLURM_ARRAY_TASK_ID}
mv reference*.csv ${d}/results/

EOL

echo "Job submitted to Slurm!"
