#!/bin/bash
###
### This bash script runs testInvPath ABC analysis of simulations in a slurm-based HPC environment
###


# Check if a parameter is provided
if [ "$#" -ne 1 ]; then
    echo "Usage: $0 <number_of_cores on bigmem001>"
    exit 1
fi

# Get the current directory name to use as the job name
d=`pwd`
job_name=`basename ${d}`
job_name=${job_name}ABC

#get the container name
container=${HOME}/bin/testinvpath.simg

# Generate the Slurm script and pipe it directly to sbatch
cat <<EOL | sbatch
#!/bin/bash

#SBATCH --partition=bigmemq
#SBATCH --job-name=$job_name
#SBATCH --time=48:00:00
#SBATCH --cpus-per-task=$1
#SBATCH --mem-per-cpu=200000	# 5.5gb memory per core
#SBATCH --export=ALL
#SBATCH --output="${d}/logs/ABClog-%A_%a.out"

scratchroot=/scratch/stranda

echo $d

cd $d/src

module load singularity
singularity run ${container} Rscript abc.R $1
mv reference*.csv ${d}/results/

EOL

echo "Job submitted to Slurm!"
