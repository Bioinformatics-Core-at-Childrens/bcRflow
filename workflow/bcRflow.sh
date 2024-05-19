#!/bin/bash
#SBATCH --job-name=<bcr-nextflow>
#SBATCH -N 1 # Ensure that all cores are on one machine
#SBATCH -t 3-00:00 #runtime in D-HH:MM
#SBATCH --cpus-per-task=4 # Request that ncpus be allocated per process.
#SBATCH --mem=16G
#SBATCH --output=./slurm_log/bcr-nf_%A.out
#SBATCH -e ./slurm_log/bcr-nf_%A.err
#SBATCH --mail-type=ALL
#SBATCH --mail-user=your_email@sample.com

unset TMPDIR

module purge
module load squashfs-tools/4.4 gcc/12.2.0 nextflow/23.04.2 singularity/3.9.6
export NXF_SINGULARITY_CACHEDIR=/path/to/bcRflow/singularity-images
export SINGULARITY_CACHEDIR=/path/to/bcRflow/singularity-images

cd /path/to/bcRflow/workflow
nextflow run ./main.nf -profile slurm -resume -work-dir ./work
