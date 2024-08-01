#!/bin/bash
#SBATCH --job-name=<bcr-covid>
#SBATCH -N 1 # Ensure that all cores are on one machine
#SBATCH -t 3-00:00 #runtime in D-HH:MM
#SBATCH --cpus-per-task=4 # Request that ncpus be allocated per process.
#SBATCH --mem=16G
#SBATCH --output=./slurm_log/bcRflow.out
#SBATCH -e ./slurm_log/bcRflow.err
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=user@email.com

unset TMPDIR

module purge
module load squashfs-tools/4.4 gcc/12.2.0 nextflow/23.04.2 singularity/3.9.6

cd ./bcRflow/workflow
nextflow run ./main.nf -profile htc -resume -work-dir /path/to/long_covid/bcRflow-work
