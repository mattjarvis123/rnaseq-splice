#!/bin/bash
# Start rnaseq-splice on a SLURM cluster. The job below is the Nextflow "head" job: it uses
# little CPU or memory but runs for the whole pipeline, and submits every task as its own job.
#
#   sbatch submit_slurm.sh
#
# Edit the #SBATCH lines, the module line and the paths for your cluster.
#SBATCH --job-name=rnaseq_splice
#SBATCH --partition=normal
#SBATCH --account=my_account
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --time=7-00:00:00
#SBATCH --output=rnaseq_splice_%j.log

set -euo pipefail

module load nextflow apptainer 2>/dev/null || true   # or however your cluster provides them

export NXF_OPTS='-Xms1g -Xmx4g'                                   # memory for Nextflow itself
export NXF_APPTAINER_CACHEDIR=/shared/containers/apptainer_cache   # shared container cache

nextflow run mattjarvis123/rnaseq-splice -r main \
    -profile apptainer \
    -c slurm_apptainer.config \
    -params-file params.yaml \
    -w /shared/scratch/$USER/rnaseq_splice/work \
    -resume
