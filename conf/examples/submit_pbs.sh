#!/bin/bash
# Start rnaseq-splice on a PBS Pro or Torque/OpenPBS cluster. The job below is the Nextflow
# "head" job: it uses little CPU or memory but runs for the whole pipeline, and submits every
# task as its own job.
#
#   qsub submit_pbs.sh
#
# Edit the #PBS lines, the module line and the paths for your cluster.
# Torque/OpenPBS: replace the select line with  #PBS -l nodes=1:ppn=2,mem=8gb
#PBS -N rnaseq_splice
#PBS -q workq
#PBS -P my_project
#PBS -l select=1:ncpus=2:mem=8gb
#PBS -l walltime=168:00:00
#PBS -j oe

set -euo pipefail
cd "${PBS_O_WORKDIR}"

module load nextflow apptainer 2>/dev/null || true   # or however your cluster provides them

export NXF_OPTS='-Xms1g -Xmx4g'                                   # memory for Nextflow itself
export NXF_APPTAINER_CACHEDIR=/shared/containers/apptainer_cache   # shared container cache

nextflow run mattjarvis123/rnaseq-splice -r 1.0.0 \
    -profile apptainer \
    -c pbs_apptainer.config \
    -params-file params.yaml \
    -w /shared/scratch/$USER/rnaseq_splice/work \
    -resume
