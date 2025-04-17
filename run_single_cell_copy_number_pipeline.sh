#!/bin/bash
#
#SBATCH --job-name=hmm_singularity_runner
#SBATCH --output=hmm_singularity_runner.txt
#
#SBATCH --time=4:00:00
#SBATCH --cpus-per-task=16
#SBATCH --mem-per-cpu=8042

#export PATH=$PWD/miniconda3/bin:$PATH
 #single_cell --helpw
 
 single_cell hmmcopy \
 --input_yaml hmm_inputs.yaml \
 --maxjobs 16 \
 --sentinel_only --submit local --loglevel DEBUG \
 --tmpdir temp --pipelinedir pipeline_allcells --out_dir output_allcells \
 --config_override '{"refdir": "path/to/DLP/refdata"}' \
 --library_id POTs_EP1