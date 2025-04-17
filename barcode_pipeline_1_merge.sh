
# Merge Forward and Reverse Reads Script
# This script merges forward and reverse reads from a NovaSeq run.

#$ -S /bin/bash
#$ -N Merge_Barcode_Reads
#$ -j y
#$ -cwd
#$ -l tmem=5G
#$ -l h_vmem=5G
#$ -l h_rt=24:00:00
#$ -t 1-14 
#$ -tc 5

# Load necessary environment variables and dependencies
source /path/to/your/environment/source_file

# Argument 1: Sample list for pre-merged FASTQ files
samp_list=$1

# Argument 2: Directory to store merged files
bc_dir=$2

# Use task ID to extract sample name from the sample list
samp=$(sed -n ${SGE_TASK_ID}'{p;q}' ${samp_list})

# Change directory to the chosen sample
test -d "${samp}" && cd "${samp}" || { echo "Sample directory not found: ${samp}"; exit 1; }

# Extract sample name from the first FASTQ file
files=( *.fastq.gz )
file=${files[0]}
samp_name=$(echo ${file} | sed -r 's/[A-Za-z0-9]+_[0-9]+_(.*)_S[0-9]+_R[0-9]+_001.fastq.gz/\1/')

# Create a directory for merged files
new_samp_dir="${bc_dir}${samp_name}/"
mkdir -p "${new_samp_dir}"

# Run NGmerge to merge reads and save outputs in the new directory
NGmerge -1 "${files[0]}" -2 "${files[1]}" -o "${new_samp_dir}${samp_name}_merged.fastq" -f "${new_samp_dir}${samp_name}_mFAIL" -l "${new_samp_dir}${samp_name}_NGmerge_log" -m 14 -p 0.2 -z

# Get current date and time
NOW=$(date '+%F_%H:%M:%S')

# Calculate number of reads pre- and post-merging
pre_reads=$(zcat "${files[0]}" | echo $((`wc -l` / 4)))
post_reads=$(zcat "${new_samp_dir}${samp_name}_merged.fastq.gz" | echo $((`wc -l` / 4)))
fail_reads=$(zcat "${new_samp_dir}${samp_name}_mFAIL_1.fastq.gz" | echo $((`wc -l` / 4)))

# Generate a sample summary file
echo "
#######################################################################
######## Barcode Analysis Pipeline - Sample Summary File ##############
########################## Summary ###############################

# Sample File Created: ${NOW}
# Sequencing Run Directory: ${bc_dir}
# Sample Name: ${samp_name}

#######################################################################
# NGMerge Results:
##################

# Reads pre-merging: ${pre_reads}
# Reads post-merging: ${post_reads}
# Reads failed-merging: ${fail_reads}

#######################################################################
" >> "${new_samp_dir}${samp_name}_bc_pipeline_summary.txt"
