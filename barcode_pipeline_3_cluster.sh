
# Cluster Barcodes Script
# This script clusters barcodes extracted using Bartender.

#$ -S /bin/bash
#$ -N Cluster_Barcodes
#$ -j y
#$ -cwd
#$ -l tmem=10G
#$ -l h_vmem=10G
#$ -l h_rt=24:00:00
#$ -t 1-69
#$ -tc 10

# Load necessary environment variables and dependencies
source /path/to/your/environment/source_file

# Argument 1: Sample list for extracted barcodes
samp_list=$1

# Use task ID to extract sample name from the sample list
samp=$(sed -n ${SGE_TASK_ID}'{p;q}' ${samp_list})

# Change directory to the chosen sample
test -d "${samp}" && cd "${samp}" || { echo "Sample directory not found: ${samp}"; exit 1; }

# Get current date and time
NOW=$(date '+%F_%H:%M:%S')

# Add barcode extraction header to sample summary file
echo "
#######################################################################
# Barcode Clustering Results:

# Barcode Clustering Performed: ${NOW}
# Extracted Clustering Files: 
${samp}_barcode.csv
${samp}_cluster.csv
${samp}_quality.csv
" >> "${samp}_bc_pipeline_summary.txt"

# Run Bartender's clustering function and save output to the summary file
bartender_single_com -f "${samp}_barcode.txt" -o "${samp}" -c 10 -d 2 -z 5 &>> "${samp}_bc_pipeline_summary.txt"

# Add a footer to indicate completion
echo "
######################################################################
" >> "${samp}_bc_pipeline_summary.txt"
