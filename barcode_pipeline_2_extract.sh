
# Extract Barcodes Script
# This script extracts barcodes from merged zipped FASTQ files.

#$ -S /bin/bash
#$ -N Extract_Barcodes
#$ -j y
#$ -cwd
#$ -l tmem=5G
#$ -l h_vmem=5G
#$ -l h_rt=24:00:00
#$ -t 1-3
#$ -tc 5

# Load necessary environment variables and dependencies
source /path/to/your/environment/source_file

# Argument 1: Sample list for post-merged FASTQ files
samp_list=$1

# Use task ID to extract sample name from the sample list
samp=$(sed -n ${SGE_TASK_ID}'{p;q}' ${samp_list})

# Change directory to the chosen sample
test -d "${samp}" && cd "${samp}" || { echo "Sample directory not found: ${samp}"; exit 1; }

# Unzip the FASTQ file
gunzip -c "${samp}_merged.fastq.gz" > "${samp}.fastq"

# Get current date and time
NOW=$(date '+%F_%H:%M:%S')

# Add barcode extraction header to sample summary file
echo "
#######################################################################
# Barcode Extractor Results:

# Barcode Extraction Performed: ${NOW}
# Extracted Barcode File: 
${samp}_extracted_barcode.txt
" >> "${samp}_bc_pipeline_summary.txt"

# Run Bartender's extractor function and save output to the summary file
bartender_extractor_com -f "${samp}.fastq" -o "${samp}" -q ? -p GACAG[30]AGCAG -m 2 -d b &>> "${samp}_bc_pipeline_summary.txt"

# Add a footer to indicate completion
echo "
######################################################################
" >> "${samp}_bc_pipeline_summary.txt"
