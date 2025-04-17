#!/bin/bash
#SBATCH --export=ALL
#SBATCH --ntasks=8
#SBATCH --time=72:00:00
#SBATCH --mem-per-cpu=8000
#SBATCH --job-name=SComatic_QR_HCTbc_Analysis
#SBATCH --output=job_out_files/%x.%j.out

# Run SComatic on the QR_HCTbc 10X scRNA outputs.

# Activate conda environment.
conda activate SComatic

# Define paths.
SCOMATIC=/path/to/scripts/SComatic
dir_10X=/path/to/data/10X_scRNA
dir_SComatic=/path/to/data/SComatic_Outputs

# Move to 10X directory and process HCTbc directories.
cd ${dir_10X}
H_dirs=(*HCTbc*)

for dir in ${H_dirs[@]}; do 
    cd ${dir_10X}
    samp_name="$(echo ${dir} | sed -re 's/run_count_//')"
    cd ${dir}/outs/

    cd filtered_feature_bc_matrix
    gunzip -k barcodes.tsv.gz
    mv barcodes.tsv ${dir_SComatic}/QR_HCTbc/${samp_name}_barcodes.tsv

done

# Define sample names.
H_samp_names=("QR_HCTbc_POTs" "QR_HCTbc_CO3_P4_EP1" "QR_HCTbc_CO4_P4_EP1" "QR_HCTbc_DT3_P4_EP1" "QR_HCTbc_DT4_P4_EP1" "QR_HCTbc_DSI_P1_EP1" "QR_HCTbc_DSII_P1_EP1")

# Generate metadata file.
Rscript /path/to/scripts/SComatic/SComatic_QR_HCTbc_make_meta_files.r

# Step 1: Splitting alignment files into cell type-specific BAMs.
output_dir1=${dir_SComatic}/QR_HCTbc/Step1_BamCellTypes
mkdir -p $output_dir1

for samp in ${H_samp_names[@]}; do 
    python $SCOMATIC/scripts/SplitBam/SplitBamCellTypes.py --bam ${dir_SComatic}/QR_HCTbc/${samp}_possorted_genome_bam.bam \
        --meta ${dir_SComatic}/QR_HCTbc/${samp}_barcodes.tsv \
        --id QR_HCTbc \
        --n_trim 5 \
        --max_nM 5 \
        --max_NH 1 \
        --outdir $output_dir1
done

# Step 2: Collecting base count information.
REF=/path/to/refdata/genome.fa
output_dir2=${dir_SComatic}/QR_HCTbc/Step2_BaseCellCounts
mkdir -p $output_dir2

for samp in ${H_samp_names[@]}; do 
    samp_short="$(echo ${samp} | sed -re 's/QR_HCTbc_//')"
    bam=${output_dir1}/QR_HCTbc.${samp_short}.bam
    temp=${output_dir2}/temp_${samp}
    mkdir -p $temp

    python $SCOMATIC/scripts/BaseCellCounter/BaseCellCounter.py --bam $bam \
        --ref $REF \
        --chrom all \
        --out_folder $output_dir2 \
        --min_dp 10 \
        --min_cc 10 \
        --min_bq 30 \
        --tmp_dir $temp \
        --nprocs 8

    rm -rf $temp
done

# Step 3: Merging base count matrices.
output_dir3=${dir_SComatic}/QR_HCTbc/Step3_BaseCellCountsMerged
mkdir -p $output_dir3

python $SCOMATIC/scripts/MergeCounts/MergeBaseCellCounts.py --tsv_folder ${output_dir2} \
  --outfile ${output_dir3}/QR_HCTbc.BaseCellCounts.AllCellTypes.tsv

# Step 4: Detection of somatic mutations.
output_dir4=${dir_SComatic}/QR_HCTbc/Step4_VariantCalling
mkdir -p $output_dir4

REF=/path/to/refdata/genome.fa

python $SCOMATIC/scripts/BaseCellCalling/BaseCellCalling.step1.py \
          --infile ${output_dir3}/QR_HCTbc.BaseCellCounts.AllCellTypes.tsv \
          --outfile ${output_dir4}/QR_HCTbc \
          --ref $REF \
          --min_cov 10 \
          --min_ac_cells 3 \
          --min_ac_reads 4 \
          --max_cell_types 7 \
          --min_cell_types 4

editing=$SCOMATIC/RNAediting/AllEditingSites.hg38.txt
PON=$SCOMATIC/PoNs/PoN.scRNAseq.hg38.tsv

python $SCOMATIC/scripts/BaseCellCalling/BaseCellCalling.step2.py \
          --infile ${output_dir4}/QR_HCTbc.calling.step1.tsv \
          --outfile ${output_dir4}/QR_HCTbc \
          --editing $editing \
          --pon $PON

bedtools intersect -header -a ${output_dir4}/QR_HCTbc.calling.step2.tsv -b $SCOMATIC/bed_files_of_interest/UCSC.k100_umap.without.repeatmasker.bed | awk '$1 ~ /^#/ || $6 == "PASS"' > ${output_dir4}/QR_HCTbc.calling.step2.pass.tsv
