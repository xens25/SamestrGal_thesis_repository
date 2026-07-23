#!/bin/bash

# This script shows the command sequence and exact parameters needed to run each SameStr tool directly from 
# the command line, to replicate the analysis without using Galaxy. The parameters match those used in 
# the original SameStr publication.


# ============================================================
# 1. Download and prepare reference database
# ============================================================


# MetaPhlAn database (mpa_vJan25_CHOCOPhlAnSGB_202503)
metaphlan --install --index mpa_vJan25_CHOCOPhlAnSGB_202503 --db_dir db_MetaPhlAn/

# mOTUs database (v3.1.0)

wget https://zenodo.org/records/7778108/files/db_mOTU_v3.1.0.tar.gz
tar -xzvf db_mOTU_v3.1.0.tar.gz


# ============================================================
# 2. Create SameStr reference database
# ============================================================


# For MetaPhlAn markers
samestr db \
    --markers-info db_MetaPhlAn/mpa_vJan25_CHOCOPhlAnSGB_202503.pkl \
    --markers-fasta db_MetaPhlAn/mpa_vJan25_CHOCOPhlAnSGB_202503.fna.bz2 \
    --db-version db_MetaPhlAn/mpa_latest \
    --output-dir samestr_db/


# For mOTUs markers
samestr db \
    --markers-info db_mOTU/db_mOTU_taxonomy_ref-mOTUs.tsv db_mOTU/db_mOTU_taxonomy_meta-mOTUs.tsv \
    --markers-fasta db_mOTU/db_mOTU_DB_CEN.fasta \
    --db-version db_mOTU/db_mOTU_versions \
    --output-dir samestr_db/


# If these files are not created automatically during database setup, mpa_latest must be a text
# file containing the exact version of the MetaPhlAn database, and db_mOTU_versions the exact version of the mOTUs database.


# ============================================================
# 3. Pre-processing with KneadData
# ============================================================

# Download human database for host decontamination

kneaddata_database \
    --download human_genome bowtie2 \
    human_genome/


# ${ID} = sample identifier

kneaddata \
    -i1 ${ID}.R1.fastq.gz \
    -i2 ${ID}.R2.fastq.gz \
    -db human_genome/ \
    --trimmomatic-options 'SLIDINGWINDOW:4:20 MINLEN:70' \
    --bypass-trf \
    --quality-scores phred33 \
    --output-prefix ${ID} \
    -o out_kneaddata/

# Combine paired and unmatched reads into a single R1 and single R2 file

cat out_kneaddata/${ID}_kneaddata_paired_1.fastq out_kneaddata/${ID}_kneaddata_unmatched_1.fastq \
    > out_kneaddata/${ID}.R1.fastq
cat out_kneaddata/${ID}_kneaddata_paired_2.fastq out_kneaddata/${ID}_kneaddata_unmatched_2.fastq \
    > out_kneaddata/${ID}.R2.fastq

gzip out_kneaddata/${ID}.R1.fastq out_kneaddata/${ID}.R2.fastq


# ============================================================
# 4. Taxonomic profiling and alignment
# ============================================================

# For MetaPhlAn

metaphlan out_kneaddata/${ID}.R1.fastq.gz,out_kneaddata/${ID}.R2.fastq.gz \
    --input_type fastq \
    --db_dir db_MetaPhlAn/ \
    --index mpa_vJan25_CHOCOPhlAnSGB_202503 \
    --nproc 30 \
    -t rel_ab \
    --mapout out_align/${ID}.mapout \
    --samout out_align/${ID}.sam.bz2 \
    -o out_align/${ID}.profile.txt

# Relevant outputs for subsequent tools are: sam.bz2 and .profile.txt

# For mOTUs
motus profile \
    -f out_kneaddata/${ID}.R1.fastq.gz -r out_kneaddata/${ID}.R2.fastq.gz \
    -g 1 \
    -t 30 \
    -y insert.raw_counts \
    -o out_align/${ID}.profile.txt

motus map_snv \
    -f out_kneaddata/${ID}.R1.fastq.gz -r out_kneaddata/${ID}.R2.fastq.gz \
    -t 30 \
    -o out_align/${ID}.bam

# ============================================================
# 5. SameStr Tools
# ============================================================

# SameStr Convert 


# Use *.sam.bz2 if profiling was done with MetaPhlAn and *.bam if using mOTUs

samestr convert \
    --input-files out_align/*.sam.bz2 \
    --marker-dir samestr_db/ \
    --nprocs 30 \
    --min-vcov 3 \
    --min-aln-identity 0.9 \
    --min-aln-len 40 \
    --min-base-qual 20 \
    --min-aln-qual 0 \
    --output-dir out_convert/

# SameStr Merge

samestr merge \
    --input-files out_convert/*.npz \
    --marker-dir samestr_db/ \
    --nprocs 30 \
    --output-dir out_merge/

# SameStr Filter

samestr filter \
    --input-files out_merge/*.npz \
    --input-names out_merge/*.names.txt \
    --marker-dir samestr_db/ \
    --clade-min-n-hcov 5000 \
    --clade-min-samples 2 \
    --marker-trunc-len 20 \
    --global-pos-min-n-vcov 2 \
    --sample-pos-min-n-vcov 1 \
    --sample-pos-min-sd-vcov 3.0 \
    --sample-var-min-n-vcov 2 \
    --sample-var-min-f-vcov 0.1 \
    --nprocs 30 \
    --output-dir out_filter/

# SameStr Stats

samestr stats \
    --input-files out_filter/*.npz \
    --input-names out_filter/*.names.txt \
    --marker-dir samestr_db/ \
    --nprocs 30 \
    --output-dir out_stats/

# SameStr Compare

samestr compare \
    --input-files out_filter/*.npz \
    --input-names out_filter/*.names.txt \
    --marker-dir samestr_db/ \
    --nprocs 30 \
    --output-dir out_compare/

# SameStr Summarize


samestr summarize \
    --input-dir out_compare/ \
    --tax-profiles-dir out_align/ \
    --marker-dir samestr_db/ \
    --aln-pair-min-overlap 5000 \
    --aln-pair-min-similarity 0.999 \
    --output-dir out_summarize/