#!/bin/bash

# This script documents the commands used to run the samestr_flow Nextflow pipeline,
# used for benchmarking SamestrGal against Nextflow on the same input datasets.

# ============================================================
# 1. Download raw sequencing data
# ============================================================

# All data can be accessed from the ENA project in accession PRJEB39023, which contains a total of 
# 8 patient samples that underwent FMT treatment. For the benchmarking only 5 cases were used:

# Case 16 (Pre-FMT: 16A, Donor: 16B, Post-FMT: 16D)

wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/run/ERR429/ERR4290868/16A.R1.fastq.gz
wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/run/ERR429/ERR4290868/16A.R2.fastq.gz
wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/run/ERR429/ERR4290678/16B.R1.fastq.gz
wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/run/ERR429/ERR4290678/16B.R2.fastq.gz
wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/run/ERR429/ERR4290679/16D.R1.fastq.gz
wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/run/ERR429/ERR4290679/16D.R2.fastq.gz

# Case 19 (Pre-FMT: 19A, Donor: 19B, Post-FMT: 19C, 19G)

wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/run/ERR429/ERR4290682/19A.R1.fastq.gz
wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/run/ERR429/ERR4290682/19A.R2.fastq.gz
wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/run/ERR429/ERR4290680/19B.R2.fastq.gz
wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/run/ERR429/ERR4290680/19B.R1.fastq.gz
wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/run/ERR429/ERR4290684/19C.R1.fastq.gz
wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/run/ERR429/ERR4290684/19C.R2.fastq.gz
wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/run/ERR429/ERR4291114/19G.R1.fastq.gz
wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/run/ERR429/ERR4291114/19G.R2.fastq.gz

# Case 20 (Pre-FMT: 20A, Donor: 20B, Post-FMT: 20E)

wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/run/ERR429/ERR4290683/20A.R2.fastq.gz
wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/run/ERR429/ERR4290683/20A.R1.fastq.gz
wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/run/ERR429/ERR4290681/20B.R2.fastq.gz
wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/run/ERR429/ERR4290681/20B.R1.fastq.gz
wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/run/ERR429/ERR4290869/20E.R1.fastq.gz
wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/run/ERR429/ERR4290869/20E.R2.fastq.gz

# Case 28 (Pre-FMT: 28A, Donor: 28B, Post-FMT: 28C)

wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/run/ERR429/ERR4291165/28A.R2.fastq.gz
wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/run/ERR429/ERR4291165/28A.R1.fastq.gz
wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/run/ERR429/ERR4291164/28B.R1.fastq.gz
wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/run/ERR429/ERR4291164/28B.R2.fastq.gz
wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/run/ERR429/ERR4291154/28C.R1.fastq.gz
wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/run/ERR429/ERR4291154/28C.R2.fastq.gz

# Case 54 (Pre-FMT: 54A, Donor: 54B, Post-FMT: 54C)

wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/run/ERR429/ERR4291167/54A.R1.fastq.gz
wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/run/ERR429/ERR4291167/54A.R2.fastq.gz
wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/run/ERR429/ERR4291171/54B.R1.fastq.gz
wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/run/ERR429/ERR4291171/54B.R2.fastq.gz
wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/run/ERR429/ERR4291954/54C.R1.fastq.gz
wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/run/ERR429/ERR4291954/54C.R2.fastq.gz

# ============================================================
# 2. Create the input_data directory
# ============================================================

# Nextflow expects one subfolder per sample, with reads renamed to exactly R1.fastq.gz / R2.fastq.gz. Example:

# input_data/
# ├── 16A/
# │   ├── R1.fastq.gz
# │   └── R2.fastq.gz
# ├── 16B/
# │   ├── R1.fastq.gz
# │   └── R2.fastq.gz
# └── ...

# ============================================================
# 3. Download reference databases
# ============================================================

# Kraken2 human database for host removal

https://zenodo.org/records/8339700

# MetaPhlAn database (mpa_vJan25_CHOCOPhlAnSGB_202503)

metaphlan --install --index mpa_vJan25_CHOCOPhlAnSGB_202503 --db_dir db_MetaPhlAn/

# Generate the SameStr reference database from the MetaPhlAn markers. samestr_flow expects this folder to be 
# named SameStr_mpa_vJan25_CHOCOPhlAnSGB_202503

samestr db \
    --markers-info db_MetaPhlAn/mpa_vJan25_CHOCOPhlAnSGB_202503.pkl \
    --markers-fasta db_MetaPhlAn/mpa_vJan25_CHOCOPhlAnSGB_202503.fna.bz2 \
    --db-version db_MetaPhlAn/mpa_latest \
    --output-dir SameStr_mpa_vJan25_CHOCOPhlAnSGB_202503/


# ============================================================
# 4. Modify SameStr parameters
# ============================================================

# samestr_flow hardcodes the parameters of the individual SameStr tools in:
# https://github.com/grp-bork/samestr_flow/blob/main/metaphlow/modules/profilers/samestr.nf
# This file must be replaced directly with the corrected version in scripts/nextflow_files/
# to be able to run nextflow with the same parameters as in SamestrGal.

# ============================================================
# 5. Configure params.yml
# ============================================================

# Use the params.yml provided in scripts/nextflow_files/ and replace the paths.

# ============================================================
# 6. Run Nextflow workflow
# ============================================================

nextflow run main.nf \
    -params-file params.yml \
    -c nextflow.config \
    -resume \
    --remove_host kraken2

# The SameStr outputs will be in results/. Files like pre-processing files and SameStr Convert outputs are 
# not copied there. To see those, look inside Nextflow's work/ directory. 