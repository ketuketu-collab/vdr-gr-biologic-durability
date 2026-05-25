#!/usr/bin/env bash
# Run this once hg38 index download is complete.
# 1. Unzip hg38 bowtie2 index
# 2. Align VDR ChIP-seq FASTQ
# 3. Generate bigWig + MACS2 peaks
# 4. Create pile-up figure
set -euo pipefail

CONDA_BASE=~/miniforge3
source $CONDA_BASE/etc/profile.d/conda.sh
conda activate chipseq

GENOMES=~/genomes
PROJ=/Volumes/M4_SSD/projects/tlr_chipseq
LOGS=$PROJ/logs

echo "=== Phase 1: Unzip hg38 bowtie2 index ==="
if [ ! -f $GENOMES/GRCh38_noalt_as/GRCh38_noalt_as.1.bt2 ]; then
    if [ ! -f $GENOMES/GRCh38_noalt_as.zip ]; then
        echo "ERROR: $GENOMES/GRCh38_noalt_as.zip not found"
        echo "Download with: curl -L https://genome-idx.s3.amazonaws.com/bt/GRCh38_noalt_as.zip -o ~/genomes/GRCh38_noalt_as.zip"
        exit 1
    fi
    echo "  Unzipping ..."
    cd $GENOMES && unzip GRCh38_noalt_as.zip
    echo "  Done. Index at: $GENOMES/GRCh38_noalt_as/"
else
    echo "  hg38 index already extracted"
fi

echo "=== Phase 2: VDR ChIP-seq alignment + bigWig ==="
bash $PROJ/scripts/06_align_vdr.sh 2>&1 | tee $LOGS/vdr_alignment.log

echo "=== Phase 3: BigWig pile-up figure ==="
python $PROJ/scripts/07_bigwig_pileup.py 2>&1 | tee $LOGS/vdr_pileup.log

echo ""
echo "=== ALL DONE: $(date) ==="
ls -lh $PROJ/results/coverage/*.bw $PROJ/results/peaks/*.narrowPeak 2>/dev/null
ls -lh $PROJ/results/figures/*.pdf 2>/dev/null
