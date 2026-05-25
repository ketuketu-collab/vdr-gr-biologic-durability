#!/usr/bin/env bash
# VDR ChIP-seq alignment pipeline: GSE89431 THP-1 cells (hg38)
# Run after hg38 bowtie2 index is extracted.
set -euo pipefail

CONDA_BASE=~/miniforge3
CONDA_ENV=chipseq
source $CONDA_BASE/etc/profile.d/conda.sh
conda activate $CONDA_ENV

PROJ=/Volumes/M4_SSD/projects/tlr_chipseq
DATA=$PROJ/data/vdr
RESULTS=$PROJ/results
LOGS=$PROJ/logs
INDEX=~/genomes/GRCh38_noalt_as/GRCh38_noalt_as
THREADS=8

mkdir -p $RESULTS/{coverage,peaks} $LOGS

step_done() { touch $LOGS/step_${1}_done.txt; }
step_skip() { [[ -f $LOGS/step_${1}_done.txt ]]; }

echo "=== VDR ChIP-seq alignment pipeline: $(date) ==="

# ── Samples ──────────────────────────────────────────────────────────────────
# SRR4828890 = Vehicle (EtOH) 2h rep1 → control
# SRR4828891 = 1,25(OH)2D3 2h rep1    → VitD treatment
# SRR4828892 = Vehicle 2h rep2
# SRR4828893 = VitD 2h rep2
declare -A SAMPLES=(
    [SRR4828890]=vehicle_rep1
    [SRR4828891]=vitd_rep1
    [SRR4828892]=vehicle_rep2
    [SRR4828893]=vitd_rep2
)

# ── Step A: Alignment ─────────────────────────────────────────────────────────
for SRR in "${!SAMPLES[@]}"; do
    NAME=${SAMPLES[$SRR]}
    FASTQ=$DATA/${SRR}.fastq
    BAM=$RESULTS/${NAME}.bam

    if [[ ! -f $FASTQ ]]; then
        echo "WARNING: $FASTQ not found, skipping $NAME"
        continue
    fi
    if [[ -f $BAM ]]; then
        echo "  $NAME BAM exists, skipping alignment"
        continue
    fi
    echo "  Aligning $NAME ($SRR) ..."
    bowtie2 -x $INDEX -U $FASTQ -p $THREADS --no-unal \
        2>$LOGS/${NAME}_bowtie2.log \
        | samtools sort -@ $THREADS -o $BAM
    samtools index $BAM
    # Log stats
    TOTAL=$(samtools flagstat $BAM | grep "in total" | cut -d' ' -f1)
    MAP_PCT=$(grep "overall alignment rate" $LOGS/${NAME}_bowtie2.log | grep -oP '\d+\.\d+(?=%)')
    echo "    $NAME: $TOTAL reads, ${MAP_PCT}% aligned"
done

step_done "align"
echo "Alignment done: $(date)"

# ── Step B: bigWig (CPM normalize) ───────────────────────────────────────────
for SRR in "${!SAMPLES[@]}"; do
    NAME=${SAMPLES[$SRR]}
    BAM=$RESULTS/${NAME}.bam
    BW=$RESULTS/coverage/${NAME}.bw

    [[ ! -f $BAM ]] && continue
    [[ -f $BW ]] && echo "  $NAME.bw exists, skipping" && continue

    echo "  bamCoverage $NAME ..."
    bamCoverage -b $BAM -o $BW \
        --normalizeUsing CPM --binSize 10 -p $THREADS \
        --ignoreDuplicates --minMappingQuality 20
done
step_done "bigwig"
echo "bigWig done: $(date)"

# ── Step C: MACS2 peak calling (VitD vs Vehicle) ─────────────────────────────
# Merge reps for peak calling
if [[ ! -f $RESULTS/coverage/vitd_merged.bam ]]; then
    samtools merge $RESULTS/coverage/vitd_merged.bam \
        $RESULTS/vitd_rep1.bam $RESULTS/vitd_rep2.bam 2>/dev/null || \
        cp $RESULTS/vitd_rep1.bam $RESULTS/coverage/vitd_merged.bam
    samtools index $RESULTS/coverage/vitd_merged.bam
    samtools merge $RESULTS/coverage/vehicle_merged.bam \
        $RESULTS/vehicle_rep1.bam $RESULTS/vehicle_rep2.bam 2>/dev/null || \
        cp $RESULTS/vehicle_rep1.bam $RESULTS/coverage/vehicle_merged.bam
    samtools index $RESULTS/coverage/vehicle_merged.bam
fi

if [[ ! -f $RESULTS/peaks/vdr_thp1_2h_peaks.narrowPeak ]]; then
    echo "  MACS2 peak calling ..."
    macs2 callpeak \
        -t $RESULTS/coverage/vitd_merged.bam \
        -c $RESULTS/coverage/vehicle_merged.bam \
        -f BAM -g hs --nomodel --extsize 200 \
        -n vdr_thp1_2h --outdir $RESULTS/peaks/ \
        2>$LOGS/macs2_vdr.log
fi
step_done "peaks"
echo "MACS2 done: $(date)"

# ── Step D: Merged bigWig for visualization ────────────────────────────────
for cond in vitd vehicle; do
    BW_MERGED=$RESULTS/coverage/${cond}_merged.bw
    [[ -f $BW_MERGED ]] && continue
    BAM=$RESULTS/coverage/${cond}_merged.bam
    [[ ! -f $BAM ]] && continue
    bamCoverage -b $BAM -o $BW_MERGED \
        --normalizeUsing CPM --binSize 10 -p $THREADS \
        --ignoreDuplicates --minMappingQuality 20
done
step_done "merged_bigwig"
echo "=== Pipeline complete: $(date) ==="
echo "Output bigWigs:"
ls -lh $RESULTS/coverage/*.bw 2>/dev/null
echo "Output peaks:"
ls -lh $RESULTS/peaks/*.narrowPeak 2>/dev/null
