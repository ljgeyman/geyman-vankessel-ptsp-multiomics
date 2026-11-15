#!/bin/bash
set -e

# ==============================================================================
# Configuration (Update paths and thread count as needed for your system)
# ==============================================================================
REF_GENOME="refs/refseq_vvul.fna"
METADATA="raw_fastq/anti_smcR_metadata.txt" # optionally anti_rpoB_metadata.txt
READ_DIR="raw_fastq"
THREADS=8

DIR_IDX="01_index"
DIR_BAM="02_alignments"
DIR_MACS="03_macs2_peaks"
DIR_BW="04_bigwigs"
DIR_OUT="05_final_output"

mkdir -p "$DIR_IDX" "$DIR_BAM" "$DIR_MACS" "$DIR_BW" "$DIR_OUT"

# ==============================================================================
# 1. Calculate Effective Genome Size
# ==============================================================================
echo "[INFO] Calculating reference genome size from $REF_GENOME..."
GENOME_SIZE=$(awk '/^[^>]/ {seqlen += length($0)} END {print seqlen}' "$REF_GENOME")
echo "[INFO] Detected Genome Size: $GENOME_SIZE bp"

# ==============================================================================
# 2. Build Bowtie2 Index
# ==============================================================================
echo "[INFO] Building Bowtie2 index..."
if [ ! -f "${DIR_IDX}/ref_index.1.bt2" ]; then
    bowtie2-build "$REF_GENOME" "${DIR_IDX}/ref_index" --threads "$THREADS"
else
    echo "[INFO] Index already exists. Skipping..."
fi

# ==============================================================================
# 3. Per-Sample Alignment, Peak Calling & BigWig Generation
# ==============================================================================
# Expects tab-delimited metadata with header: Strain_ID <tab> IP_Base <tab> Control_Base
# Strain_ID must end in _rep1, _rep2, etc. (supports 2 or 3 replicates per condition)
tail -n +2 "$METADATA" | while IFS=$'\t' read -r STRAIN IP_BASE CTRL_BASE; do
    [ -z "$STRAIN" ] && continue

    echo "------------------------------------------------"
    echo "[INFO] Processing Sample: $STRAIN"
    echo "------------------------------------------------"

    IP_R1="${READ_DIR}/${IP_BASE}_R1.fastq.gz"
    IP_R2="${READ_DIR}/${IP_BASE}_R2.fastq.gz"
    CTRL_R1="${READ_DIR}/${CTRL_BASE}_R1.fastq.gz"
    CTRL_R2="${READ_DIR}/${CTRL_BASE}_R2.fastq.gz"

    for FILE in "$IP_R1" "$IP_R2" "$CTRL_R1" "$CTRL_R2"; do
        if [ ! -f "$FILE" ]; then
            echo "[ERROR] Missing file: $FILE"
            exit 1
        fi
    done

    # Align IP sample (proper pairs only, MAPQ >= 20)
    echo "[INFO] Aligning IP sample for $STRAIN..."
    bowtie2 -p "$THREADS" -x "${DIR_IDX}/ref_index" -1 "$IP_R1" -2 "$IP_R2" \
        | samtools view -@ "$THREADS" -bS -f 2 -q 20 - \
        | samtools sort -@ "$THREADS" -o "${DIR_BAM}/${STRAIN}_IP.bam" -
    samtools index -@ "$THREADS" "${DIR_BAM}/${STRAIN}_IP.bam"

    # Align Control sample (proper pairs only, MAPQ >= 20)
    echo "[INFO] Aligning Control sample for $STRAIN..."
    bowtie2 -p "$THREADS" -x "${DIR_IDX}/ref_index" -1 "$CTRL_R1" -2 "$CTRL_R2" \
        | samtools view -@ "$THREADS" -bS -f 2 -q 20 - \
        | samtools sort -@ "$THREADS" -o "${DIR_BAM}/${STRAIN}_Control.bam" -
    samtools index -@ "$THREADS" "${DIR_BAM}/${STRAIN}_Control.bam"

    # Call paired-end peaks with MACS2
    echo "[INFO] Calling peaks with MACS2 for $STRAIN..."
    macs2 callpeak -t "${DIR_BAM}/${STRAIN}_IP.bam" \
                   -c "${DIR_BAM}/${STRAIN}_Control.bam" \
                   -f BAMPE \
                   -g "$GENOME_SIZE" \
                   -n "$STRAIN" \
                   --outdir "$DIR_MACS" \
                   -q 0.05 \
                   --keep-dup auto

    # Generate CPM-normalized and log2FC BigWigs
    echo "[INFO] Generating normalized BigWig for $STRAIN IP..."
    bamCoverage -b "${DIR_BAM}/${STRAIN}_IP.bam" \
                -o "${DIR_BW}/${STRAIN}_IP.bw" \
                --normalizeUsing CPM \
                --binSize 10 \
                -p "$THREADS"

    echo "[INFO] Generating normalized BigWig for $STRAIN Control..."
    bamCoverage -b "${DIR_BAM}/${STRAIN}_Control.bam" \
                -o "${DIR_BW}/${STRAIN}_Control.bw" \
                --normalizeUsing CPM \
                --binSize 10 \
                -p "$THREADS"

    echo "[INFO] Generating log2FC BigWig (IP vs Control) for $STRAIN..."
    bamCompare -b1 "${DIR_BAM}/${STRAIN}_IP.bam" \
               -b2 "${DIR_BAM}/${STRAIN}_Control.bam" \
               -o "${DIR_BW}/${STRAIN}_IP_lfc.bw" \
               --scaleFactorsMethod readCount \
               --operation log2 \
               --binSize 10 \
               -p "$THREADS"
done

# ==============================================================================
# 4. Replicate Summit Intersection per Condition (Supports 2 or 3 Replicates)
# ==============================================================================
echo "------------------------------------------------"
echo "[INFO] Intersecting replicate peaks for each experimental condition..."

CONDITIONS=$(tail -n +2 "$METADATA" | awk -F'\t' '{print $1}' | sed 's/_rep[0-9]*$//' | sort -u)

for COND in $CONDITIONS; do
    shopt -s nullglob
    SUMMIT_FILES=("${DIR_MACS}/${COND}"_rep*_summits.bed)
    shopt -u nullglob

    NUM_REPS=${#SUMMIT_FILES[@]}

    if [ "$NUM_REPS" -lt 2 ] || [ "$NUM_REPS" -gt 3 ]; then
        echo "[WARN] Condition '$COND' has $NUM_REPS replicate(s) (expected 2 or 3). Skipping intersection."
        continue
    fi

    echo "[INFO] Generating reproducible summits for condition: $COND ($NUM_REPS replicates detected)"
    TMP_SUPPORTED="${DIR_MACS}/${COND}_all_supported_summits.tmp.bed"
    : > "$TMP_SUPPORTED"

    # Retain summits from each replicate that overlap called peaks across all other replicates
    for i in "${!SUMMIT_FILES[@]}"; do
        CURR_SUMMITS="${SUMMIT_FILES[$i]}"
        CURR_FILTERED="${DIR_MACS}/${COND}_rep$((i+1))_supported.tmp.bed"
        cp "$CURR_SUMMITS" "$CURR_FILTERED"

        for j in "${!SUMMIT_FILES[@]}"; do
            [ "$i" -eq "$j" ] && continue
            OTHER_PEAKS="${SUMMIT_FILES[$j]/_summits.bed/_peaks.narrowPeak}"
            if [ ! -f "$OTHER_PEAKS" ]; then
                echo "[ERROR] Missing peak file: $OTHER_PEAKS"
                exit 1
            fi
            bedtools intersect -u -a "$CURR_FILTERED" -b "$OTHER_PEAKS" > "${CURR_FILTERED}.next"
            mv "${CURR_FILTERED}.next" "$CURR_FILTERED"
        done

        cat "$CURR_FILTERED" >> "$TMP_SUPPORTED"
        rm -f "$CURR_FILTERED"
    done

    sort -k1,1 -k2,2n "$TMP_SUPPORTED" | bedtools merge > "${DIR_MACS}/${COND}_reproducible_summits.bed"
    rm -f "$TMP_SUPPORTED"
done

# ==============================================================================
# 5. Master Consensus Peak Sets & Aggregate Signal Profiling
# ==============================================================================
echo "------------------------------------------------"
echo "[INFO] Creating master consensus peak sets for R and deepTools..."

CONSENSUS_SUMMITS="${DIR_MACS}/consensus_summits_merged.bed"
cat "${DIR_MACS}"/*_reproducible_summits.bed | sort -k1,1 -k2,2n | bedtools merge > "$CONSENSUS_SUMMITS"

CONSENSUS_PEAKS="${DIR_MACS}/consensus_peaks_merged.bed"
cat "${DIR_MACS}"/*_peaks.narrowPeak | sort -k1,1 -k2,2n | bedtools merge > "$CONSENSUS_PEAKS"

# Consolidate all narrowPeak calls into a single tab-delimited table for R
FINAL_TABLE="${DIR_OUT}/combined_macs2_peaks.txt"
echo -e "chrom\tchromStart\tchromEnd\tname\tscore\tstrand\tfold_change\tlog10_pval\tlog10_qval\tpeak_summit\tStrain" > "$FINAL_TABLE"

for PEAK_FILE in "${DIR_MACS}"/*.narrowPeak; do
    if [ -f "$PEAK_FILE" ]; then
        STRAIN_NAME=$(basename "$PEAK_FILE" _peaks.narrowPeak)
        awk -v strain="$STRAIN_NAME" 'BEGIN{OFS="\t"} {print $0, strain}' "$PEAK_FILE" >> "$FINAL_TABLE"
    fi
done

# Compute coverage matrix and profile plot centered on reproducible summits
computeMatrix reference-point \
    --referencePoint center \
    -b 1000 -a 1000 \
    -R "$CONSENSUS_SUMMITS" \
    -S "${DIR_BW}"/*_IP.bw \
    -o "${DIR_OUT}/matrix_center_1kb.gz" \
    --skipZeros \
    -p "$THREADS"

plotProfile \
    -m "${DIR_OUT}/matrix_center_1kb.gz" \
    -out "${DIR_OUT}/aggregate_peak_profile.pdf" \
    --plotType se \
    --perGroup \
    --plotTitle "Aggregate ChIP-seq Signal" \
    --regionsLabel "Peaks" \
    --yAxisLabel "Signal (CPM)" \
    --refPointLabel "center"

echo "[SUCCESS] Pipeline complete. Replicate-filtered consensus files and plots are in: $DIR_MACS and $DIR_OUT"
