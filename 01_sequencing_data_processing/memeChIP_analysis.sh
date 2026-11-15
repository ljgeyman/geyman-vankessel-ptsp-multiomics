#!/usr/bin/env bash
set -euo pipefail

# ==============================================================================
# ChIP-Seq Replicate Summit Intersection & MEME-ChIP Motif Discovery (2 or 3 Reps)
# ==============================================================================

# --- Configuration ---
GENOME_FASTA="refs/refseq_vvul.fna"
COLLECTF_DB="refs/motif_databases/PROKARYOTE/collectf.meme"

PEAK_DIR="03_macs2_peaks"
CONDITION="wt_PTSP"   # Automatically detects ${CONDITION}_rep*_summits.bed (2 or 3 reps)

OUT_DIR="PTSP_meme_chip_analysis"
TOP_N=1000
WINDOW_HALF_SIZE=170  # Expands 1 bp summit +/- 170 bp (341 bp total window)
THREADS=10

mkdir -p "${OUT_DIR}"

# --- 1. Generate Chromosome Sizes ---
echo "[1/5] Generating chromosome sizes file..."
samtools faidx "${GENOME_FASTA}"
cut -f1,2 "${GENOME_FASTA}.fai" > "${OUT_DIR}/genome.chrom.sizes"

# --- 2. Identify Reproducible Summits Across 2 or 3 Replicates (+/- 100 bp tolerance) ---
shopt -s nullglob
SUMMIT_FILES=("${PEAK_DIR}/${CONDITION}"_rep*_summits.bed)
shopt -u nullglob

NUM_REPS=${#SUMMIT_FILES[@]}
if [ "$NUM_REPS" -lt 2 ] || [ "$NUM_REPS" -gt 3 ]; then
    echo "[ERROR] Expected 2 or 3 replicate summit files for '${CONDITION}', found ${NUM_REPS}."
    exit 1
fi

echo "[2/5] Identifying consensus summits across ${NUM_REPS} replicates for ${CONDITION}..."
cp "${SUMMIT_FILES[0]}" "${OUT_DIR}/reproducible_summits.bed"

for ((i=1; i<NUM_REPS; i++)); do
    REP_IDX=$((i+1))
    bedtools slop \
        -i "${SUMMIT_FILES[$i]}" \
        -g "${OUT_DIR}/genome.chrom.sizes" \
        -b 100 > "${OUT_DIR}/rep${REP_IDX}_expanded.bed"

    bedtools intersect \
        -u \
        -a "${OUT_DIR}/reproducible_summits.bed" \
        -b "${OUT_DIR}/rep${REP_IDX}_expanded.bed" > "${OUT_DIR}/reproducible_summits.tmp.bed"

    mv "${OUT_DIR}/reproducible_summits.tmp.bed" "${OUT_DIR}/reproducible_summits.bed"
done

# --- 3. Select Top N Summits by MACS2 Score & Expand Windows ---
echo "[3/5] Selecting top ${TOP_N} summits and expanding windows..."
sort -k5,5nr "${OUT_DIR}/reproducible_summits.bed" | \
    head -n "${TOP_N}" > "${OUT_DIR}/top_${TOP_N}_summits.bed"

bedtools slop \
    -i "${OUT_DIR}/top_${TOP_N}_summits.bed" \
    -g "${OUT_DIR}/genome.chrom.sizes" \
    -b "${WINDOW_HALF_SIZE}" > "${OUT_DIR}/top_${TOP_N}_windows.bed"

# --- 4. Extract Genomic FASTA Sequences ---
echo "[4/5] Extracting FASTA sequences..."
bedtools getfasta \
    -fi "${GENOME_FASTA}" \
    -bed "${OUT_DIR}/top_${TOP_N}_windows.bed" \
    -fo "${OUT_DIR}/top_${TOP_N}_peaks.fasta"

# --- 5. Build 1st-Order Markov Background & Run MEME-ChIP ---
echo "[5/5] Building background model and running MEME-ChIP..."
fasta-get-markov -m 1 "${GENOME_FASTA}" "${OUT_DIR}/bg_model.txt"

meme-chip \
    -oc "${OUT_DIR}/meme_chip_results" \
    -dna \
    -bfile "${OUT_DIR}/bg_model.txt" \
    -db "${COLLECTF_DB}" \
    -minw 21 \
    -maxw 21 \
    -meme-mod zoops \
    -meme-nmotifs 3 \
    -meme-p "${THREADS}" \
    "${OUT_DIR}/top_${TOP_N}_peaks.fasta"

echo "Workflow complete. Results available in: ${OUT_DIR}/meme_chip_results/index.html"
