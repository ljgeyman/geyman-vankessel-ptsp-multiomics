#!/bin/bash
set -e

# ==============================================================================
# Configuration (Update paths and thread count as needed for your system)
# ==============================================================================
METADATA="metadata.txt"
THREADS=14
RAW_DIR="raw_fastq"
OUT_DIR="analysis_results"

R1_SUFFIX="_R1.fastq.gz"
R2_SUFFIX="_R2.fastq.gz"

export _JAVA_OPTIONS="-Xmx16G"

# Reference paths (relative to repository root; see README for RefSeq accessions)
ADAPTERS="${CONDA_PREFIX}/share/trimmomatic/adapters/TruSeq3-PE.fa"
GENOME_FASTA="refs/refseq_vvul.fna"
GTF_FILE="refs/refseq_vvul.gtf"
BOWTIE_IDX="refs/bowtie_indexes/vvul_refseq"
# ==============================================================================

# --- Preflight Checks ---
if [[ ! -f "$METADATA" ]]; then
    echo "ERROR: Cannot find $METADATA in the current directory. Exiting."
    exit 1
fi

if ! command -v bedtools &> /dev/null; then
    echo "ERROR: bedtools is missing. Please activate the RNA-seq conda environment."
    exit 1
fi

if ! command -v bamCoverage &> /dev/null; then
    echo "ERROR: bamCoverage (deepTools) is missing. Please activate the RNA-seq conda environment."
    exit 1
fi

echo "------------------------------------------------"
echo "Starting pipeline using: $METADATA"
echo "------------------------------------------------"

# --- 1. Build Bowtie2 Index (if not present) ---
echo "Checking for Bowtie2 index at: ${BOWTIE_IDX}"
if [[ ! -f "${BOWTIE_IDX}.1.bt2" ]]; then
    echo "Index not found. Building Bowtie2 index from $GENOME_FASTA..."
    mkdir -p "$(dirname "$BOWTIE_IDX")"
    bowtie2-build --threads "$THREADS" "$GENOME_FASTA" "$BOWTIE_IDX"
    echo "Bowtie2 index built successfully!"
else
    echo "Bowtie2 index found. Moving on."
fi

mkdir -p "${OUT_DIR}"/{fastqc_raw,trimmed,fastqc_trimmed,bam,counts,multiqc,bigwig}
echo "------------------------------------------------"

# --- 2. Per-Sample QC, Trimming, Alignment & Coverage ---
while read -r SAMPLE; do
    [[ -z "$SAMPLE" ]] && continue
    SAMPLE=$(echo "$SAMPLE" | tr -d '\r')
    
    echo ">>> Processing sample: $SAMPLE"
    
    READ1="${RAW_DIR}/${SAMPLE}${R1_SUFFIX}"
    READ2="${RAW_DIR}/${SAMPLE}${R2_SUFFIX}"

    if [[ ! -f "$READ1" ]] || [[ ! -f "$READ2" ]]; then
        echo "ERROR: Cannot find $READ1 or $READ2. Skipping $SAMPLE."
        continue
    fi

    BAM_OUT="${OUT_DIR}/bam/${SAMPLE}_sorted.bam"
    
    if [[ ! -f "$BAM_OUT" ]]; then
        # Raw FastQC
        fastqc -t "$THREADS" -o "${OUT_DIR}/fastqc_raw" -q "$READ1" "$READ2"

        # Conditional Trimmomatic (applies ILLUMINACLIP only if FastQC flags adapter content)
        TRIM_R1="${OUT_DIR}/trimmed/${SAMPLE}_paired_R1.fq.gz"
        TRIM_R2="${OUT_DIR}/trimmed/${SAMPLE}_paired_R2.fq.gz"
        UNP_R1="${OUT_DIR}/trimmed/${SAMPLE}_unpaired_R1.fq.gz"
        UNP_R2="${OUT_DIR}/trimmed/${SAMPLE}_unpaired_R2.fq.gz"

        R1_BASE=$(basename "$READ1" .fastq.gz)
        FASTQC_ZIP="${OUT_DIR}/fastqc_raw/${R1_BASE}_fastqc.zip"
        ADAPTER_STATUS=$(unzip -p "$FASTQC_ZIP" "*/fastqc_data.txt" | grep ">>Adapter Content" | cut -f2)

        if [[ "$ADAPTER_STATUS" == "fail" || "$ADAPTER_STATUS" == "warn" ]]; then
            echo ">>> Adapters detected in $SAMPLE ($ADAPTER_STATUS). Running Trimmomatic WITH ILLUMINACLIP."
            trimmomatic PE -threads "$THREADS" -quiet \
                "$READ1" "$READ2" \
                "$TRIM_R1" "$UNP_R1" "$TRIM_R2" "$UNP_R2" \
                ILLUMINACLIP:"${ADAPTERS}":2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:18
        else
            echo ">>> No significant adapters detected in $SAMPLE ($ADAPTER_STATUS). Running quality trim only."
            trimmomatic PE -threads "$THREADS" -quiet \
                "$READ1" "$READ2" \
                "$TRIM_R1" "$UNP_R1" "$TRIM_R2" "$UNP_R2" \
                LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:18
        fi

        # Post-trimming FastQC
        fastqc -t "$THREADS" -o "${OUT_DIR}/fastqc_trimmed" -q "$TRIM_R1" "$TRIM_R2"

        # Bowtie2 alignment & coordinate sorting
        echo "Aligning reads to genome..."
        bowtie2 -p "$THREADS" -x "$BOWTIE_IDX" \
            --sensitive --end-to-end \
            -X 1000 --dovetail \
            --no-unal --no-mixed \
            -1 "$TRIM_R1" -2 "$TRIM_R2" \
            2> "${OUT_DIR}/bam/${SAMPLE}_bowtie2.log" | \
            samtools view -@ "$THREADS" -bS - | \
            samtools sort -@ "$THREADS" -o "$BAM_OUT"
            
        samtools index -@ "$THREADS" "$BAM_OUT"
    else
        echo ">>> BAM file already exists for $SAMPLE. Skipping alignment steps..."
    fi

    # Stranded CPM-normalized BigWigs
    BW_FWD="${OUT_DIR}/bigwig/${SAMPLE}_fwd.bw"
    BW_REV="${OUT_DIR}/bigwig/${SAMPLE}_rev.bw"
    
    if [[ ! -f "$BW_FWD" ]] || [[ ! -f "$BW_REV" ]]; then
        echo ">>> Generating stranded CPM-normalized BigWigs for $SAMPLE..."
        bamCoverage -b "$BAM_OUT" -o "$BW_FWD" -p "$THREADS" --normalizeUsing CPM --filterRNAstrand forward
        bamCoverage -b "$BAM_OUT" -o "$BW_REV" -p "$THREADS" --normalizeUsing CPM --filterRNAstrand reverse
    else
        echo ">>> Stranded BigWigs already exist for $SAMPLE. Skipping..."
    fi

    echo ">>> Finished $SAMPLE"
    echo "------------------------------------------------"
    
done < "$METADATA"

# --- 3. Joint FeatureCounts (Reverse-stranded, -s 2) ---
echo ">>> Running FeatureCounts on all samples..."
featureCounts -T "$THREADS" -p -s 2 \
    -t gene \
    -g locus_tag \
    -a "$GTF_FILE" \
    -o "${OUT_DIR}/counts/all_samples_counts.txt" \
    "${OUT_DIR}"/bam/*_sorted.bam

# --- 4. Compile MultiQC Report ---
echo ">>> Compiling final MultiQC report..."
multiqc "${OUT_DIR}" -o "${OUT_DIR}/multiqc" -n final_multiqc_report.html -q

echo "Pipeline fully complete! Check the ${OUT_DIR} directory for results."
