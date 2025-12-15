#!/bin/bash
#SBATCH -p long
#SBATCH -c 1
#SBATCH -N 1
#SBATCH --mem-per-cpu 4000
#SBATCH -t 02:00:00
#SBATCH -o cnv_launcher.out
#SBATCH -e cnv_launcher.err

# Paths (Anonymized)
ROOT="/path/to/project_root" # Base path for project data and singularity images
BASEDIR="${ROOT}/project_name" # Main project directory (Anonymized)
DEDUP="${BASEDIR}/deduplicated_rg" # Input directory containing deduplicated and RG-tagged BAMs
IMAGES="${ROOT}/images" # Directory containing Singularity container images
REF="${ROOT}/path/to/reference/Homo_sapiens.GRCh37.75.dna.primary_assembly.fa" # Reference genome FASTA
BED="${BASEDIR}/resources/S000021_hg19_targets.bed" # Target regions BED file
OUTDIR="${BASEDIR}/cnv_gatk" # Output directory for GATK CNV results
LOGDIR="${BASEDIR}/logs/gatk_cnv" # Directory for SLURM job logs

# Create output and log directories if they do not exist
mkdir -p ${OUTDIR}
mkdir -p ${LOGDIR}

# --- Step 1: Pre-processing Intervals (Run Once) ---
# This step defines and annotates the genomic intervals based on the target BED file.
if [[ ! -f ${OUTDIR}/preprocessed_intervals.interval_list ]]; then
    echo "Creating and Annotating Intervals (Run Once)..."
    
    # 1a. PreprocessIntervals
    singularity exec -B ${ROOT}:${ROOT} ${IMAGES}/gatk4_v4.3.0.sif \
        gatk PreprocessIntervals \
        -R ${REF} \
        -L ${BED} \
        --bin-length 0 \
        --interval-merging-rule OVERLAPPING_ONLY \
        -O ${OUTDIR}/preprocessed_intervals.interval_list

    # 1b. AnnotateIntervals (for potential GC bias correction, etc.)
    singularity exec -B ${ROOT}:${ROOT} ${IMAGES}/gatk4_v4.3.0.sif \
        gatk AnnotateIntervals \
        -R ${REF} \
        -L ${OUTDIR}/preprocessed_intervals.interval_list \
        -O ${OUTDIR}/annotated_intervals.tsv
fi

# --- Step 2: Launch CNV Analysis per Sample ---
cd ${DEDUP}
for BAM in *.bam; do
    SAMPLE=$(basename ${BAM} .bam)

    echo "Submitting CNV job for sample: ${SAMPLE}"

    # Submit SLURM job for per-sample CNV steps
    sbatch --export=ALL,ROOT=${ROOT},DEDUP=${DEDUP},IMAGES=${IMAGES},OUTDIR=${OUTDIR},SAMPLE=${SAMPLE},REF=${REF} \
        --output=${LOGDIR}/gatkcnv_${SAMPLE}.out \
        --error=${LOGDIR}/gatkcnv_${SAMPLE}.err << 'EOF'
#!/bin/bash
#SBATCH -p fast
#SBATCH -c 2
#SBATCH -N 1
#SBATCH --mem-per-cpu 5000
#SBATCH -t 04:00:00

# Step 2a: CollectReadCounts (Quantification of reads per interval)
singularity exec -B ${ROOT}:${ROOT} ${IMAGES}/gatk4_v4.3.0.sif gatk CollectReadCounts \
    -I ${DEDUP}/${SAMPLE}.bam \
    -L ${OUTDIR}/preprocessed_intervals.interval_list \
    --interval-merging-rule OVERLAPPING_ONLY \
    -O ${OUTDIR}/${SAMPLE}.counts.hdf5

# Step 2b: DenoiseReadCounts (Normalization and noise reduction)
singularity exec -B ${ROOT}:${ROOT} ${IMAGES}/gatk4_v4.3.0.sif gatk DenoiseReadCounts \
    -I ${OUTDIR}/${SAMPLE}.counts.hdf5 \
    --standardized-copy-ratios ${OUTDIR}/${SAMPLE}.standardizedCR.tsv \
    --denoised-copy-ratios ${OUTDIR}/${SAMPLE}.denoisedCR.tsv

# Step 2c: ModelSegments (Segmentation and calling of copy number events)
singularity exec -B ${ROOT}:${ROOT} ${IMAGES}/gatk4_v4.3.0.sif gatk ModelSegments \
    --denoised-copy-ratios ${OUTDIR}/${SAMPLE}.denoisedCR.tsv \
    --output ${OUTDIR}/${SAMPLE}_segments \
    --output-prefix ${SAMPLE}

EOF

done
