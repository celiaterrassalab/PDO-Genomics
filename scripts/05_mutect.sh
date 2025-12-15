#!/bin/bash
#SBATCH -p long
#SBATCH -c 4
#SBATCH -N 1 
#SBATCH --mem-per-cpu 6000
#SBATCH -t 06:00:00 
#SBATCH -o mutect2_launcher.out
#SBATCH -e mutect2_launcher.err

# Paths (Anonymized)
ROOT="/path/to/project_root" # Base path for project data and singularity images
BASEDIR="${ROOT}/project_name" # Main project directory (Anonymized)
DEDUP="${BASEDIR}/deduplicated_rg" # Input directory containing deduplicated and Read Group-tagged BAMs
IMAGES="${ROOT}/images" # Directory containing Singularity container images
REF="${ROOT}/path/to/reference/Homo_sapiens.GRCh37.75.dna.primary_assembly.fa" # Reference genome FASTA
BED="${BASEDIR}/resources/S000021_hg19_targets.bed" # Target regions BED file for focused calling
OUTDIR="${BASEDIR}/variants" # Output directory for VCF files
LOGDIR="${BASEDIR}/logs/mutect2" # Directory for SLURM job logs

# Create output and log directories if they do not exist
mkdir -p ${OUTDIR}
mkdir -p ${LOGDIR}

cd ${DEDUP}

# Iterate over BAM files for variant calling
for BAM in *.bam; do
    FILE=${DEDUP}/${BAM}
    
    # Extract the sample name (SM) from the @RG header line in the BAM file
    SAMPLE_RG=$(singularity exec -B ${ROOT}:${ROOT} ${IMAGES}/samtools_v1.15.sif \
        samtools view -H ${FILE} | grep '^@RG' | sed -n 's/.*SM:\([^\t]*\).*/\1/p' | head -n 1)

    # Sanity check: verify that the SM tag was found
    if [[ -z "$SAMPLE_RG" ]]; then
        echo "❌ SM: tag not found in ${BAM}, skipping sample."
        continue
    fi

    echo "✅ Launching variant calling for: ${SAMPLE_RG} (File: ${BAM})"

    # Submit SLURM job for Mutect2 and FilterMutectCalls
    sbatch --export=ALL,ROOT=${ROOT},DEDUP=${DEDUP},IMAGES=${IMAGES},REF=${REF},BED=${BED},OUTDIR=${OUTDIR},BAM=${BAM},SAMPLE_RG=${SAMPLE_RG} \
      --output=${LOGDIR}/mutect2_${SAMPLE_RG}.out \
      --error=${LOGDIR}/mutect2_${SAMPLE_RG}.err << 'EOF'
#!/bin/bash
#SBATCH -p fast
#SBATCH -c 4
#SBATCH -N 1 
#SBATCH --mem-per-cpu 6000
#SBATCH -t 06:00:00

# Step 1: Run Mutect2 in tumor-only mode
singularity exec -B ${ROOT}:${ROOT} ${IMAGES}/gatk4_v4.3.0.sif \
gatk Mutect2 \
-R ${REF} \
-I ${DEDUP}/${BAM} \
-tumor ${SAMPLE_RG} \
-L ${BED} \
--native-pair-hmm-threads 4 \
-O ${OUTDIR}/${SAMPLE_RG}.raw.vcf.gz

# Step 2: Filter the raw variant calls (critical GATK step for quality)
singularity exec -B ${ROOT}:${ROOT} ${IMAGES}/gatk4_v4.3.0.sif \
gatk FilterMutectCalls \
-V ${OUTDIR}/${SAMPLE_RG}.raw.vcf.gz \
-R ${REF} \
-O ${OUTDIR}/${SAMPLE_RG}.filtered.vcf.gz
EOF

done
