#!/bin/bash
#SBATCH -p long
#SBATCH -c 4
#SBATCH -N 1 
#SBATCH --mem-per-cpu 5000
#SBATCH -t 06:00:00 
#SBATCH -o align_launcher.out
#SBATCH -e align_launcher.err

# Paths (Anonymized)
ROOT="/path/to/project_root" # Base path for project data and singularity images
BASEDIR="${ROOT}/project_name" # Main project directory (Anonymized)
TRIMMED="${BASEDIR}/trimmed" # Directory containing quality-trimmed FASTQ files
IMAGES="${ROOT}/images" # Directory containing Singularity container images
# Reference genome path (e.g., GRCh37/hg19)
REF="/path/to/reference/genome/Homo_sapiens.GRCh37.75.dna.primary_assembly.fa"
OUTDIR="${BASEDIR}/alignment" # Output directory for sorted BAM files
LOGDIR="${BASEDIR}/logs/align" # Directory for SLURM job logs

# Create output and log directories if they do not exist
mkdir -p ${OUTDIR}
mkdir -p ${LOGDIR}

cd ${TRIMMED}

# Detect unique samples based on trimmed file names (*_val_1.fq.gz)
for base in $(ls *_val_1.fq.gz | sed 's/\.1_val_1.fq.gz//' | sed 's/\.2_val_1.fq.gz//' | sort -u); do
    R1="${base}.1_val_1.fq.gz" # Forward read file
    R2="${base}.2_val_2.fq.gz" # Reverse read file

    SAMPLE="${base}" # Sample name (e.g., S000021_S14266Nr1)

    echo "Launching BWA alignment for sample: ${SAMPLE}"

    # Submit SLURM job
    sbatch --export=ALL,ROOT=${ROOT},TRIMMED=${TRIMMED},IMAGES=${IMAGES},REF=${REF},OUTDIR=${OUTDIR},SAMPLE=${SAMPLE},R1=${R1},R2=${R2} \
      --output=${LOGDIR}/align_${SAMPLE}.out \
      --error=${LOGDIR}/align_${SAMPLE}.err << 'EOF'
#!/bin/bash
#SBATCH -p fast
#SBATCH -c 4
#SBATCH -N 1 
#SBATCH --mem-per-cpu 5000
#SBATCH -t 06:00:00

# Perform BWA alignment and pipe output directly to Samtools for sorting
singularity exec -B ${ROOT}:${ROOT} ${IMAGES}/bwa_v0.7.17.sif \
  bwa mem -t 4 ${REF} ${TRIMMED}/${R1} ${TRIMMED}/${R2} | \
  singularity exec -B ${ROOT}:${ROOT} ${IMAGES}/samtools_v1.15.sif \
  samtools sort -@ 4 -o ${OUTDIR}/${SAMPLE}.sorted.bam -

# Index the resulting BAM file
singularity exec -B ${ROOT}:${ROOT} ${IMAGES}/samtools_v1.15.sif \
  samtools index ${OUTDIR}/${SAMPLE}.sorted.bam
EOF

done
