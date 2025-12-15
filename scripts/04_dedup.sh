#!/bin/bash
#SBATCH -p long
#SBATCH -c 2
#SBATCH -N 1 
#SBATCH --mem-per-cpu 6000
#SBATCH -t 03:00:00 
#SBATCH -o dedup_launcher.out
#SBATCH -e dedup_launcher.err

# Paths (Anonymized)
ROOT="/path/to/project_root" # Base path for project data and singularity images
BASEDIR="${ROOT}/project_name" # Main project directory (Anonymized)
ALIGN="${BASEDIR}/alignment" # Input directory containing sorted BAM files
IMAGES="${ROOT}/images" # Directory containing Singularity container images
OUTDIR="${BASEDIR}/deduplicated" # Output directory for deduplicated BAM files
LOGDIR="${BASEDIR}/logs/dedup" # Directory for SLURM job logs

# Create output and log directories if they do not exist
mkdir -p ${OUTDIR}
mkdir -p ${LOGDIR}

cd ${ALIGN}

# Iterate over aligned BAM files
for BAM in *.sorted.bam; do
SAMPLE=$(basename ${BAM} .sorted.bam)

echo "Launching Picard MarkDuplicates for sample: ${SAMPLE}"

# Submit SLURM job
sbatch --export=ALL,ROOT=${ROOT},ALIGN=${ALIGN},IMAGES=${IMAGES},OUTDIR=${OUTDIR},SAMPLE=${SAMPLE},BAM=${BAM} \
--output=${LOGDIR}/dedup_${SAMPLE}.out \
--error=${LOGDIR}/dedup_${SAMPLE}.err << 'EOF'
#!/bin/bash
#SBATCH -p fast
#SBATCH -c 2
#SBATCH -N 1 
#SBATCH --mem-per-cpu 6000
#SBATCH -t 03:00:00

# Execute Picard MarkDuplicates via Singularity
singularity exec -B ${ROOT}:${ROOT} ${IMAGES}/picard_v2.27.4.simg \
picard MarkDuplicates \
I=${ALIGN}/${BAM} \
O=${OUTDIR}/${SAMPLE}.dedup.bam \
M=${OUTDIR}/${SAMPLE}.metrics.txt \
VALIDATION_STRINGENCY=LENIENT \
REMOVE_DUPLICATES=false \
CREATE_INDEX=true
EOF

done
