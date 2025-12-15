#!/bin/bash
#SBATCH -p long
#SBATCH -c 2
#SBATCH -N 1
#SBATCH --mem-per-cpu 6000
#SBATCH -t 03:00:00
#SBATCH -o add_rg_launcher.out
#SBATCH -e add_rg_launcher.err

# Paths (Anonymized)
ROOT="/path/to/project_root" # Base path for project data and singularity images
BASEDIR="${ROOT}/project_name" # Main project directory (Anonymized)
DEDUP="${BASEDIR}/deduplicated" # Input directory containing deduplicated BAM files
IMAGES="${ROOT}/images" # Directory containing Singularity container images
OUTDIR="${BASEDIR}/deduplicated_rg" # Output directory for BAM files with Read Groups
LOGDIR="${BASEDIR}/logs/add_rg" # Directory for SLURM job logs

# Create output and log directories if they do not exist
mkdir -p ${OUTDIR}
mkdir -p ${LOGDIR}

cd ${DEDUP}

# Iterate over deduplicated BAM files
for BAM in *.dedup.bam; do
SAMPLE=$(basename ${BAM} .dedup.bam)

echo "Adding Read Group information to sample: ${SAMPLE}"

# Submit SLURM job
sbatch --export=ALL,ROOT=${ROOT},DEDUP=${DEDUP},IMAGES=${IMAGES},OUTDIR=${OUTDIR},SAMPLE=${SAMPLE},BAM=${BAM} \
--output=${LOGDIR}/add_rg_${SAMPLE}.out \
--error=${LOGDIR}/add_rg_${SAMPLE}.err << 'EOF'
#!/bin/bash
#SBATCH -p fast
#SBATCH -c 2
#SBATCH -N 1
#SBATCH --mem-per-cpu 6000
#SBATCH -t 03:00:00

# Execute Picard AddOrReplaceReadGroups via Singularity
singularity exec -B ${ROOT}:${ROOT} ${IMAGES}/picard_v2.27.4.simg \
picard AddOrReplaceReadGroups \
I=${DEDUP}/${BAM} \
O=${OUTDIR}/${SAMPLE}.rg.bam \
RGID=${SAMPLE} \
RGLB=lib1 \
RGPL=illumina \
RGPU=unit1 \
RGSM=${SAMPLE} \
CREATE_INDEX=true \
VALIDATION_STRINGENCY=LENIENT
EOF

done
