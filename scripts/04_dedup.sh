#!/bin/bash
#SBATCH -p long
#SBATCH -c 2
#SBATCH -N 1 
#SBATCH --mem-per-cpu 6000
#SBATCH -t 03:00:00 
#SBATCH -o dedup_launcher.out
#SBATCH -e dedup_launcher.err

# Paths
ROOT=""
BASEDIR="${ROOT}/"
ALIGN="${BASEDIR}/alignment"
IMAGES="${ROOT}/images"
OUTDIR="${BASEDIR}/deduplicated"
LOGDIR="${BASEDIR}/logs/dedup"

mkdir -p ${OUTDIR}
mkdir -p ${LOGDIR}

cd ${ALIGN}

# Iterar sobre BAMs alineados
for BAM in *.sorted.bam; do
SAMPLE=$(basename ${BAM} .sorted.bam)

echo "Lanzando Picard MarkDuplicates para: ${SAMPLE}"

sbatch --export=ALL,ROOT=${ROOT},ALIGN=${ALIGN},IMAGES=${IMAGES},OUTDIR=${OUTDIR},SAMPLE=${SAMPLE},BAM=${BAM} \
--output=${LOGDIR}/dedup_${SAMPLE}.out \
--error=${LOGDIR}/dedup_${SAMPLE}.err << 'EOF'
#!/bin/bash
#SBATCH -p fast
#SBATCH -c 2
#SBATCH -N 1 
#SBATCH --mem-per-cpu 6000
#SBATCH -t 03:00:00

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
