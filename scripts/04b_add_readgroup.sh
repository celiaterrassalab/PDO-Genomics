#!/bin/bash
#SBATCH -p long
#SBATCH -c 2
#SBATCH -N 1
#SBATCH --mem-per-cpu 6000
#SBATCH -t 03:00:00
#SBATCH -o add_rg_launcher.out
#SBATCH -e add_rg_launcher.err

# Paths
ROOT=""
BASEDIR="${ROOT}/"
DEDUP="${BASEDIR}/deduplicated"
IMAGES="${ROOT}/images"
OUTDIR="${BASEDIR}/deduplicated_rg"
LOGDIR="${BASEDIR}/logs/add_rg"

mkdir -p ${OUTDIR}
mkdir -p ${LOGDIR}

cd ${DEDUP}

for BAM in *.dedup.bam; do
SAMPLE=$(basename ${BAM} .dedup.bam)

echo "Añadiendo Read Group a: ${SAMPLE}"

sbatch --export=ALL,ROOT=${ROOT},DEDUP=${DEDUP},IMAGES=${IMAGES},OUTDIR=${OUTDIR},SAMPLE=${SAMPLE},BAM=${BAM} \
--output=${LOGDIR}/add_rg_${SAMPLE}.out \
--error=${LOGDIR}/add_rg_${SAMPLE}.err << 'EOF'
#!/bin/bash
#SBATCH -p fast
#SBATCH -c 2
#SBATCH -N 1
#SBATCH --mem-per-cpu 6000
#SBATCH -t 03:00:00

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
