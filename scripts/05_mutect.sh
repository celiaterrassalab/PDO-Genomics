#!/bin/bash
#SBATCH -p long
#SBATCH -c 4
#SBATCH -N 1 
#SBATCH --mem-per-cpu 6000
#SBATCH -t 06:00:00 
#SBATCH -o mutect2_launcher.out
#SBATCH -e mutect2_launcher.err

# Paths
ROOT="/projects/cancer"
BASEDIR="${ROOT}/HelenaBrunel/Chiara_DNA_CEGAT"
DEDUP="${BASEDIR}/deduplicated_rg"
IMAGES="${ROOT}/images"
REF="${ROOT}/db_files/Genomes/Ensembl/human/hg19/release-75/Homo_sapiens.GRCh37.75.dna.primary_assembly.fa"
BED="${BASEDIR}/input_data/S000021_hg19_targets.bed"
OUTDIR="${BASEDIR}/variants"
LOGDIR="${BASEDIR}/logs/mutect2"

mkdir -p ${OUTDIR}
mkdir -p ${LOGDIR}

cd ${DEDUP}

for BAM in *.bam; do
    FILE=${DEDUP}/${BAM}
    
    # Extraer nombre de muestra desde @RG línea (SM:)
    SAMPLE_RG=$(singularity exec -B ${ROOT}:${ROOT} ${IMAGES}/samtools_v1.15.sif \
        samtools view -H ${FILE} | grep '^@RG' | sed -n 's/.*SM:\([^\t]*\).*/\1/p' | head -n 1)

    # Verificación adicional
    if [[ -z "$SAMPLE_RG" ]]; then
        echo "❌ No se encontró SM: en ${BAM}, saltando muestra."
        continue
    fi

    echo "✅ Llamando variantes para: ${SAMPLE_RG} (archivo: ${BAM})"

    sbatch --export=ALL,ROOT=${ROOT},DEDUP=${DEDUP},IMAGES=${IMAGES},REF=${REF},BED=${BED},OUTDIR=${OUTDIR},BAM=${BAM},SAMPLE_RG=${SAMPLE_RG} \
      --output=${LOGDIR}/mutect2_${SAMPLE_RG}.out \
      --error=${LOGDIR}/mutect2_${SAMPLE_RG}.err << 'EOF'
#!/bin/bash
#SBATCH -p fast
#SBATCH -c 4
#SBATCH -N 1 
#SBATCH --mem-per-cpu 6000
#SBATCH -t 06:00:00

# Mutect2 tumor-only
singularity exec -B ${ROOT}:${ROOT} ${IMAGES}/gatk4_v4.3.0.sif \
gatk Mutect2 \
-R ${REF} \
-I ${DEDUP}/${BAM} \
-tumor ${SAMPLE_RG} \
-L ${BED} \
--native-pair-hmm-threads 4 \
-O ${OUTDIR}/${SAMPLE_RG}.raw.vcf.gz

# Filtrado de llamadas
singularity exec -B ${ROOT}:${ROOT} ${IMAGES}/gatk4_v4.3.0.sif \
gatk FilterMutectCalls \
-V ${OUTDIR}/${SAMPLE_RG}.raw.vcf.gz \
-R ${REF} \
-O ${OUTDIR}/${SAMPLE_RG}.filtered.vcf.gz
EOF

done
