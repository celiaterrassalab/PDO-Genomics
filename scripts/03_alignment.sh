#!/bin/bash
#SBATCH -p long
#SBATCH -c 4
#SBATCH -N 1 
#SBATCH --mem-per-cpu 5000
#SBATCH -t 06:00:00 
#SBATCH -o align_launcher.out
#SBATCH -e align_launcher.err

# Paths
ROOT=""
BASEDIR="${ROOT}/"
TRIMMED="${BASEDIR}/trimmed"
IMAGES="${ROOT}/images"
REF=".../hg19/release-75/Homo_sapiens.GRCh37.75.dna.primary_assembly.fa"
OUTDIR="${BASEDIR}/alignment"
LOGDIR="${BASEDIR}/logs/align"

mkdir -p ${OUTDIR}
mkdir -p ${LOGDIR}

cd ${TRIMMED}

# Detectar las muestras únicas ignorando .1/.2
for base in $(ls *_val_1.fq.gz | sed 's/\.1_val_1.fq.gz//' | sed 's/\.2_val_1.fq.gz//' | sort -u); do
    R1="${base}.1_val_1.fq.gz"
    R2="${base}.2_val_2.fq.gz"

    SAMPLE="${base}"  # Esto dejará por ejemplo: S000021_S14266Nr1

    echo "Lanzando alineamiento para: ${SAMPLE}"

    sbatch --export=ALL,ROOT=${ROOT},TRIMMED=${TRIMMED},IMAGES=${IMAGES},REF=${REF},OUTDIR=${OUTDIR},SAMPLE=${SAMPLE},R1=${R1},R2=${R2} \
      --output=${LOGDIR}/align_${SAMPLE}.out \
      --error=${LOGDIR}/align_${SAMPLE}.err << 'EOF'
#!/bin/bash
#SBATCH -p fast
#SBATCH -c 4
#SBATCH -N 1 
#SBATCH --mem-per-cpu 5000
#SBATCH -t 06:00:00

# Alineamiento con BWA + sort con samtools
singularity exec -B ${ROOT}:${ROOT} ${IMAGES}/bwa_v0.7.17.sif \
  bwa mem -t 4 ${REF} ${TRIMMED}/${R1} ${TRIMMED}/${R2} | \
  singularity exec -B ${ROOT}:${ROOT} ${IMAGES}/samtools_v1.15.sif \
  samtools sort -@ 4 -o ${OUTDIR}/${SAMPLE}.sorted.bam -

# Indexado del BAM
singularity exec -B ${ROOT}:${ROOT} ${IMAGES}/samtools_v1.15.sif \
  samtools index ${OUTDIR}/${SAMPLE}.sorted.bam
EOF

done
