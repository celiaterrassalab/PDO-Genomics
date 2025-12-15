#!/bin/bash
#SBATCH -p long
#SBATCH -c 1
#SBATCH -N 1
#SBATCH --mem-per-cpu 4000
#SBATCH -t 02:00:00
#SBATCH -o cnv_launcher.out
#SBATCH -e cnv_launcher.err

# Paths
ROOT="/projects/cancer"
BASEDIR="${ROOT}/HelenaBrunel/Chiara_DNA_CEGAT"
DEDUP="${BASEDIR}/deduplicated_rg"
IMAGES="${ROOT}/images"
REF="${ROOT}/db_files/Genomes/Ensembl/human/hg19/release-75/Homo_sapiens.GRCh37.75.dna.primary_assembly.fa"
BED="${BASEDIR}/input_data/S000021_hg19_targets.bed"
OUTDIR="${BASEDIR}/cnv_gatk"
LOGDIR="${BASEDIR}/logs/gatk_cnv"

mkdir -p ${OUTDIR}
mkdir -p ${LOGDIR}

# Paso 1: Generar una vez los intervalos
if [[ ! -f ${OUTDIR}/preprocessed_intervals.interval_list ]]; then
    singularity exec -B ${ROOT}:${ROOT} ${IMAGES}/gatk4_v4.3.0.sif     gatk PreprocessIntervals         -R ${REF}         -L ${BED}         --bin-length 0         --interval-merging-rule OVERLAPPING_ONLY         -O ${OUTDIR}/preprocessed_intervals.interval_list

    singularity exec -B ${ROOT}:${ROOT} ${IMAGES}/gatk4_v4.3.0.sif     gatk AnnotateIntervals         -R ${REF}         -L ${OUTDIR}/preprocessed_intervals.interval_list         -O ${OUTDIR}/annotated_intervals.tsv
fi

# Paso 2: Lanzar un job por muestra
cd ${DEDUP}
for BAM in *.bam; do
    SAMPLE=$(basename ${BAM} .bam)

    echo "Enviando job para muestra: ${SAMPLE}"

    sbatch --export=ALL,ROOT=${ROOT},DEDUP=${DEDUP},IMAGES=${IMAGES},OUTDIR=${OUTDIR},SAMPLE=${SAMPLE},REF=${REF}       --output=${LOGDIR}/gatkcnv_${SAMPLE}.out       --error=${LOGDIR}/gatkcnv_${SAMPLE}.err << 'EOF'
#!/bin/bash
#SBATCH -p fast
#SBATCH -c 2
#SBATCH -N 1
#SBATCH --mem-per-cpu 5000
#SBATCH -t 04:00:00

# Step 1: CollectReadCounts
singularity exec -B ${ROOT}:${ROOT} ${IMAGES}/gatk4_v4.3.0.sif gatk CollectReadCounts     -I ${DEDUP}/${SAMPLE}.bam     -L ${OUTDIR}/preprocessed_intervals.interval_list     --interval-merging-rule OVERLAPPING_ONLY     -O ${OUTDIR}/${SAMPLE}.counts.hdf5

# Step 2: DenoiseReadCounts
singularity exec -B ${ROOT}:${ROOT} ${IMAGES}/gatk4_v4.3.0.sif gatk DenoiseReadCounts     -I ${OUTDIR}/${SAMPLE}.counts.hdf5     --standardized-copy-ratios ${OUTDIR}/${SAMPLE}.standardizedCR.tsv     --denoised-copy-ratios ${OUTDIR}/${SAMPLE}.denoisedCR.tsv

# Step 3: ModelSegments (sin -I)
singularity exec -B ${ROOT}:${ROOT} ${IMAGES}/gatk4_v4.3.0.sif gatk ModelSegments     --denoised-copy-ratios ${OUTDIR}/${SAMPLE}.denoisedCR.tsv     --output ${OUTDIR}/${SAMPLE}_segments     --output-prefix ${SAMPLE}

EOF

done
