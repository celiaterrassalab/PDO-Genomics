#!/bin/bash
#SBATCH -p long
#SBATCH -c 2
#SBATCH -N 1
#SBATCH --mem-per-cpu 5000
#SBATCH -t 02:00:00
#SBATCH -o fastqc_launcher.out
#SBATCH -e fastqc_launcher.err

# Paths
ROOT="/projects/cancer"
BASEDIR="${ROOT}/HelenaBrunel/Chiara_DNA_CEGAT"
RAWDATA="${BASEDIR}/input_data"
IMAGES="${ROOT}/images"
OUTPUT_DIR="${BASEDIR}/qc/fastqc_raw"

mkdir -p ${OUTPUT_DIR}

cd ${RAWDATA}

# Loop sobre todos los FASTQ y lanza un job por archivo
for fq in *.fastq.gz; do
echo "Launching FastQC for: ${fq}"

sbatch --export=ROOT=${ROOT},RAWDATA=${RAWDATA},IMAGES=${IMAGES},OUTPUT_DIR=${OUTPUT_DIR},FQFILE=${fq} << 'EOF'
#!/bin/bash
#SBATCH -p fast
#SBATCH -c 2
#SBATCH -N 1
#SBATCH --mem-per-cpu 5000
#SBATCH -t 01:00:00
#SBATCH -o fastqc_%j.out
#SBATCH -e fastqc_%j.err

singularity exec -B ${ROOT}:${ROOT} ${IMAGES}/fastqc_v0.11.9.sif \
fastqc ${RAWDATA}/${FQFILE} -o ${OUTPUT_DIR}
EOF

done
