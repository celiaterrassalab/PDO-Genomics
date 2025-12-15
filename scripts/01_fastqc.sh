#!/bin/bash
#SBATCH -p long
#SBATCH -c 2
#SBATCH -N 1
#SBATCH --mem-per-cpu 5000
#SBATCH -t 02:00:00
#SBATCH -o fastqc_launcher.out
#SBATCH -e fastqc_launcher.err

# Paths (Anonymized)
ROOT="/path/to/project_root" # Base path for project data and singularity images
BASEDIR="${ROOT}/project_name" # Main project directory (Anonymized)
RAWDATA="${BASEDIR}/input_data" # Directory containing raw FASTQ files
IMAGES="${ROOT}/images" # Directory containing Singularity container images
OUTPUT_DIR="${BASEDIR}/qc/fastqc_raw" # Output directory for raw FastQC reports

# Create output directory if it does not exist
mkdir -p ${OUTPUT_DIR}

cd ${RAWDATA}

# Loop over all FASTQ files and launch one SLURM job per file
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

# Execute FastQC via Singularity
singularity exec -B ${ROOT}:${ROOT} ${IMAGES}/fastqc_v0.11.9.sif \
fastqc ${RAWDATA}/${FQFILE} -o ${OUTPUT_DIR}
EOF

done
