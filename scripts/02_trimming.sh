#!/bin/bash
#SBATCH -p long
#SBATCH -c 2
#SBATCH -N 1 
#SBATCH --mem-per-cpu 5000
#SBATCH -t 05:00:00 
#SBATCH -o trim_launcher.out
#SBATCH -e trim_launcher.err

# Variables (Anonymized)
ROOT="/path/to/project_root" # Base path for project data and singularity images
BASEDIR="${ROOT}/project_name" # Main project directory (Anonymized)
RAWDATA="${BASEDIR}/input_data" # Directory containing raw FASTQ files
IMAGES="${ROOT}/images" # Directory containing Singularity container images
OUTPUT_DIR="${BASEDIR}/trimmed" # Output directory for trimmed FASTQ files

# Create output directory if it does not exist
mkdir -p ${OUTPUT_DIR}

cd ${RAWDATA}

# Detect paired-end samples (based on *.1.fastq.gz)
for R1 in *.1.fastq.gz; do
SAMPLE=$(basename ${R1} .1.fastq.gz)
R2="${SAMPLE}.2.fastq.gz"

echo "Launching trimming job for sample: ${SAMPLE}"

sbatch --export=ROOT=${ROOT},IMAGES=${IMAGES},RAWDATA=${RAWDATA},OUTPUT_DIR=${OUTPUT_DIR},SAMPLE=${SAMPLE},R1=${R1},R2=${R2} << 'EOF'
#!/bin/bash
#SBATCH -p fast
#SBATCH -c 2
#SBATCH -N 1 
#SBATCH --mem-per-cpu 5000
#SBATCH -t 05:00:00 
#SBATCH -o trim_%j.out
#SBATCH -e trim_%j.err

# Execute Trim Galore via Singularity to remove adapters and low quality bases
singularity exec -B ${ROOT}:${ROOT} ${IMAGES}/trimgalore_v0.6.6.sif trim_galore \
-q 30 --paired -o ${OUTPUT_DIR} ${RAWDATA}/${R1} ${RAWDATA}/${R2}
EOF

done
