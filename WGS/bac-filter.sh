#!/bin/bash
#SBATCH -p bioinfo                      #partition/queue name
#SBATCH --job-name=bac_decont                #Job name
#SBATCH -N 1                   #Un nodo
#SBATCH -n 12                   #Cores
#SBATCH -t 48:00:00             #Time limit hrs:min:sec
#SBATCH --output=/home/compubio/log-%j.o #Log de salida
#SBATCH --error=/home/compubio/log-%j.e #Log de errores
#SBATCH --mail-user=
#SBATCH --mail-type=ALL

export TMPDIR=/scratch/compubio

#Module load
module load bowtie2/2.5.4

# Directories
input_dir="/path/to/data/filtered_data"
output_dir="/path/to/data/nobac"
bacteria_db="/path/to/bacteria_db"

# Create output dir
mkdir -p "${output_dir}"

# Read files
for fastq_file1 in "${input_dir}"/*_1_filtered.fastq.gz; do
    base_name=$(basename "${fastq_file1}" _1_filtered.fastq.gz)
    fastq_file2="${input_dir}/${base_name}_2_filtered.fastq.gz"

    # Bowtie2 decontamination
    bowtie2 -x "${bacteria_db}/uhgg" \
        -1 "${fastq_file1}" \
        -2 "${fastq_file2}" \
        --un-conc-gz "${output_dir}/${base_name}_nobac.fastq.gz" \
        -S /dev/null \
        --threads 15
done
