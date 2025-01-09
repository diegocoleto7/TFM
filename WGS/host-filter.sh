#!/bin/bash
#SBATCH -p bioinfo                      #partition/queue name
#SBATCH --job-name=Decont-host                #Job name
#SBATCH -N 1                   #Un nodo
#SBATCH -n 20                   #Cores
#SBATCH -t 48:00:00             #Time limit hrs:min:sec
#SBATCH --output=/home/compubio/log-%j.o #Log de salida
#SBATCH --error=/home/compubio/log-%j.e #Log de errores
#SBATCH --mail-user=
#SBATCH --mail-type=ALL

# Module load
module load bowtie2/2.5.4

# Directories
input_dir="/path/to/data/WGS"  
output_dir="/path/to/data/filtered_data"  
genome_dir="/path/to/genome"
mkdir -p ${output_dir}

# Filter FASTQ´s
for fastq_file in ${input_dir}/*.fastq.gz; do
    base_name=$(basename ${fastq_file} .fastq.gz)
    output_fastq="${output_dir}/${base_name}_filtered.fastq"
    bowtie2 -x ${genome_dir}/GRCh38_latest.bowtie-index -U ${fastq_file} --un-gz ${output_fastq}.gz -S /dev/null
done
