#!/bin/bash
#SBATCH -p bioinfo                      #partition/queue name
#SBATCH --job-name=kraken2_no                #Job name
#SBATCH -N 1                   #Un nodo
#SBATCH -n 10                   #Cores
#SBATCH -t 03:00:00             #Time limit hrs:min:sec
#SBATCH --output=/home/compubio/log-%j.o #Log de salida
#SBATCH --error=/home/compubio/log-%j.e #Log de errores
#SBATCH --mail-user=
#SBATCH --mail-type=ALL

# Module load
module load kraken2/2.1.3

threads=$SLURM_NTASKS

# Directories
input_dir="/path/to/data/nobac"
output_dir="/path/to/out/nobac"
kraken_output_dir="${output_dir}/kraken_outputs_confidence"
kraken_report_dir="${output_dir}/kraken_reports_confidence"
fungidb="/path/to/fungi_db"
# Create directories
mkdir -p $output_dir
mkdir -p $kraken_output_dir
mkdir -p $kraken_report_dir

# Initialize an array to hold sample names
samples=()

# Process each file to get sample names
for file in "${input_dir}"/*_nobac.fastq.1.gz; do
    base=$(basename "$file")
    sample=${base%_nobac.fastq.1.gz}
    samples+=("$sample")
done

# Iterate over each sample and process with Kraken2
for sample in "${samples[@]}"; do
    # Define the input paired-end files
    read1="${input_dir}/${sample}_nobac.fastq.1.gz"
    read2="${input_dir}/${sample}_nobac.fastq.2.gz"
    
    # Define the output files
    kraken_output="${kraken_output_dir}/${sample}_kraken_output.txt"
    kraken_report="${kraken_report_dir}/${sample}.txt"

    # Run Kraken2
    kraken2 --db "$fungidb" --paired --confidence 0.1 --output "$kraken_output" --report "$kraken_report" "$read1" "$read2" --threads "$threads"
done
