#!/bin/bash
#SBATCH -p bioinfo                      #partition/queue name
#SBATCH --job-name=Genome-index                #Job name
#SBATCH -N 1                   #Un nodo
#SBATCH -n 15                   #Quince cores
#SBATCH -t 24:00:00             #Time limit hrs:min:sec
#SBATCH --output=/home/compubio/log-%j.o #Log de salida
#SBATCH --error=/home/compubio/log-%j.e #Log de errores
#SBATCH --mail-user=diegocoleto7@gmail.com
#SBATCH --mail-type=ALL

# Cargar el modulo
module load bowtie2/2.5.4

# Directorios
genome_dir="home/proyectos/imdeaalim/compubio/genome"  
mkdir -p ${genome_dir}

# Descargar e indexar el genoma del host.
wget -O ${genome_dir}/GRCh38_latest_genomic.fna.gz https://ftp.ncbi.nlm.nih.gov/refseq/H_sapiens/annotation/GRCh38_latest/refseq_identifiers/GRCh38_latest_genomic.fna.gz

gunzip -c ${genome_dir}/GRCh38_latest_genomic.fna.gz > ${genome_dir}/GRCh38_latest_genomic.fna

bowtie2-build ${genome_dir}/GRCh38_latest_genomic.fna.gz ${genome_dir}/GRCh38_latest.bowtie-index
