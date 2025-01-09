#!/bin/bash
#SBATCH -p bioinfo                      #partition/queue name
#SBATCH --job-name=qiime2024.5              #Job name
#SBATCH -N 1                   #Un nodo
#SBATCH -n 50                   #Cores
#SBATCH -t 100:00:00             #Time limit hrs:min:sec
#SBATCH --output=/home/compubio/log-%j.o #Log de salida
#SBATCH --error=/home/compubio/log-%j.e #Log de errores
#SBATCH --mail-user=··········
#SBATCH --mail-type=END,FAIL

# Module load
module load qiime/amplicon-2024.5


# Temporal dir configuration
export TMPDIR=/scratch/$USER/$SLURM_JOB_ID

# Output dir
output_dir="/path/to/qiime2"

# Import merged data
qiime tools import \
  --type 'SampleData[SequencesWithQuality]' \
  --input-path /path/to/manifest.tsv \
  --output-path ${output_dir}/merged-demux.qza \
  --input-format SingleEndFastqManifestPhred33V2

# Visualization of demux
qiime demux summarize \
  --i-data ${output_dir}/merged-demux.qza \
  --o-visualization ${output_dir}/demux.qzv

# Dereplication of  sequences
qiime vsearch dereplicate-sequences \
  --i-sequences ${output_dir}/merged-demux.qza \
  --o-dereplicated-table ${output_dir}/table.qza \
  --o-dereplicated-sequences ${output_dir}/rep-seqs.qza
  --threads 50
# Clasification  UNITE trained

qiime feature-classifier classify-sklearn \
  --i-classifier /path/to/unite_ver10_99_04.04.2024-Q2-2024.5.qza \
  --i-reads ${output_dir}/rep-seqs.qza \
  --o-classification ${output_dir}/taxonomy.qza
  --p-n-jobs 50

# Export feature table
qiime tools export \
  --input-path ${output_dir}/table.qza \
  --output-path ${output_dir}/exported-feature-table

# Convert feature table to biom format
biom convert \
  -i ${output_dir}/exported-feature-table/feature-table.biom \
  -o ${output_dir}/exported-feature-table/feature-table.tsv \
  --to-tsv
  
# Export representative sequences
qiime tools export \
  --input-path ${output_dir}/rep-seqs.qza \
  --output-path ${output_dir}/exported-rep-seqs

# Export taxonomy
qiime tools export \
  --input-path ${output_dir}/taxonomy.qza \
  --output-path ${output_dir}/exported-taxonomy
