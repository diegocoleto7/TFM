#!/bin/bash

# Display help
Help() {
    echo "Usage: script.sh -p <path to the database> -f <species_file.txt> -t <number of threads>"
    exit 1
}

# Variables
while getopts "p:f:t:h" flag; do
    case "$flag" in
        p) dbpath=$OPTARG;; # path to db
        f) file_in=$OPTARG;; # path of species_file.txt
        t) threads=$OPTARG;; # number of threads
        h) Help;;
        *) Help;;
    esac
done

# Check that all necessary variables are defined
if [[ -z "$dbpath" || -z "$file_in" || -z "$threads" ]]; then
    echo "Missing required parameters"
    Help
fi

mkdir -p "${dbpath}"

#Get the species TaxID
python speciestoid.py $file_in $dbpath

output_file="$dbpath/species_taxid_list.txt"

#Download the sequences
python download_refseq.py $output_file $dbpath

fasta_output_file="${dbpath}/refseq_sequences.fasta"

# Apply dustmasker to the FASTA sequence file
output_fasta_masked="${dbpath}/masked_sequences.fasta"
dustmasker -in "$fasta_output_file" -out "$output_fasta_masked" -outfmt fasta

echo "Repetitive sequences have been masked..."

# Download the taxonomy
kraken2-build --download-taxonomy --db "$dbpath" --threads "$threads"
echo "Taxonomy downloaded..."

# Add masked sequences to the database
kraken2-build --add-to-library "$output_fasta_masked" --db "$dbpath"
kraken2-build --build --db "$dbpath" --threads "$threads"

echo "Kraken2 database created in $dbpath"
