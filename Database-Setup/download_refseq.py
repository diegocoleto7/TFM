import sys
from Bio import Entrez, SeqIO
import os

#Mail required
Entrez.email = "###############"

# Function to download RefSeq sequences by taxID
def download_refseq_by_taxid(taxid_list, output_file):
    with open(output_file, "w") as out_handle:
        for taxid in taxid_list:
            taxid = taxid.strip() 
            print(f"Downloading sequences for taxID: {taxid}")

            # Search for nucleotide sequences in RefSeq for the given taxID
            handle = Entrez.esearch(db="nucleotide", term=f"txid{taxid}[Organism] AND refseq[filter]", retmax=10000)
            record = Entrez.read(handle)
            ids = record["IdList"]
            if ids:
                # Fetch the sequences in FASTA format
                handle = Entrez.efetch(db="nucleotide", id=",".join(ids), rettype="fasta", retmode="text")
                out_handle.write(handle.read())

# Function to load taxIDs from a file
def load_taxids_from_file(filename):
    with open(filename, "r") as file:
        taxids = file.readlines()  
    return taxids


# Get the command-line arguments
taxid_file = sys.argv[1]
db_path = sys.argv[2]



# Define the output file path using db_path
output_file = os.path.join(db_path, "refseq_sequences.fasta")

# Load taxIDs from the file
taxids = load_taxids_from_file(taxid_file)

# Download sequences and save them to the output file
download_refseq_by_taxid(taxids, output_file)

print(f"Sequences have been saved to {output_file}")
