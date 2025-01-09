import sys
from Bio import Entrez
import os

# Set email, required by NCBI
Entrez.email = "#############"

# Function to search for the taxID of a species
def search_taxid(species):
    print(f"Searching taxid for: {species}")
    handle = Entrez.esearch(db="taxonomy", term=species)
    record = Entrez.read(handle)
    if record['IdList']:
        return record['IdList'][0]  # The first ID is the taxID
    else:
        return None

# Read the species list from a file
def load_species_from_file(filename):
    with open(filename, "r") as file:
        species_list = [line.strip() for line in file.readlines()]
    return species_list

# Save taxIDs to a file
def save_taxids_to_file(taxids, output_file):
    with open(output_file, "w") as file:
        for taxid in taxids:
            if taxid:
                file.write(f"{taxid}\n")

# Check if the correct number of command-line arguments were passed
if len(sys.argv) != 3:
    print("Usage: python fetch_taxids.py <species_file> <db_path>")
    sys.exit(1)

# Get the command-line arguments
species_file = sys.argv[1]
db_path = sys.argv[2]

# Check if the db_path exists, if not, create it
if not os.path.exists(db_path):
    os.makedirs(db_path)

# Define the output file path using db_path
output_file = os.path.join(db_path, "species_taxid_list.txt")

# Load species from the file
species_list = load_species_from_file(species_file)

# Search for taxIDs and save them
taxids = []
for species in species_list:
    taxid = search_taxid(species)
    if taxid:
        print(f"TaxID found: {taxid}")
        taxids.append(taxid)
    else:
        print(f"TaxID not found for: {species}")

# Save the taxIDs to a file
save_taxids_to_file(taxids, output_file)

print(f"TaxIDs have been saved to {output_file}")
