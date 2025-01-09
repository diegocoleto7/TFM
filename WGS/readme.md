# Workflow Details

1. **Host DNA Filtering**
   
   - **Script:** `host-filter.sh`
   - **Purpose:** Filters out host DNA from the samples in both datasets.
   - **Requirements:**
     - **GRCh38 Genome:** Download and index the GRCh38 genome using Bowtie2.
     - **Bowtie2 Indexing:**
       ```bash
       bowtie2-build GRCh38.fasta GRCh38_index
       ```


2. **Prokaryotic Sequence Filtering**
   
   - **Script:** `bac-filter.sh`
   - **Purpose:** Filters prokaryotic sequences from the samples in both datasets.
   - **Requirements:**
     - **Unified Human Gastrointestinal Genome (UHGG) v2.0.2:** Download and index the UHGG version 2.0.2 using Bowtie2.
     - **Download Link:** [UHGG v2.0.2 Library](https://ftp.ebi.ac.uk/pub/databases/metagenomics/mgnify_genomes/human-gut/v2.0.2/kraken2_db_uhgg_v2.0.2/library/library.fna)
     - **Bowtie2 Indexing:**
       ```bash
       bowtie2-build library.fna UHGG_v2.0.2_index
       ```


3. **Main Analysis Script**
   
   - **Script:** `kraken2.sh`
   - **Purpose:** Performs the taxonomic classification using Kraken2.


4. **Transforming Kraken2 Reports to BIOM Format**
   
   - **Command:**
     ```bash
     kraken-biom "{input_report_dir}" -o "{output_file}" --fmt json
     ```
   - **Purpose:** Converts Kraken2 reports into `.biom` format for downstream analysis.

5. **Post-Analysis and Graph Generation**
   
   - **Scripts:**
     - `AI4food-WGS.R`: Post-analysis script for the AI4food dataset.
     - `ECNR-WGS.R`: Post-analysis script for the ECNR dataset.
   - **Purpose:** Perform statistical analyses and generate visualizations for both WGS datasets.
