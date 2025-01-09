# Workflow Details

1. **Main Analysis Script**
   
   - **Script:** `Qiime2.sh`
   - **Purpose:** Performs the main ITS analysis using QIIME 2.
   - **Requirements:**
     - **Preprocessed and Effective Samples:** Ensure that your samples are preprocessed and ready for analysis.
     - **Pre-trained Classifier:** Download the pre-trained classifier from [this link](https://github.com/colinbrislawn/unite-train/releases/download/v10.0-v04.04.2024-qiime2-2024.5/unite_ver10_99_04.04.2024-Q2-2024.5.qza).


2. **Post-Analysis and Graph Generation**
   
   - **Scripts:**
     - `AI4food-ITS.R`: Post-analysis script for the AI4food ITS dataset.
     - `ECNR-ITS.R`: Post-analysis script for the ECNR ITS dataset.
   - **Purpose:** Perform statistical analyses and generate visualizations for ITS datasets.
