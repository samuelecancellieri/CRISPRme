# CRISPRme 
[![install with bioconda](https://img.shields.io/badge/install%20with-bioconda-brightgreen.svg?style=flat)](http://bioconda.github.io/recipes/crispritz/README.html)

CRISPRme is a web-based tool created to perform predictive analysis and result assessement on CRISPR/Cas experiments and visualizing everything in a simple and complete way using a web interface.

We exploited CRISPRtiz (Cancellieri, Samuele, et al. "CRISPRitz: rapid, high-throughput and variant-aware in silico off-target site identification for CRISPR genome editing." Bioinformatics 36.7 (2020): 2001-2008.) a powerful search engine we developed to fast search and account for targets and off-targets on any genome.
CRISPRme integrates CRISPRitz with a graphical interface completely dedicated to visualize and enlighten important targets and off-targets in reference and variant-constructed genome. 

# CRISPRme Installation and Usage
The two fastest way to use CRISPRme is through the installation of Docker or Conda.

##CRISPRme Conda installation

- In a conda environment, set the channels necessary to install all the packages and dependencies
    ```
    conda config --add channels defaults
    conda config --add channels bioconda
    conda config --add channels conda-forge
    conda install python=3.8
    ```
- Now, you can install CRISPRitz by typing the command:
    ```
    conda install crisprme -y
    ```
- To test your installation, type the command:
    ```
    crisprme.py
    ```
Now the software is installed and ready to be used.

## Post installation test
**Conda:**
- Download and run this script if you have installed CRISPRitz with Conda:
    ```
    curl https://raw.githubusercontent.com/pinellolab/CRISPRitz/master/test_scripts/auto_test_crispritz_conda.sh --output auto_test_crispritz_conda.sh
    ```
- Write this command to execute the script:
    ```
    bash auto_test_crispritz_conda.sh
    ```
- Wait until this confirmation message appears:  
“EVERY TEST PASSED!!! ENJOY CRISPRitz”

**Docker:**
- Download and run this script if you have installed CRISPRitz with Docker:
    ```
    curl https://raw.githubusercontent.com/pinellolab/CRISPRitz/master/test_scripts/auto_test_crispritz_docker.sh --output auto_test_crispritz_docker.sh
    ```
- Write this command to execute the script:
    ```
    bash auto_test_crispritz_docker.sh
    ```
- Wait until this confirmation message appears:  
“EVERY TEST PASSED!!! ENJOY CRISPRitz”

## Usage (Phase 3):
Here is a brief guide to help use CRISPRme, **if you already execute the post installation test
(Phase [2](#phase2)), and you obtain a positive result, you have all the necessary file in the
test_crispritz directory and you can skip this list of steps.**  
If you did not execute the test, follow these few steps to download the necessary files to try
CRISPRme.  
Download test files (ONE TIME STEP):
- The script will download the chr22 from UCSC (hg19), the correspondent VCF file from
the 1000 Genome Project, a directory containing some test guides, a directory
containing some PAM sequences and a directory of pre-computed genomic annotations
for the hg19 genome.
- Download the script with this command:
    ```
    curl https://raw.githubusercontent.com/pinellolab/CRISPRitz/master/test_scripts/download_test_files.sh --output download_test_files.sh
    ```
- Write this command to execute the script:
    ```
    bash download_test_files.sh
    ```
- The script will download every necessary file to test the software, we download only one
chromosome and one vcf file, to save time. All the examples can be run on an entire
genome, if you want to use the entire hg19 genome, you only need to add chromosomes
into the `hg19_ref` directory.
