# CRISPRme
[![install with bioconda](https://img.shields.io/badge/install%20with-bioconda-brightgreen.svg?style=flat)](http://bioconda.github.io/recipes/crispritz/README.html)

CRISPRme is a web based tool dedicated to perform predictive analysis and result assessement on CRISPR/Cas experiments with a user-friendly GUI and the precise scope of searching individual variant in VCF dateset.

With this aim in mind we created a simple package that takes care of any step, from downloading the necessary data, to execute a complete search and present to the user an exhaustive result report with images and tabulated targets to navigate with the included web-based GUI.

The software is composed of two principal function:

- complete-search, function to perform a search starting from scratch input data, like a reference genome, a set of VCF data and list of sgRNAs targets.
- targets-integration, function to perform the integration of results targets with a gencode proximity file to identify genes near targets and collect only the top scoring targets in term of CFD score.
- web-interface, function to start the web-interface accessible locally from a web browser

# CRISPRme Installation and Usage
The two fastest way to use CRISPRitz is through the installation of Docker or Conda.
Here we summarize the steps to install CRISPRitz with Docker and Conda.

## Installation (Phase 1)
**Conda installation (Linux and MacOS):**
- Open a terminal window
- Paste this command into the terminal (Linux):
    ```
    curl https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh --output Miniconda3-latest-Linux-x86_64.sh
    ```
- Paste this command into the terminal (MacOS):
    ```
    curl https://repo.anaconda.com/miniconda/Miniconda3-latest-MacOSX-x86_64.sh --output Miniconda3-latest-MacOSX-x86_64.sh
    ```
- If the file is correctly downloaded you now need to execute it to complete the installation, so paste this command into the terminal:
    - Linux
    ```
    bash Miniconda3-latest-Linux-x86_64.sh
    ```
    - MacOS
    ```    
    bash Miniconda3-latest-MacOSX-x86_64.sh
    ```
- Press ENTER when requested and yes when an answer is requested, in this way you allow conda to set all the directories in your HOME path for an easy use
- After the complete installation you will receive this message “Thank you for installing Miniconda3!” to certify the correct installation.
- Now you need to close the terminal window you are using and open a new one, to allow the system to start conda.
- In the new terminal window you should see something like this:
    ```
    (base) user@nameofPC:~$
    ```
    If you read the "(base)" like this, conda is loaded correctly and you can start using it.
- Now you need to set the channels to allow conda to access different repositories and set the default version of python to version 3.8, so paste these commands into the terminal you just opened:
    ```
    conda config --add channels defaults
    conda config --add channels bioconda
    conda config --add channels conda-forge
    conda install python=3.8
    ```
- Now, you can install CRISPRme by typing the command:
    ```
    conda install crisprme -y
    ```
- To test your installation, type the command:
    ```
    crisprme.py
    ```
- After the execution of the command you should see a list of CRISPRitz tools.

Now the software is installed and ready to be used.

**Docker installation:  
Note: if you are using MasOS or Windows, you just need to download the installer file
and follow the on screen instructions.  
https://docs.docker.com/docker-for-windows/install/ (Windows)  
https://docs.docker.com/docker-for-mac/install/ (MacOS)**

**Ubuntu installation guide:**
- Open a terminal window
- Paste this command to update the index repository:
    ```
    sudo apt-get update
    ```
- Paste this command to allow package installation over HTTPS:
    ```
    sudo apt-get install \
    apt-transport-https \
    ca-certificates \
    curl \
    gnupg-agent \
    software-properties-common
    ```
- Paste this command to add the docker key:
    ```
    curl -fsSL https://download.docker.com/linux/ubuntu/gpg | sudo apt-key add -
    ```
- Paste this command to set the correct version of docker for your system:
    ```
    sudo add-apt-repository \
    "deb [arch=amd64] https://download.docker.com/linux/ubuntu \
    $(lsb_release -cs) \
    Stable"
    ```
- Paste this command to update the index repository another time, to make sure everything is ready and set to install docker:
    ```
    sudo apt-get update
    ```
- Then paste this command to finally install docker:
    ```
    sudo apt-get install docker-ce docker-ce-cli containerd.io
    ```
- Paste this last command to check if the installation is complete and functional:
    ```
    sudo docker run hello-world
    ```
- If this message is printed, everything is perfectly installed
![docker hello world](https://user-images.githubusercontent.com/40895152/63214349-769e3880-c117-11e9-8ee2-d754096b3aca.png)
- Now, we need to do some more steps to complete the settings. Paste this command to create a user group for docker user:
    ```
    sudo groupadd docker
    ```
- Paste this command to add your current user to the created group:
    ```
    sudo usermod -aG docker $USER
    ```
- Now you need to restart your machine or the virtual environment, to re-evaluate the user groups.
- One last command to test if the group is well configured. Paste this command:
    ```
    docker run hello-world
    ```
- If the previous “hello from docker” message is printed, everything is perfectly set.

## Post installation test (Phase <a name="phase2">2</a>):
- Download and run this package if you have installed CRISPRme or Docker with Conda:
    ```
    https://drive.google.com/file/d/11wn9pg6AWzDYZ7V_LjBIjGvgx95bnVJ1/view?usp=sharing
    ```
- Write this command to execute the script if you installed with Conda:
    ```
    bash crisprme_auto_test_conda.sh
    ```
- Write this command to execute the script if you installed with Docker:
    ```
    bash crisprme_auto_test_docker.sh
    ```
- You will see the processing starting, first with the download of all the necessary data and then with the analysis, depending on the system hardware and internet connection this may take very different time

## Usage (Phase 3):
**<a name="Add-Variant">3.1</a> CRISPRme complete-search function**  
This tool is created to perform a complete search from scratch.  
Input:
- Directory containing a genome in fasta format, need to be separated into single
chromosome files.
- Text file containing path to VCF directories [OPTIONAL]
- Text file with a list of guide (1 to N)
- Text file with a single PAM sequence
- Bed file with annotation, containing a list of genetic regions with a function association
- Text file containing a list of path to samplesID file (1 to N) equal to the number of VCF dataset used [OPTIONAL]
- Bed file extracted from Gencode data to find gene proximity of targets
- Maximal number of allowed bulges of any kind to compute in the index genome
- Threshold of mismatches allowed
- Size of DNA bulges allowed
- Size of RNA bulges allowed
- Merge range, necessary to reduce the inflation of targets due to bulges, it's the window of bp necessary to merge one target into another maintaining the highest scoring one
- Output directory, in which all the data will be produced
- Number of threads to use in computation
Output:
- bestMerge targets file, containing all the highest scoring targets, in terms of CFD and targets with the lowest combination of mismatches and bulges (with preference to lowest mismatches count), each genomic position is represent by one target
- altMerge targets file, containing all the discarded targets from the bestMerge file, each genomic position can be represented by more than target
- Parameters data file, containing all the parameters used in the search
- Count and distribution files, containing all the data count file useful in the web-tool representation to generate main tables and view
- Integrated results and database, containing all the tabulated data with genes proximity analysis and database representation to rapid querying on the web-tool GUI
- Directory with raw targets, containing the un-processed results from the search, useful to recover any possible target found during the search
- Directory with images, containing all the images generated to be used with the web-tool

Example call:
- Conda
    ```
    crisprme.py complete-search --genome Genomes/hg38/ --vcf list_vcf.txt/ --guide sg1617.txt --pam PAMs/20bp-NGG-spCas9.txt --annotation Annotations/gencode_encode.hg38.bed --samplesID list_samplesID.txt --gene_annotation Gencode/gencode.protein_coding.bed --bMax 2 --mm 6 --bDNA 2 --bRNA 2 --merge 3 --output Results/sg1617/ --thread 4
    ```
- Docker
    ```
    docker run -v ${PWD}:/DATA -w /DATA -i scancellieri/crisprme crisprme.py complete-search --genome Genomes/hg38/ --vcf list_vcf.txt/ --guide sg1617.txt --pam PAMs/20bp-NGG-spCas9.txt --annotation Annotations/gencode_encode.hg38.bed --samplesID list_samplesID.txt --gene_annotation Gencode/gencode.protein_coding.bed --bMax 2 --mm 6 --bDNA 2 --bRNA 2 --merge 3 --output Results/sg1617/ --thread 4
    ```

**<a name="Index-Genome">3.2</a> CRISPRitz Index-Genome Tool**  
This tool is created to generate an index genome (similar to the bwa-index step). This step is
time consuming (from 30 to 60 minutes) but helps to save a lot of execution time while
searching with lot of guides and with the support of bulges (DNA and RNA). If do not need to
search with bulges, skip this passage.  
Input:
- Name of the genome to create (e.g. `hg19_ref`).
- Directory containing a genome in fasta format, need to be separated into single
chromosome files.
- Text file containing the PAM (including a number of Ns equal to the guide length) and a
space separated number indicating the length of the PAM sequence (e.g. Cas9 PAM is
NNNNNNNNNNNNNNNNNNNNNGG 3). The sequence is composed by 20 Ns and
NGG, followed by number 3, representing the length of the PAM sequence.
- Number of bulges to include in the database to perform the following search (i.e. the max
number bulges allowed for DNA and RNA when searching on the database)
- Number of threads to use for the analysis (Optional)

Output:
- Directory containing an index genome in .bin format, separated into single chromosome
files, containing all the candidate targets for a selected PAM, adding also characters to
perform bulge search.

Example call:
- Conda
    ```
    crispritz.py index-genome hg19_ref hg19_ref/ pam/pamNGG.txt -bMax 2 -th 4
    ```
- Docker
    ```
    docker run -v ${PWD}:/DATA -w /DATA -i pinellolab/crispritz crispritz.py index-genome hg19_ref hg19_ref/ pam/pamNGG.txt -bMax 2 -th 4
    ```
