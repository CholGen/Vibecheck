<p align="center">
    <picture>
      <source media="(prefers-color-scheme: dark)" srcset=".github/logo_light.png">
      <source media="(prefers-color-scheme: light)" srcset=".github/logo_dark.png">
      <img alt="bacpage" src=".github/logo_dark.png" width=500>
    </picture>
</p>

This repository contains an easy-to-use pipeline for the assembly and analysis of bacterial genomes using ONT long-read or Illumina short-read technology. 
Read the complete documentation and instructions for bacpage and each of its functions [here](https://cholgen.github.io/sequencing-resources/bacpage-command.html)

# Introduction
Advances in sequencing technology during the COVID-19 pandemic has led to massive increases in the generation of sequencing data. Many bioinformatics tools have been developed to analyze this data, but very few tools can be utilized by individuals without prior bioinformatics training.

This pipeline was designed to encapsulate pre-existing tools to automate analysis of whole genome sequencing of bacteria. 
Installation is fast and straightfoward. 
The pipeline is easy to setup and contains rationale defaults, but is highly modular and configurable by more advance users.
Bacpage has individual commands to generate consensus sequences, perform *de novo* assembly, construct phylogenetic tree, and generate quality control reports.

# Features
We anticipate the pipeline will be able to perform the following functions:
- [x] Reference-based assembly of Illumina paired-end reads
- [x] *De novo* assembly of Illumina paired-end reads
- [ ] *De novo* assembly of ONT long reads
- [x] Run quality control checks
- [x] Variant calling using [bcftools](https://github.com/samtools/bcftools)
- [x] Maximum-likelihood phylogenetic inference of processed samples and background dataset using [iqtree](https://github.com/iqtree/iqtree2) 
- [x] MLST profiling and virulence factor detection
- [x] Antimicrobial resistance genes detection
- [ ] Plasmid detection

# Installation
1. Install `mamba` by running the following two command:
```commandline
curl -L -O "https://github.com/conda-forge/miniforge/releases/latest/download/Mambaforge-$(uname)-$(uname -m).sh"
bash Mambaforge-$(uname)-$(uname -m).sh
```

2. Clone the bacpage repository:
```commandline
git clone https://github.com/CholGen/bacpage.git
```

3. Switch to the development branch of the pipeline:
```commandline
cd bacpage/
git checkout -b split_into_command
```

3. Install and activate the pipeline's conda environment:
```commandline
mamba env create -f environment.yaml
mamba activate bacpage
```

4. Install the `bacpage` command:
```commandline
pip install .
```

5. Test the installation:
```commandline
bacpage -h
bacpage version
```
These command should print the help and version of the program. Please create an issue if this is not the case.

# Updating

1. Navigate to the directory where you cloned the bacpage repository on the command line:
```commandline
cd bacpage/
```
2. Activate the bacpage conda environment:
```commandline
mamba activate bacpage
```
3. Pull the lastest changes from GitHub:
```commandline
git pull
```
4. Update the bacpage conda environemnt:
```commandline
mamba env update -f environment.yaml
```
5. Reinstall the `bacpage` command:
```commandline
pip install .
```

# Usage
0. Activate the bacpage conda environment:
```commandline
mamba activate bacpage
```
1. Create a directory specifically for the batch of samples you would like to analyze (called a project directory).
```commandline
bacpage setup [your-project-directory-name]
```
2. Place paired sequencing reads in the `input/` directory of your project directory.
3. From the pipeline's directory, run the reference-based assembly pipeline on your samples using the following command:
```commandline
bacpage assemble [your-project-directory-name]
```
This will generate a consensus sequence in FASTA format for each of your samples and place them in 
`<your-project-directory-name>/results/consensus_sequences/<sample>.masked.fasta`. An HTML report containing alignment 
and quality metrics for your samples can be found at `<your-project-directory-name>/results/reports/qc_report.html`.
