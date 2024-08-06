[Installation](#installation) | [Usage](#usage) | [Citation](#citation)
# Snakemake workflow to process the Illumina reads for Methylated Fraction Analysis

This software is a [snakemake](https://snakemake.readthedocs.io/en/stable/index.html) workflow that processes the raw fastq files from Illumina sequencing. It first generated SLURM scripts to 1) align the raw fastq files to reference genome using Bowtie2; 2) to calculate fragment coverage and counts of read 5' ends using Bedtools. 

After running the first step on SLURM, the workflow will process the raw coverage and read counts to obtain the methylated fraction at GATC sites to measure the DAM methylation (cut by DpnI). The final result will be stored both in a SQLite database for all the samples and in separate BigWig files for individual samples.

Two metadata files are required to execute this script: 1) a sample table that contains the prefix of the fastq files (the pair-end reads must share the same prefix) and the label of the samples. The exmpale sample table is in `resources/example_sample_table.csv`. 2) a configuration file that contains all the other required information. The example configuration file is in `config/example_config.yaml`.

The script provides the chromosome size table for T2T v1.1 and Hg38 assemlby in `resources/`. The three columns are Chromosome Name (chosen by users for downstream analysis), Sequence Name (in reference fasta files), and the size of the chromosome.

The workflow will generate one SQLite data for all samples and BigWig files for each sample. The SQLite database is stored in {gatc_db} specified in the configuration file. The BigWig files are stored in the {bw_prefix}.Filter_{sample-name}.bw specified in the configuration file. 

The BigWig files could be processed by [methylFracAnalyzer](https://github.com/zhuweix/methylFracAnalyzer) for downstream analysis.

## Installation

On the SLURM server, the slurm scripts required that Bowtie 2 and bedtools are installed and could be loaded by module load. The scripts are generated based on [NIH Biowulf HPC](https://hpc.nih.gov/). The users could also modify the `workflow/scripts/make_slurm_script.py` to modify the master script on their own server.

To install this package:
```
# First clone this repository
git clone https://github.com/zhuweix/snakemakeMethylFrac.git


# Then install the package
cd snakemakeMethylFrac

pip install .
```


## Usage

This workflow could be executed in two steps.

1. Peform the alignment and coverage calculation on SLURM
```
snakemake -c 1 run_slurm_script --configfile {path/to/config-file.yaml}
```
For the first time to run this script, it is recommended to check the SLURM scripts first:
```
snakemake -c 1 make_slurm_script --configfile {path/to/config-file.yaml}
```
It is also recommended to perform a dry-run with `-np` argument:
```
snakemake -np run_slurm_script --configfile {path/to/config-file.yaml}
```

The output of the bam files are stored in {output_dir} specified in the configuration file. The coverage and read counts are stored in the same directory as the bam files.

2. Calculate the Methylated Fractions
```
snakemake -c 1 --configfile {path/to/config-file.yaml}
```
Currently, this package doesn't support run with multiple cores. Using multiple cores (e.g. `-c 4`) will cause the script to fail.
We also recommend to perform a dry-run with `-np` argument:
```
snakemake -np --configfile {path/to/config-file.yaml}
```

The output SQLite database and BigWig files are stored in {gatc_db} and {bw_prefix}.Filter_{sample-name}.bw specified in the configuration file, respectively (sample-name is specificed in the sample table file, which is also specified in the configuration file).

Furthermore, this workflow also generated the histogram of Fragment Occupancy and Half-site Methylated Fraction (measured by DpnI Digestion) in `results/figures`.

The workflow could also be executed step by step. Please refer to te `workflow/Snakefile` for more details.

## Citation
TBD
