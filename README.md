# FARmiR: Framework for Analysis and Refinement of miRNA profiles

FARmiR is a Snakemake-based pipeline for processing and analyzing miRNAs/isomiRs from small RNA sequencing data. It includes steps for quality control, adapter prediction and trimming, annotation and quantification with [isoMiRmap](https://github.com/TJU-CMC-Org/isoMiRmap), evaluation of miRNA candidates using non-coding RNA databases and expression thresholds, and refinement of miRNA profiles. The workflow supports high-performance computing environments with Slurm integration and provides structured outputs for downstream analysis and visuazliation with [AIMEE](https://github.com/jonahcullen/AIMEE).

<p align="center">
  <img src="figures/pipeline.png" width="700" height="240" title="farmir pipeline">
</p>

## Dependencies

- [Python](https://www.python.org/)
- [Mamba](https://github.com/mamba-org/mamba)
- [Snakemake](https://snakemake.readthedocs.io/)
- [Snakemake-Profiles](https://github.com/Snakemake-Profiles)
- [Apptainer/Singularity](https://apptainer.org/)

## Initial setup

**1. Install dependencies**

If you do not have `conda` already installed, the Snakemake developers **recommend** to install via [Mambaforge](https://github.com/conda-forge/miniforge#mambaforge).

**NOTE:** *In January of 2024, Snakamake underwent a major update to version 8 which introduced a number of breaking changes. Until these are addressed, we suggest using version 7 which can be specified using version number as below.*

```
# download Mambaforge installer (assuming Unix-like platform)
wget "https://github.com/conda-forge/miniforge/releases/latest/download/Mambaforge-$(uname)-$(uname -m).sh"
bash Mambaforge-$(uname)-$(uname -m).sh

# updata Mamba
mamba update mamba -c conda-forge

# create a Snakemake environment which includes all Snakemake dependencies in
# addition to the miscellaneous modules above
mamba create \
    -c conda-forge -c bioconda \
    -n snakem \ # name of the environment
    python=3.10 snakemake=7.19
```

Alternatively, if you are already familiar with `conda` and creating environments, it is suggested to install `mamba` in your base environment and use that to build your environment.

```
# install mamba
conda install -n base -c conda-forge mamba

# create a Snakemake environment
mamba create \
    -c conda-forge -c bioconda \
    -n snakemake \ # name of the environment
    python=3.10 snakemake=7.19
```

**2. Download the container**

Due to the size of the included reference genomes and index files (and depending on your internet speed) this should take ~5 minutes.

```
wget https://s3.msi.umn.edu/se-smallseq/public/sif/hoof.sif
```

**3. Clone this repo**

```
git clone git@github.com:jonahcullen/FARmiR.git
```

