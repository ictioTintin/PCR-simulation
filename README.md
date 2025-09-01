# In silico design and validation of fish-specific primers

>Versión en español aquí:\
>[PCR-simulation](https://github.com/ictioTintin/PCR-simulation)

**Author:** Martin Holguin Osorio  
**Date:** January 2020  

> These experiments are based on the work by Eric Coissac (2015):  
> [Designing a new DNA metabarcode for fish](https://metabarcoding.org/IMG/html/primerdesign.html)

---

## Overview

This repository contains scripts and instructions to design and validate fish-specific primers **in silico** using mitochondrial genomes, ecoPCR, and OBITools.  
The workflow includes:  

1. Downloading mitochondrial genome sequences and taxonomy from NCBI.  
2. Formatting the data for OBITools and ecoPCR.  
3. Designing primer pairs for **Teleostei** using *ecoPrimers*.  
4. Testing primer specificity with simulated PCRs (*ecoPCR*).  
5. Evaluating primer performance and taxonomic resolution in **R** with ROBITools.  

---

## Dependencies

To run the pipeline, you need:  

- A computer with a **UNIX-based OS** (Linux, macOS) or a server/VM capable of running command-line tools.  
- Installed software:  
  - [OBITools](https://pythonhosted.org/OBITools/welcome.html#installing-the-obitools)  
  - [ecoPCR](https://git.metabarcoding.org/obitools/ecopcr/-/wikis/home)  
  - [ecoPrimers](https://git.metabarcoding.org/obitools/ecocebadores/-/wikis/home)  
  - [R Project](https://www.r-project.org/) with packages: **ROBITools, ROBITaxonomy, ROBIBarcodes**  

---

## Usage

### Download data

```bash
mkdir mitochondria2
cd mitochondria2

# Download mitochondrial genomes
wget "ftp://ftp.ncbi.nlm.nih.gov/genomes/refseq/mitochondrion/mitochondrion.1.genomic.gbff.gz"
wget "ftp://ftp.ncbi.nlm.nih.gov/genomes/refseq/mitochondrion/mitochondrion.2.genomic.gbff.gz"

# Download taxonomy
wget "ftp://ftp.ncbi.nlm.nih.gov/pub/taxonomy/taxdump.tar.gz"

# Extract taxonomy
mkdir ncbi20150518
cd ncbi20150518/
tar xf ../taxdump.tar.gz
cd ..
```

### Format data for OBITools

```
# Activate OBITools
cd /home/user/programs
./obitools

# Return to working directory
cd /home/user/mitochondria2

# Format taxonomy
obitaxonomy -t ncbi20150518 -d ncbi20150518

# Merge genomes
obiconvert mitochondria/* > mito.all.fasta
head -5 mito.all.fasta

# Select vertebrates
ecofind -d ncbi20150518 '^vertebrata$'
obigrep -d ncbi20150518 -r 7742 mito.all.fasta > mito.vert.fasta

# Convert to ecoPCR database
obiconvert -d ncbi20150518 --ecopcrdb-output=mito.vert mito.vert.fasta
```







