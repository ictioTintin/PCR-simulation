# In silico design and validation of fish-specific primers

>Versión en español aquí:\
>[PCR-simulation](https://github.com/ictioTintin/PCR-simulation)

_Martin Holguin Osorio_\
_January 2020_

> These experiments are based on the work done by Eric Coissac in 2015:\
> [Designing a new DNA metabarcode for fish](https://metabarcoding.org/IMG/html/primerdesign.html)

## Dependencies

To run and analyze the downloaded data, the following software is required:\
* A computer with a UNIX operating system, or a server or virtual machine capable of running command line code.  
* The computer or machine must have the following programs installed:  
  * Software package [OBITools](https://pythonhosted.org/OBITools/welcome.html#installing-the-obitools).  
  * Program [ecoPCR](https://git.metabarcoding.org/obitools/ecopcr/-/wikis/home)  
  * Program [ecoPrimers](https://git.metabarcoding.org/obitools/ecocebadores/-/wikis/home)  
  * R [r-project](https://www.r-project.org/) with the packages ROBITools, ROBITaxonomy, and ROBIBarcodes.  

Once everything is ready, we begin:

### Downloading the data

```bash
# Create a folder to store all data from my home directory
mkdir mitocondria2

# Enter the folder
cd mitochondria2

# Download complete mitochondrial DNA sequences from NCBI
wget "ftp://ftp.ncbi.nlm.nih.gov/genomes/refseq/mitochondrion/mitochondrion.1.genomic.gbff.gz"

wget "ftp://ftp.ncbi.nlm.nih.gov/genomes/refseq/mitochondrion/mitochondrion.2.genomic.gbff.gz"

# Download taxonomy associated with these sequences
wget 'ftp://ftp.ncbi.nlm.nih.gov://pub/taxonomy/taxdump.tar.gz'

# Extract taxonomy into a folder
mkdir ncbi20150518
cd ncbi20150518/
tar xf ../taxdump.tar.gz
cd ..
```

### Formatting the data for OBITools

```bash
# Activate obitools from where it is installed
cd /home/user/programs
./obitools

# Return to the working directory
cd /home/user/mitochondria2

# Format the taxonomy data for obitools
obitaxonomy -t ncbi20150518 -d ncbi20150518

# Merge all genomes into a single fasta file and verify
obiconvert mitochondria/* > mito.all.fasta
head -5 mito.all.fasta

# Search for vertebrates in the taxonomy files based on their taxid
ecofind -d ncbi20150518 '^vertebrata$'

# Check for existing errors "like a genus named vertebrata"
ecofind -d ncbi20150518 -p 1261581

# Fix errors "re-annotate and select the genomes"
obiannotate -d ncbi20150518             --with-taxon-at-rank=species             mito.all.fasta | obiannotate -S 'ori_taxid=taxid' | obiannotate -S 'taxid=species' | obiuniq -c taxid | obiselect -c taxid -n 1 -f count -M > mito.one.fasta

# Verify number of sequences after filtering
obicount mito.all.fasta

# Select vertebrate genomes
obigrep -d ncbi20150518 -r 7742 mito.one.fasta > mito.vert.fasta

# Format fasta file into an ecoPCR database
obiconvert -d ncbi20150518 --ecopcrdb-output=mito.vert mito.vert.fasta

# Search for teleosts based on their taxid
ecofind -d ncbi20150518 '^Teleostei$'
```

### Estimating the best primer pairs for the group of interest (Teleostei)

```bash
# Estimate the best primer sets with ecoPrimers and select the best pair from the generated file

# From the src folder of the program "in my case /home/user/programs/ecopcr"
./ecocebadores -d /home/user/mitochondria2/mito.vert  -e 3 -3 2  -l 30 -L 150 -r 32443 -c > Teleostei.ecocebadores
```

### Testing the new primer pair

```bash
# Simulate PCR with the two best primers provided in the file Teleostei.ecocebadores 

# From the src folder of the program "in my case /home/user/programs/ecopcr/src"
./ecoPCR -d /home/user/mitochondria2/mito.vert -e 5 -l 30 -L 300 -c ACACCGCCCGTCACTCTC ACCTTCCGGTACACTTAC > Teleostei.ecoPCR
```

### Start an R session

```r
# Open R and load packages
R
library(ROBITools)
library(ROBITaxonomy)
library(ROBIBarcodes)

# Load the results from ecoPCR along with taxonomy
fish = read.ecopcr.result('Teleostei.ecoPCR')
taxo = read.taxonomy('ncbi20150518')
# Verify
head(fish,n = 2)

# Identify which sequences are fish based on taxonomy data
teleo.taxid = ecofind(taxo,'^Teleostei$')
teleo.taxid

# Identify which sequences are fish among the ecoPCR results
is_a_fish=is.subcladeof(taxo,fish$taxid,teleo.taxid)
table(is_a_fish)

# Test primer binding site conservation
Fish.forward = ecopcr.forward.shanon(ecopcr = fish,
                                     group = is_a_fish)
Fish.reverse = ecopcr.reverse.shanon(ecopcr = fish,
                                     group = is_a_fish)
									 
# Plot results
pdf("primers_Teleostei.pdf")
par(mfcol=c(2,2))
dnalogoplot(Fish.forward$'TRUE',
            primer = "ACACCGCCCGTCACTCTC",
            main='Forward Fish')
dnalogoplot(Fish.forward$'FALSE',
            primer = "ACACCGCCCGTCACTCTC",
            main='Forward not Fish')

dnalogoplot(Fish.reverse$'TRUE',
            primer = "ACCTTCCGGTACACTTAC",
            main='Reverse Fish')
dnalogoplot(Fish.reverse$'FALSE',
            primer = "ACCTTCCGGTACACTTAC",
            main='Reverse not Fish')
dev.off()			

# Plot mismatches in the primers
pdf("mismatches.primers_Teleostei.pdf")
par(mfcol=c(1,1))
mismatchplot(fish,group = is_a_fish,
             legend=c('2722 non-fish vertebrates','2617 fish'))+title(xlab="Number of mismatches in the forward primer", 
                                   ylab="Number of mismatches in the reverse primer",
                                   main = 'Distribution of mismatches in the primer pair')
dev.off()			 
			 
# Check the taxonomic resolution of the primers	
only.fish=fish[is_a_fish,]

res = resolution(taxo,only.fish)
resolution = with(only.fish,
                  unique(data.frame(species_name,taxid,rank=res))
                 )
t(t(sort(table(resolution$rank)/length(resolution$rank),decreasing = TRUE)))


## Same analysis but for families

# Check amplification in Cyprinidae
# Create a group containing only Cyprinidae (instead of teleosts) based on its ID
is_a_cyprinidae = is.subcladeof(taxo,fish$taxid,7953)
only.cyprinidae=fish[is_a_cyprinidae,]

# Check primer resolution within Cyprinidae (instead of Teleostei)
res = resolution(taxo,only.cyprinidae)
resolution = with(only.cyprinidae,
                  unique(data.frame(species_name,taxid,rank=res))
                 )
t(t(sort(table(resolution$rank)/length(resolution$rank),decreasing = TRUE)))

# Check primer resolution for non-Cyprinidae fish
no.cyprinidae=fish[is_a_fish & !is_a_cyprinidae,]
res = resolution(taxo,no.cyprinidae)
resolution = with(no.cyprinidae,
                  unique(data.frame(species_name,taxid,rank=res))
                 )
t(t(sort(table(resolution$rank)/length(resolution$rank),decreasing = TRUE)))
```




