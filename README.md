# Rhinotyper: A package for rhinovirus ghenotyping

## Table of Contents
1. [Background](#Background)
2. [Test-Data ](#Test-Data )
3. [Package](#Package)
4. [Installation](#Installation)
5. [Functions](#Functions)
6. [License](#License)
7. [Contact](#Contact)
8. [Citation](#Citation)
9. [Contributors](#Contributors)


## Background

Rhinoviruses (RV), common respiratory pathogens, are positive-sense, single-stranded RNA viruses characterized by a high antigenic diversity and mutation rate. With their genome approximately 7.2 kb in length, RVs exhibit mutation rates between 10^-3 and 10^-5 mutations per nucleotide per replication event. These viruses are classified into 169 types across three species: RV-A, RV-B, and RV-C. Genotype assignment, a critical aspect of RV research, is based on pairwise genetic distances and phylogenetic clustering with prototype strains, a process currently executed manually and laboriously. 

## Test-Data 
The project will utilize VP4/2 sequences available in the public domain from GenBank and reference prototype strains from www.picornaviridae.com  

## Package
Our project aims to develop an R package to automate RV genotype assignment, facilitating genomic scientists in efficiently genotyping RV infections.

Our methodology involves:
1. Parsing and preprocessing of VP4/2 sequence data.
2. Implementation of algorithms to calculate pairwise genetic distances.
3. Integration of methods for constructing Maximum Likelihood phylogenetic 
trees.


### Package installation 

The package can be installed as follows:

    install.packages("devtools")

    library(devtools)

    devtools::install_github("omicscodeathon/rhinotyper/Package")

    library(rhinotyper)
    
## Functions
The package will encompass functions to compute genetic distances, perform phylogenetic clustering, and compare sequences against RV prototype strains. 
These functionalities will be designed to be user-friendly and adaptable to various research needs.

The package consists of 6 data analysis functions (processfastq(); counts(); expression(); callsnp(); callcnv(); and callindel()), 3 data visualization functions (vizexp(), vizsnp(), and vizcnv()), and a requirement() function is also provided to install all the dependencies. The list of functions and their roles are represented in Table 1.

| Function        | Role                                   | Input                    | Output                     |
|-----------------|----------------------------------------|--------------------------|-----------------------------|
| `requirement()` | Install required packages              | -                        | -                           |
| `processfastq()`| Preprocess fastq files                 | Fastq files              | BAM files                   |
| `counts()`      | Gene count analysis                    | BAM files                | CSV file                    |
| `expression()`  | Identify DEGs                          | BAM files                | CSV file                    |
| `callsnp()`     | SNP calling                            | BAM files                | VCF files                   |
| `callcnv()`     | CNV calling                            | BAM files                | CSV file                    |
| `callindel()`   | Indel calling                          | BAM files                | VCF files                   |
| `vizexp()`      | Analyze and visualize gene expression  | CSV file                 | Interactive interface       |
| `vizsnp()`      | Analyze and visualize SNP data          | VCF files                | Interactive interface       |
| `vizcnv()`      | Analyze and visualize CNV data          | CSV file                 | Interactive interface       |

#### Programming language: R 4.2.1

#### Requirements: 

## License  

## Citation

## Contributors

   - Martha Luka

   - Winfred Gatua

   - Wafaa

   - Ruth Nanjala
