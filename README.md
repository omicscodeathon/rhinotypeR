# rhinotypeR: A package for rhinovirus genotyping

```

           /**       /**                       /**                                   /*** *** 
          | **      |**/                      | **                                  | **__  **
  /****** | *******  /** /*******   /******  /******   /**   /**  /******   /****** | **  \ **
 /**__  **| **__  **| **| **__  ** /**__  **|_  **_/  | **  | ** /**__  ** /**__  **| *******/
| **  \__/| **  \ **| **| **  \ **| **  \ **  | **    | **  | **| **  \ **| ********| **__  **
| **      | **  | **| **| **  | **| **  | **  | ** /**| **  | **| **  | **| **_____/| **  \ **
| **      | **  | **| **| **  | **|  ******/  |  *****/|  *******/| ********|  *******/| **  | 
|__/      |__/  |__/|__/|__/  |__/ \______/    \___/   \____  **| **____/  \_______/|__/  |__/
                                                       /**  | **| **                          
                                                      |  ******/| **                          
                                                       \______/ |__/                          

```

## Table of Contents
1. [Background](#Background)
2. [Test-Data](#Test-Data)
3. [Workflow](#Workflow)
4. [Package](#Package)
5. [Citation](#Citation)
6. [Contributors](#Contributors)

## Background
Rhinoviruses (RV), common respiratory pathogens, are positive-sense, single-stranded RNA viruses characterized by a high antigenic diversity and mutation rate. With their genome approximately 7.2 kb in length, RVs exhibit mutation rates between 10^-3 and 10^-5 mutations per nucleotide per replication event. These viruses are classified into 169 types across three species: RV-A, RV-B, and RV-C. Genotype assignment, a critical aspect of RV research, is based on pairwise genetic distances and phylogenetic clustering with prototype strains, a process currently executed manually and laboriously. 

## Test-Data

The project utilizes VP4/2 sequences available in the public domain from
GenBank and reference prototype strains from www.picornaviridae.com
The input datasets (target, reference and prototype) are fasta files.
Here’s an example of a FASTA file: 
![fasta file](https://github.com/omicscodeathon/rhinotyper/blob/main/man/figures/example_fasta_file.png)

## Workflow
![workflow](https://github.com/omicscodeathon/rhinotyper/blob/main/man/figures/workflow.png)
## Package
Our project aims to develop an R package to automate RV genotype assignment, facilitating genomic scientists in efficiently genotyping RV infections.

Our methodology involves:
1. Parsing and preprocessing of VP4/2 sequence data.
2. Implementation of algorithms to calculate pairwise genetic distances.
3. Integration of methods for constructing Maximum Likelihood phylogenetic 
trees.

### Installation

You can install the development version of rhinotypeR from
[GitHub](https://github.com) with:

``` r
devtools::install_github("omicscodeathon/rhinotypeR")
```

##### Load Library

``` r
library("rhinotypeR")
```

### Functions

The package encompasses functions to compute genetic distances, perform phylogenetic clustering, and compare sequences against RV prototype strains. 
These functionalities are designed to be user-friendly and adaptable to various research needs.

- The package (summarized in Table 1) does the following:
  - Assigns genotypes to query sequences
  - Computes for pairwise distance among query sequences
  - Calculates pairwise distance between query and prototype sequences
  - Calculates overall genetic distance of query sequences
 
#### Table 1. A summary of the functions

| Function        | Role                                   | Input                    | Output                     
|-----------------|----------------------------------------|--------------------------|-----------------------------
| `getPrototypeSeqs()`| Creates an output directory for the prototype files | fasta file | output directory
| `readFasta`()` | Assigns genotypes to query sequence | fasta file | CSV file
| `SNPeek()` | Assigns genotypes to query sequence | fasta file | CSV file
| `assignTypes()` | Assigns genotypes to query sequence | fasta file | CSV file  
| `allPrototypeDistances()` | Generates pairwise distance between query and prototype sequences | fasta file |  CSV file
| `pairwiseDistances()` | Assigns genotypes to query sequence | fasta file | CSV file
| `overallMeanDistance()` | Calculates overall genetic distance of query sequences  |  fasta file |  CSV file
| `countSNPs`()` | Assigns genotypes to query sequence | fasta file | CSV file
| `PlotFrequency()` | Plot the assigned genotypes | output from assignTypes | Barplot
| `PlotPrototypeDistances()` | Plots prototype distances | distance matrix from prototype distance function | Heatmap
| `PlotTree()` | Plot the pairwise distance | output from pairwiseDistances | Barplot


### Running the functions

##### Function 1: getPrototypeSeqs 
- Create an output directory for the prototype files
- Ensure you create an output directory in your current directory before running this function, otherwise it will throw an error saying "output" folder is not found.

Example
```{r }
getPrototypeSeqs(destinationFolder = "./output")
```

Own data
```{r }
#getPrototypeSeqs(destinationFolder = "path to an output folder")
```

##### Function 2: readFasta 
- Read sequences from a FASTA file

Example
```{r }
readFasta(RVAPrototype)
```

Own data
```{r }
#readFasta("path to fastaFile")
```

##### Function 3: SNPeek 
- Visualizing Single Nucleotide Polymorphisms (SNPs) in the fasta file

Example
```{r }
SNPeek(readFasta(RVAPrototype))
```

Own data
```{r }
#SNPeek(readFasta("path to fasta file"))
```

##### Function 4: assignTypes - Assign genotypes to the query sequence
- The model could be either p-distance, JC, Kimura2p, Tamura3p based on your preference.
- The threshhold can also be changed based on your needs.

Example
```{r }
assignTypes(RVBPrototype, readFasta(target_fasta_file), "p-distance", 0.105)
```

Own data
```{r }
#assignTypes("path to reference sequence e.g RVBPrototype", readFasta("path to query sequence"), model, threshold = 0.105)
```

##### Function 4: allPrototypeDistances 
- Calculates pairwise distance between query and prototype sequences
- The model could be either p-distance, JC, Kimura2p, Tamura3p based on your preference.

Example
```{r }
allPrototypeDistances(RVBPrototype, readFasta(target_fasta_file), "p-distance")
```

Own data
```{r }
#allPrototypeDistances("path to reference sequence e.g RVBPrototype", readFasta("path to query sequence"), "Tamura3p")
```

##### Function 5: pairwiseDistances 
- Calculate pairwise distance among input sequences

Example
```{r }
pairwiseDistances(RVBPrototype, "p-distance")
```

Own data
```{r }
#pairwiseDistances("path to prototype file", model = "p-distance")
```

##### Function 6: overallMeanDistance 
- Calculates the overall genetic distance of query sequences

Example
```{r }
overallMeanDistance(readFasta(RVAPrototype), model="p-distance")
```

Own data
```{r }
#overallMeanDistance(readFasta("path to reference sequence e.g RVAPrototype"),  model="p-distance")
```

##### Function 7: countSNPs 
- Count SNPs among input sequences

Example
```{r }
countSNPs(inputSequencesPath = RVBPrototypeAligned)
```

Own data
```{r }
#countSNPs(inputSequencesPath = "path to prototype file") 
```
##### Function 8: Plot frequency
- Create a barplot of genotype frequencies

Example
```{r }
plotFrequency(target_fasta_2, refSeq, "Tamura3p")
```

## Citation

## Contributors

   - Ruth Nanjala

   - Martha Luka

   - Winfred Gatua

   - Wafaa Rashed

 
