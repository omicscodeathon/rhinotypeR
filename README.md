# rhinotypeR: A package for rhinovirus genotyping

## Table of Contents
1. [Background](#Background)
2. [Test-Data](#Test-Data)
3. [Package](#Package)
4. [Citation](#Installation)
5. [Contributors](#Contributors)

## Background
Rhinoviruses (RV), common respiratory pathogens, are positive-sense, single-stranded RNA viruses characterized by a high antigenic diversity and mutation rate. With their genome approximately 7.2 kb in length, RVs exhibit mutation rates between 10^-3 and 10^-5 mutations per nucleotide per replication event. These viruses are classified into 169 types across three species: RV-A, RV-B, and RV-C. Genotype assignment, a critical aspect of RV research, is based on pairwise genetic distances and phylogenetic clustering with prototype strains, a process currently executed manually and laboriously. 

## Test-Data

The project utilizes VP4/2 sequences available in the public domain from
GenBank and reference prototype strains from www.picornaviridae.com

#### Data Format

The input datasets (target, reference and prototype) are fasta files.
Here’s an example of a FASTA file: 
![fasta file](https://github.com/omicscodeathon/rhinotyper/blob/main/man/figures/example_fasta_file.png)

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
devtools::install_github("omicscodeathon/rhinotyper")
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

| Function                  | Role                                                              | Input      | Output |
|---------------------------|-------------------------------------------------------------------|------------|--------|
| `assignTypes()`           | Assigns genotypes to query sequence                               | fasta file | CSV file |
|`PlotFrequency()`          | Plot the assigned genotypes                                       | output from assignTypes  | Barplot
| `getPrototypeSeqs()`      | Creates an output directory for the prototype files               | fasta file |        |
| `allPrototypeDistances()` | Generates pairwise distance between query and prototype sequences | fasta file |        |
| `PlotPrototypeDistances()` | Plots prototype distances | distance matrix from prototype distance function | Heatmap
| `overallMeanDistance()`   | Calculates overall genetic distance of query sequences            | fasta file |        |


### Running the functions using our test data

##### Function 1: Assign genotypes to query sequence

``` r
assignTypes(pathToRef, pathToQuery, model = "p-distance", threshold = 0.105)
#>              query    assigned_type distance
#> 1 AY040238.1_RVB92 AY040238.1_RVB92        0
```

##### Function 2: Plot frequency
```{r }
PlotFrequency(genotypeassigned)
```

##### Function 3: Create an output directory for the prototype files

``` r
#Create an output directory in your current directory before running this function
getPrototypeSeqs(destinationFolder = "./output")
#> [1] "The reference sequences have been downloaded to ./output"
```

##### Function 4: Generate pairwise distance between query and prototype sequences

``` r
allPrototypeDistances(RVBPrototype, pathToQuery, "p-distance")
#>                  HQ123444.1_RVB100genome JF781500.1_RVB101genome
#> AY040238.1_RVB92               0.1690476               0.1761905
#>                  JX074053.1_RVB102genome JN614996.1_RVB103genome
#> AY040238.1_RVB92               0.1714286               0.1785714
#>                  FJ445137.1_RVB104genome L05355.1_RVB14genome AF343645.1_RVB17
#> AY040238.1_RVB92               0.1952381            0.1738095        0.2071429
#>                  AF343653.1_RVB26 AF343654.1_RVB27 AY016403.1_RVB3
#> AY040238.1_RVB92        0.2261905        0.2238095       0.1928571
#>                  AY040241.1_RVB35 AY016401.1_RVB37 AF343655.1_RVB4
#> AY040238.1_RVB92        0.1714286        0.1833333       0.2333333
#>                  AY016404.1_RVB42 AY016400.1_RVB48 AF343651.1_RVB5
#> AY040238.1_RVB92        0.2095238        0.2285714       0.2119048
#>                  AY016398.1_RVB52 AY016402.1_RVB6 AY016399.1_RVB69
#> AY040238.1_RVB92        0.2047619       0.1904762        0.2119048
#>                  AF343646.1_RVB70 AF343650.1_RVB72 AF343649.1_RVB79
#> AY040238.1_RVB92        0.2214286              0.2        0.1785714
#>                  AF343647.1_RVB83 AY040240.1_RVB84 AF343648.1_RVB86
#> AY040238.1_RVB92        0.1357143        0.2071429        0.1642857
#>                  AY040237.1_RVB91 AY040238.1_RVB92 AY040239.1_RVB93
#> AY040238.1_RVB92        0.2214286                0        0.2214286
#>                  AY040242.1_RVB97 AF343652.1_RVB99
#> AY040238.1_RVB92        0.1952381        0.2214286
```

##### Function 5: Plot Prototype distances
```{r }
#PlotPrototypeDistances(pairwisedist)
```

##### Function 6: Calculate overall genetic distance of query sequences

``` r
overallMeanDistance(RVAPrototype)
#> [1] 0.2203165
```

### Running the functions on your own dataset

##### Function 1: Assign genotypes to query sequence

The model could be either p-distance, JC, Kimura2p, Tamura3p based on
your preference. The threshhold can also be changed based on your needs.

``` r
#assignTypes("path to reference sequence e.g RVBPrototype", "path to query sequence", model = "Tamura3p", threshold = 0.105)
```

##### Function 2: Plot frequency
The output of the assign genotype (function 1) is the input.
```{r }
#genotypeassigned <- assignTypes(pathToRef, pathToQuery, model = "p-distance", threshold = 0.105)
#PlotFrequency(genotypeassigned)
```

##### Function 3: Create an output directory for the prototype files

Ensure you create an output directory in your current directory before
running this function, otherwise it will throw an error saying “output”
folder is not found.

``` r
getPrototypeSeqs(destinationFolder = "./output")
#> [1] "The reference sequences have been downloaded to ./output"
```

##### Function 4: Generate pairwise distance between query and prototype sequences

The model could be either p-distance, JC, Kimura2p, Tamura3p based on
your preference.

``` r
#allPrototypeDistances("path to reference sequence e.g RVBPrototype", "path to query sequence", "Tamura3p")
```

##### Function 5: Plot Prototype distances
The output of the generate pairwise distance between query and prototype sequences (function 4) is the input data.
```{r }
#PlotPrototypeDistances(pairwisedist)
```

##### Function 6: Calculate overall genetic distance of query sequences

``` r
#overallMeanDistance("path to reference sequence e.g RVAPrototype")
```


## Citation

## Contributors

   - Martha Luka

   - Winfred Gatua

   - Wafaa Rashed

   - Ruth Nanjala
