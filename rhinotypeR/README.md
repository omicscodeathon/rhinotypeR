
# rhinotypeR

<!-- badges: start -->
<!-- badges: end -->

The goal of rhinotypeR is to perform genotyping on rhinovirus sequences.

- The package does the following:
  - Assigns genotypes to query sequences
  - Computes for pairwise distance among query sequences
  - Calculates pairwise distance between query and prototype sequences
  - Calculates overall genetic distance of query sequences

## Data Format

The input data is a reference fasta file.

## Installation

You can install the development version of rhinotypeR from
[GitHub](https://github.com) with:

``` r
devtools::install_github("omicscodeathon/rhinotyper/rhinotypeR")
```

## Example

How to run the functions using our test data:

``` r
library(Rhinotyper)
```

``` r
assignTypes(pathToRef, pathToQuery, model = "p-distance", threshold = 0.105)
#>              query    assigned_type distance
#> 1 AY040238.1_RVB92 AY040238.1_RVB92        0
```

``` r
#Create an output directory in your current directory before running this function
getPrototypeSeqs(destinationFolder = "./output")
#> [1] "The reference sequences have been downloaded to ./output"
```

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

``` r
overallMeanDistance(RVAPrototype)
#> [1] 0.2203165
```
