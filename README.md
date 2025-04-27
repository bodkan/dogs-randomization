This repository contains code for the bootstrapping procedure supporting the finding of the
sharing of ROH-depleted "deserts" (genomic windows) between ancient and modern dogs, as
presented in [Bougiouri _et al_., 2024](https://www.biorxiv.org/content/10.1101/2024.03.15.585179v3).

## Reproducing the results

1. Clone the repository with `git clone https://github.com/bodkan/dogs-bootstrapping` (or just
download the zip file through the green "Code" button above).
2. Run `cd` into the directory with the cloned project and start R in it.
3. R will automatically download the [_renv_](https://rstudio.github.io/renv/) package for creating
a reproducible R environment.
4. Call `renv::restore()` which should download the package dependencies used.

If the _renv_ setup fails for whatever reason, you can also do this:

```
install.packages(c("data.table", "ggplot2", "dplyr", "tidyr", "tibble", "readr", "BiocManager"))
BiocManager::install("GenomicRanges")
```

to get everything manually without renv, or simply install install individual packages
from [this](bootstrap.R) script. All are standard R packages, so you shouldn't run into
any issues (although this repository contains a fixed, reproducible _renv_ R environment, exact
versions of R packages pretty much won't matter at all).

The log summary of the entire pipeline is in [`bootstrap.txt`](bootstrap.txt). It
was produced by the following command (which also generates all files and figures described below):

```
Rscript bootstrap.R > bootstrap.txt
```

The results of this script are:

- the full log information in [`bootstrap.txt`](bootstrap.txt),
- some intermediate files which can be used to reproduce the results without having to run
  all data processing steps and bootstrap procedure (`*.qs` and `*.rds` files) as well as
  the coordinates of deserts in [`deserts.tsv`](deserts.tsv) (reading the [script](bootstrap.R)
  will make it obvious at which point is which intermediate data set used),
- [`ecdf.pdf`](ecdf.pdf) and [`bootstrap.pdf`](bootstrap.pdf) which show the
[ECDF](https://en.wikipedia.org/wiki/Empirical_distribution_function) of the bootstrap desert
counts and their histogram with a p-value of the observed count, respectively.

## Description of the procedure

The entire workflow is implemented in a single R script:
[`bootstrap.R`](bootstrap.R).

Briefly, the script first finds
["ROH desert" windows](https://github.com/bodkan/dogs-bootstrapping/blob/main/bootstrap.txt#L18-L183)
in a set of ancient and modern dogs, using a cutoff "ROH frequency" for each site in a window of 5%
(determining which window is a desert or not) in each of those two groups of dogs.

It then reshuffles genomic windows in each individual (each window carrying SNPs
with an assigned `TRUE` or `FALSE` state, depending on whether or not it overlaps
an ROH in that individual) and performs the same desert-detection procedure on
each reshuffled data set to establish a bootstrapping distribution of expected
numbers of shared deserts. The p-value of the observed number of deserts is then
computed from the [ECDF](https://en.wikipedia.org/wiki/Empirical_distribution_function)
of the bootstrap values.

### Bit-encoding of ROH states at each genomic locus

Note that storing the entire data set with `TRUE`/`FALSE` ROH states for each
of the millions of sites in each of the 552 individuals is incredibly memory intensive,
and bootstrapping the tabular data in memory (to avoid repeated, unnecessary disk
access) would involve a prohibitively large amounts of memory copying when speeding
up the process via parallelization. In fact, a naive in-memory parallelized bootstrapping
approach was crashing a 750 GB RAM 100 CPU computing server.

To make the bootstrapping procedure fast and memory efficient, we store the `TRUE`/`FALSE`
ROH states at every site in each individual not as standard logical R vectors
(this would imply storing each logical state as a 32-bit integer, requiring 4 bytes of memory
for each site), but as literal bits, where `TRUE` is represented by bit 1, and `FALSE`
is represented by bit 0 in a binary form. In this more efficient representation, we can store
the ROH states of 32 genomic sites in the space of one site in the traditional integer-based
representation of logical values. As a result, a parallelized bootstrapping procedure
which otherwise crashed a huge computing server can now fit on a standard MacBook.

## Primary data

Here is a brief description of all primary [`data`](data/) used in the pipeline. For more details on
how this data was generated, please refer to
[paper](https://www.biorxiv.org/content/10.1101/2024.03.15.585179v3), the code base supporting
the study linked in the paper, or contact the corresponding authors.

- coordinates of ROH regions in modern dogs:
  - `data/ref-panel_allchrom_sample-snp_filltags_filter_MAF_0.01_all_sites_hom_win_het_1_dogs.hom`
  - used [here](https://github.com/bodkan/dogs-bootstrapping/blob/6f2069394334468b6ff060b9e04337a1937060f3/bootstrap.R#L24)

- coordinates of ROH regions in ancient dogs:
  - `data/merged_phased_annotated.allchrom_MAF_0.01_recalibrated_INFO_0.8_all_sites_hom_win_het_1_dogs.hom`
  - used [here](https://github.com/bodkan/dogs-bootstrapping/blob/6f2069394334468b6ff060b9e04337a1937060f3/bootstrap.R#L29)
  
- coordinates of "modern dog sites":
  - `data/ref-panel_allchrom_sample-snp_filltags_filter_MAF_0.01_all_sites_hom_win_het_1_dogs.hom.summary.gz`
  - used [here](https://github.com/bodkan/dogs-bootstrapping/blob/6f2069394334468b6ff060b9e04337a1937060f3/bootstrap.R#L52)
  
- coordinates of "ancient dog imputed sites" (a perfect subset of the modern dog sites):
  -`data/merged_phased_annotated.allchrom_MAF_0.01_recalibrated_INFO_0.8_all_sites_hom_win_het_1_dogs.hom.summary.gz`
  - used [here](https://github.com/bodkan/dogs-bootstrapping/blob/6f2069394334468b6ff060b9e04337a1937060f3/bootstrap.R#L48)
  
- coordinates of genomic windows to be determined as "ROH deserts" or not:
  - `data/dogs_allchrom_windows_cov_500kb.txt`
  - used [here](https://github.com/bodkan/dogs-bootstrapping/blob/6f2069394334468b6ff060b9e04337a1937060f3/bootstrap.R#L68)

