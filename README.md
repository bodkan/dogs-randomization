## Reproducing the results

1. Clone the repository with `git clone https://github.com/bodkan/dogs-bootstrapping`
2. Run `cd` into the directory with the cloned project and start R in it.
3. R will automatically download renv.
4. Call `renv::restore()` which should download the package dependencies used.

If the renv setup fails for whatever reason, you can also do this:

```
install.packages(c("data.table", "ggplot2", "dplyr", "tidyr", "tibble", "readr", "BiocManager"))
BiocManager::install("GenomicRanges")
```

to get everything manually without renv, or simply install install individual packages
from [this](bootstrap.R) script.

The log summary of the entire pipeline is in [`bootstrap.txt`](bootstrap.txt). It
was produced by the following command (which also generates all associated output
files and figures):

```
Rscript bootstrap.R > bootstrap.txt
```

The results of this script are:

- the full log information in [`bootstrap.txt`](bootstrap.txt),
- some intermediate files which can be used to reproduce the results without having to run
  all data processing steps and bootstrap procedure (`*.qs` and `*.rds` files) as well as
  the coordinates of deserts in [`deserts.tsv`](deserts.tsv),
- [`ecdf.pdf`](ecdf.pdf) and [`bootstrap.pdf`](bootstrap.pdf) which show the [ECDF](https://en.wikipedia.org/wiki/Empirical_distribution_function)
  of the bootstrap desert counts and their histogram with a p-value of the observed count,
  respectively.

## Description of the procedure

The entire workflow is implemented in a single R script:
[`bootstrap.R`](bootstrap.R).

Briefly, the R script first finds "ROH desert" windows in a set of ancient and
modern dogs, using a cutoff "ROH frequency" for each site in a window of 5% in
each of those two groups.

It then reshuffles genomic windows in each individual (each window carrying SNPs
with an assigned `TRUE` or `FALSE` state, depending on whether or not it overlaps
an ROH in that individual) and performs the same desert-detection procedure on
each reshuffled data set to establish a bootstrapping distribution of expected
numbers of shared deserts. The p-value of the observed number of deserts is then
computed from the [ECDF](https://en.wikipedia.org/wiki/Empirical_distribution_function)
of the bootstrap values.

### Bit-encoding of ROH states at each genomic locus

Note that storing the entire data set with `TRUE`/`FALSE` ROH states for each
of the millions of sites in each of the 552 individuals is quite memory intensive,
and bootstrapping the tabular data in memory (to avoid repeated, unnecessary disk
access) would involve a prohibitively large amounts of memory copying when speeding
up the process via parallelization.

To make the bootstrapping procedure fast and memory efficient, we store the `TRUE`/`FALSE`
ROH states at every site in each individual not as standard logical R vectors
(this would imply storing each logical state as a 32-bit integer, requiring 4 bytes of memory
for each site), but as literal bits, where `TRUE` is represented by bit 1, and `FALSE`
is represented by bit 0. In this more efficient representation, we can store the ROH
states of 32 genomic sites in the space of one site in the traditional integer-based
representation of logical values. As a result, a parallelized bootstrapping procedure
which otherwise crashed a 750 GB computing server can now fit on a standard MacBook.
