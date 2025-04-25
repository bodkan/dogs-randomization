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
from [this](01_detect_deserts.R) or [this](02_bootstrap_deserts.R) script.

## Description of the procedure

The entire workflow is implemented in a single R script:
[`bootstrap.R`](bootstrap.R).

Briedfly, the R script first finds "ROH desert" windows in a set of ancient and
modern dogs, using a cutoff "ROH frequency" for each site in a window of 5% in
each of those two groups.

It then reshuffles genomic windows in each individual (each window carrying SNPs
with an assigned `TRUE` or `FALSE` state, depending on whether or not it overlaps
an ROH in that individual) and performs the same desert-detection procedure on
each reshuffled data set to establish a bootstrapping distribution of expected
numbers of shared deserts.

### Bit-encoding of ROH states at each genomic locus

Note that storing the entire data set with `TRUE`/`FALSE` ROH states for each
of the millions of sites in each of the 552 individuals is quite memory intensive,
and bootstrapping the tabular data in memory (to avoid extremely expensive disk
access) would involve a prohibitively large amounts of memory copying when speeding
up the process via parallelization.

To make the bootstrapping procedure fast and memory efficient, the program stores
the `TRUE`/`FALSE` ROH states in each individual not as standard logical R vectors
(this would imply storing each state as a 32-bit integer, requiring 4 bytes of memory
for each stat), but as literal bits, where `TRUE` is represented by bit 1, and `FALSE`
is represented by bit 0. In this more efficient representation, we can store the ROH
states of 32 genomic sites in the space of one site in the traditional integer-based
representation of logical values. As a result, a parallelized bootstrapping procedure
which otherwise crashed a 750 GB computing server can now fit on a standard MacBook.
