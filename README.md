1. Clone the repository with `git clone https://github.com/bodkan/dogs-bootstrapping`
2. Run `cd` into the directory with the cloned project and start R in it.
3. R will automatically download renv.
4. Call `renv::restore()` which should download the package dependencies used.

If renv setup fails for whatever reason, just do this:

```
install.packages(c("data.table", "ggplot2", "dplyr", "tidyr", "tibble", "readr", "BiocManager"))
BiocManager::install("GenomicRanges")
```

to get everything manually without renv, or just install individual packages you're missing.
You can see they're all fairly standard R data science packages so there hopefully won't be any issues.
