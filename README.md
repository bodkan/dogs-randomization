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
