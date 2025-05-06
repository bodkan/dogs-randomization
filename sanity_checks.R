# install.packages("data.table")
library(data.table)

# read names of ancient and modern samples
ancient_samples <- fread("https://raw.githubusercontent.com/bodkan/dogs-bootstrapping/refs/heads/main/data/merged_phased_annotated.allchrom_MAF_0.01_recalibrated_INFO_0.8_all_sites_hom_win_het_1_dogs.hom")[, unique(FID)]
modern_samples <- fread("https://raw.githubusercontent.com/bodkan/dogs-bootstrapping/refs/heads/main/data/ref-panel_allchrom_sample-snp_filltags_filter_MAF_0.01_all_sites_hom_win_het_1_dogs.hom")[, unique(FID)]

# read ROH status of each window in each sample, i.e. the proportion of sites in each
# window (we have 4300 windows in total) covered by some ROH in that sample
#   - rows are windows
#   - columns are individuals
#   - each element of the matrix is the proportion of sites in a given window,
#     in a given individual, covered by an ROH in that individual
# (I think this is like Fernando's 'matrix' above, except we have a proportion of each window
# covered by an ROH, not a 0/1 value)
mean_win_df <- fread("https://raw.githubusercontent.com/bodkan/dogs-bootstrapping/refs/heads/main/mean_win_df.tsv")
# for instance, this is a subset of this data for the first five ancient samples
mean_win_df[, ..ancient_samples][, 1:5]

# for starters, let's focus on ancient dogs because it's much smaller (uncomment
# the second line to get numbers for modern dogs)
samples <- ancient_samples
# samples <- modern_samples

# get what Fernando calls 'propwins'
# -- for each sample, compute mean ROH frequency across all of its windows
(propwins <- colMeans(mean_win_df[, ..samples]))
# get what Fernando calls 'propinds'
# -- for each window, compute mean ROH frequency across all individuals
(propinds <- rowMeans(mean_win_df[, ..samples]))

# Fernando's plots, now using the empirical data

# Plot to match Katia's histogram
hist(propwins, main = "Histogram of individuals", xlab = "Proportion of individuals with X percent ROH", xlim = c(0, 1), breaks = 50)

hist(propinds, main = "Histogram of windows", xlab = "Proportion of windows with X percent individuals who have ROH there", xlim = c(0, 1), breaks = 50)
abline(v = 0.05, col = "red", lty = 2)

# Our definition of a desert (before the intersection)
sum(propinds < 0.05)
# This should give the following numbers:
# Number of deserts in ancient individuals: 2623
# Number of deserts in modern individuals: 164
