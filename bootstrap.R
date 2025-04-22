suppressPackageStartupMessages({
library(data.table)
library(R.utils)
library(qs)
library(ggplot2)
library(dplyr)
library(tidyr)
library(tibble)
library(readr)
library(parallel)
library(GenomicRanges)
library(testthat)
library(pryr)
})

###############################################################
# read all input data
###############################################################

cat("Reading input data... ")

# coordinates of windows
windows_gr <- fread("cat data/dogs_allchrom_windows_cov_500kb.txt | cut -f1-3 | sort | uniq") %>%
  setNames(c("chrom", "start", "end")) %>%
  mutate(chrom = as.integer(gsub("chr", "", chrom))) %>%
  arrange(chrom, start, end) %>%
  makeGRangesFromDataFrame()

# ROH in modern samples
modern_gr <- fread("data/ref-panel_allchrom_sample-snp_filltags_filter_MAF_0.01_all_sites_hom_win_het_1_dogs.hom")
modern_gr <- tibble(modern_gr) %>% select(sample = FID, chrom = CHR, start = POS1, end = POS2)
modern_gr$set <- "modern"

# ROH in ancient samples
ancient_gr <- fread("data/merged_phased_annotated.allchrom_MAF_0.01_recalibrated_INFO_0.8_all_sites_hom_win_het_1_dogs.hom")
ancient_gr <- tibble(ancient_gr) %>% select(sample = FID, chrom = CHR, start = POS1, end = POS2)
ancient_gr$set <- "ancient"

# combine ROH sets into a single GRanges object
roh_gr <- rbind(modern_gr, ancient_gr)
roh_gr <- makeGRangesFromDataFrame(roh_gr, keep.extra.columns = TRUE)

# extract sample names for downstream use
ancient_samples <- unique(roh_gr[roh_gr$set == "ancient"]$sample)
ancient_samples
modern_samples <- unique(roh_gr[roh_gr$set == "modern"]$sample)

all_samples <- c(ancient_samples, modern_samples)
qsave(all_samples, "all_samples.qs")

# ancient sites
ancient_sites <- fread("data/merged_phased_annotated.allchrom_MAF_0.01_recalibrated_INFO_0.8_all_sites_hom_win_het_1_dogs.hom.summary.gz")
ancient_sites <- tibble(ancient_sites) %>% select(chrom = CHR, pos = BP) %>% mutate(ancient = TRUE)

# modern sites
modern_sites <- fread("data/ref-panel_allchrom_sample-snp_filltags_filter_MAF_0.01_all_sites_hom_win_het_1_dogs.hom.summary.gz")
modern_sites <- tibble(modern_sites) %>% select(chrom = CHR, pos = BP) %>% mutate(modern = TRUE)

# apparently ancient sites are a perfect superset of modern sites (due to
# stricter filtering on some original set of sites)
expect_equal(ancient_sites, inner_join(ancient_sites, modern_sites, by = c("chrom", "pos"))[c("chrom", "pos", "ancient")])

# combine both sets to make the downstream code simpler (switching between
# ancient and modern column subsets as needed for each individual sample)
merged_sites <- left_join(modern_sites, ancient_sites, by = c("chrom", "pos")) %>%
  mutate(
    ancient = if_else(is.na(ancient), FALSE, ancient)
  )
merged_sites_gr <- makeGRangesFromDataFrame(merged_sites, keep.extra.columns = TRUE, start.field = "pos", end.field = "pos")

cat("done.\n")

###############################################################
# detecting sites overlapping ROH in each individual
###############################################################

# assign each site to a window (adding a window index column)
assign_sites <- function(sites_gr, windows_gr) {
  # find which window overla each SNP
  hits <- findOverlaps(sites_gr, windows_gr)

  # remove SNPs not hitting a window
  sites_gr <- sites_gr[queryHits(hits)]

  # annotate each SNP with its respective window
  sites_gr$win_i <- subjectHits(hits)

  # return the modified SNP table
  return(sites_gr)
}

# add 'dummy' NA sites to each window so that they have all the same length
# and can be reshuffled freely during bootstrap within each individual
add_padding <- function(sites_gr) {
  # split sites into window chunks
  snp_windows <- split(sites_gr, sites_gr$win_i)
  # what's the longest window in terms of the number of sites?
  max_length <- snp_windows %>% sapply(length) %>% max

  # loop over all windows with consecutive sites and add the required
  # padding of the required number of dummy sites as needed
  snp_windows_padded <- lapply(snp_windows, function(win_gr) {
    # how many dummy sites should we add to pad this window?
    padding_length <- max_length - length(win_gr)

    # add dummy sites if needed
    if (padding_length > 0) {
      padding_gr <- GRanges(seqnames = seqnames(win_gr)[1],
                            ranges = IRanges(start = rep(-1, padding_length), width = 1),
                            modern = NA, ancient = NA, win_i = win_gr$win_i[1])
    } else {
      padding_gr <- NULL
    }
    win_gr <- do.call(c, list(win_gr, padding_gr))
  }) %>% GRangesList(compress = FALSE)

  # turn the list of GRanges object into a single GRanges table again
  unlist(snp_windows_padded)
}

# Score each site in a given individual as TRUE or FALSE depending on
# whether or not it overlaps with ROH. This code needs a lot of RAM
# (I ran it on a 750 Gb RAM number cruncher machine).
sites_coverage <- function(samples, sites_gr, roh_gr) {
  hits_inds <- mclapply(seq_along(samples), function(i) {
    ind <- samples[i]

    # extract the flag for whether that individual is ancient or modern
    set <- roh_gr[roh_gr$sample == ind]$set[1]

    cat(sprintf("Processing %s (%s) [%s/%s]\n", ind, set, i, length(samples)))

    # take a respective subset of all sites (either ancient set or modern set)
    subset <- mcols(sites_gr)[[set]]
    # and only those sites which are not 'dummy sites' used for padding windows
    # to the same length
    subset <- subset & !is.na(subset)

    # detect which sites in this individual overlap with an ROH -- a site
    # in the overall data set will be either TRUE or FALSE (depending on whether
    # or not it overlaps or ROH) or NA (if its a site only available in the modern
    # set but an ancient individual is being currently processed)
    roh_hits <- rep(NA, length(subset))
    hits <- findOverlaps(sites_gr[subset], roh_gr[roh_gr$sample == ind])
    roh_hits[subset][queryHits(hits)] <- TRUE
    roh_hits[subset][-queryHits(hits)] <- FALSE

    roh_hits
  }, mc.cores = detectCores())
  names(hits_inds) <- samples

  sites_df <- as.data.table(sites_gr)[, .(chrom = seqnames, pos = start, win_i, modern, ancient)]
  hits_df <- cbind(sites_df, as.data.table(hits_inds))

  return(hits_df)
}

cat("Processing sites... ")

# assign a window number to each site
sites_gr <- assign_sites(merged_sites_gr, windows_gr)

# pad each window with dummy NA sites so that each window is of the same length
# (and we can reshuffle windows in each individual during bootstrapping)
sites_gr <- add_padding(sites_gr)

# sanity check -- all windows now must be of the same length
expect_true(length(unique(table(sites_gr$win_i))) == 1)

cat("done.\n")

cat("Assigning sites to ROH segments... ")

# assign TRUE or FALSE to each site in each individual depending on whether or
# not a given site overlaps an ROH in that individual
if (!file.exists("cov_df.qs")) {
  #cov_df <- sites_coverage(c(ancient_samples[1:3], modern_samples[1:3]), sites_gr, roh_gr)
  cov_df <- sites_coverage(all_samples, sites_gr, roh_gr)
  # cache the table to make things a little faster to iterate on
  qsave(cov_df, "cov_df.qs", preset = "high")
} else {
  cov_df <- qread("cov_df.qs")
}

# only ancient samples have NA values
expect_true(all(cov_df[!is.na(modern), lapply(.SD, function(x) any(is.na(x))), .SDcols = ancient_samples]))
# no modern sample has a NA value
expect_true(all(cov_df[!is.na(modern), lapply(.SD, function(x) all(!is.na(x))), .SDcols = modern_samples]))
# all ancient samples are NA at modern-only sites
expect_true(all(cov_df[!is.na(ancient) & (!ancient), lapply(.SD, function(x) all(is.na(x))), .SDcols = ancient_samples]))

cat("done.\n")

# 3. just a visual test that the # of ROHs in ancient vs modern match our expectation
# TODO: adapt to data.table in a single merged form
#roh_counts <- lapply(coverages, function(ind) {
#  data.frame(set = ind$set[1], ind = ind$sample[1], roh_coverage = mean(ind$coverage))
#}) %>% do.call(rbind, .)
#ggplot(roh_counts) +
#  geom_jitter(aes(set, roh_coverage, color = set))

#ancient_freq <- cov_df[(ancient), .(chrom, pos, win_i, count = rowSums(.SD, na.rm = TRUE)), .SDcols = ancient_samples]
#ancient_freq <- ancient_freq[pos != -1]
#ancient_sites <- fread("data/merged_phased_annotated.allchrom_MAF_0.01_recalibrated_INFO_0.8_all_sites_hom_win_het_1_dogs.hom.summary.gz")


###############################################################
# Detecting desert status of each window
###############################################################

# compute the proportion of sites covered by ROH for each window in each sample
windows_coverage <- function(samples, cov_df) {
  cov_df[, lapply(.SD, mean, na.rm = TRUE), by = win_i, .SDcols = samples]
}

# detect deserts in a given set of individuals (df here is a subset of
# the merged data frame, either a set of ancient individuals or a set of
# modern individuals, specified upon calling the function)
detect_deserts <- function(df, cutoff) {
  avg_coverages <- rowMeans(df)

  # flag each window as a desert or not
  cutoff <- 0.05
  deserts <- avg_coverages < cutoff

  return(deserts)
}

###############################################################
# Testing desert detection on the original (unshuffled) data
###############################################################

## compute ROH/SNP coverage in each window
#mean_win_df <- windows_coverage(all_samples, cov_df)
#
## assign desert status to each window (TRUE or FALSE)
#deserts_ancient <- detect_deserts(mean_win_df[, ..ancient_samples], cutoff = 0.05)
#deserts_modern <- detect_deserts(mean_win_df[, ..modern_samples], cutoff = 0.05)
#
## count the deserts
#sum(deserts_ancient)
## 2677
#sum(deserts_modern)
## 212
#shared_deserts <- deserts_ancient & deserts_modern
#sum(shared_deserts)
## 159
#
## for debugging purposes, add ROH frequencies for each window (in ancient and
## modern individuals) to the original table of windows for easier reference
#original_win <- windows_gr
## remove the window which is missing any SNPs
#original_win <- original_win[sort(unique(cov_df$win_i))]
#
#original_win$mean_ancient <- rowMeans(mean_win_df[, .SD, .SDcols = ancient_samples])
#original_win$mean_modern <- rowMeans(mean_win_df[, .SD, .SDcols = modern_samples])
#original_win$desert <- shared_deserts
#original_win
#
## just the deserts in ancient individuals
#original_win[original_win$mean_ancient < 0.05]
#
## just the deserts in modern individuals
#original_win[original_win$mean_modern < 0.05]
#
## deserts in both ancient and modern individuals
#original_win[original_win$desert]
#
## finally, save coordinates of shared ancient & modern deserts
#new_deserts <- original_win[original_win$desert]
#as.data.table(new_deserts)[, .(chrom = seqnames, start, end, mean_ancient, mean_modern, desert)] %>%
#  fwrite("new_deserts.tsv", sep = "\t", row.names = FALSE)

###############################################################
# Comparison of pre-review and post-review desert windows
###############################################################

# read coordinates of new deserts
new_deserts <-
  fread("new_deserts.tsv")[(desert)] %>%
  makeGRangesFromDataFrame(keep.extra.columns = TRUE)

# read coordinates of deserts from the pre-review manuscript
old_deserts <-
  fread("data/imputed_modern_window_bed_0.01_INFO_0.8_all_sites_hom_win_het_1_dogs_deserts.bed",
        col.names = c("chrom", "start", "end")) %>%
  makeGRangesFromDataFrame(keep.extra.columns = TRUE)

new_deserts
length(new_deserts)
# 159

old_deserts
length(old_deserts)
# 124

# detect deserts which were present both pre- and post-review
overlapping_deserts <- findOverlaps(new_deserts, old_deserts, type = "equal")
length(overlapping_deserts)
# 124 -- old desert windows are all within the new desert windows

###############################################################
# Bootstrapping
###############################################################

# shuffle sites in a given individual (keeping the sites sitting on the window together),
# then return a corresponding shuffled vector of TRUE/FALSE/NA ROH status of each site
# in this individual
shuffle_one <- function(ind_cov, win_list) {
  # shuffle the list of windows/rows to get new site indices in order to...
  shuffled_sites <- do.call(rbind, sample(win_list))
  # ... shuffle the sites in this individual sample
  shuffled_ind_cov <- ind_cov[shuffled_sites$row]
  shuffled_ind_cov
}

# shuffle windows (and, therefore, sites) in all given individuals
shuffle_samples <- function(samples, sites_df, cov_list) {
  # get a list of rows of the sites table corresponding to each window
  win_list <- sites_df[, .(win_i, row = 1:.N)] %>% { split(., .$win_i) }

  shuffled_cov_df <- mclapply(cov_list, shuffle_one, win_list = win_list, mc.cores = 20) %>% as.data.table

  shuffled_cov_df
}

suppressPackageStartupMessages({
    library(data.table)
    library(R.utils)
    library(qs)
    library(ggplot2)
    library(dplyr)
    library(tidyr)
    library(tibble)
    library(readr)
    library(parallel)
    library(GenomicRanges)
    library(testthat)
    library(pryr)
})


cat("Splitting table of ROH sites into a per-individual list... ")
if (!file.exists("cov_list.qs") || !file.exists("sites_df.qs")) {
  sites_df <- cov_df[, .(chrom, pos, win_i)]
  qsave(sites_df, "sites_df.qs", preset = "high")

  # split the combined ROH-sites-coverage data table into a list of sites in each
  # individual's window
  cov_list <- mclapply(all_samples, function(s) split(cov_df[[s]], cov_df$win_i), mc.cores = 50)
  names(cov_list) <- all_samples
  qsave(cov_list, "cov_list.qs", preset = "high")

  # free up memory before the bootstrap
  rm(cov_df)
  gc()
} else {
  all_samples <- qread("all_samples.qs")
  sites_df <- qread("sites_df.qs")
  cov_list <- qread("cov_list.qs")
}
cat("done.\n")

stop("asdf")

shuffled_cov_df <- shuffle_samples(all_samples, sites_df, cov_list)

# run a single replicate of the desert inference
run_replicate <- function(rep_i, sites_df, cov_list) {
  cat("Running replicate #", rep_i, "\n")

  # 1. reshuffle windows in each individual
  shuffled_cov_df <- shuffle_samples(all_samples, sites_df, cov_list)
  shuffled_cov_df <- cbind(sites_df, shuffled_cov_df)

  # 2. compute SNP-ROH coverage in each window
  mean_win_df <- windows_coverage(all_samples, shuffled_cov_df)

  # 3. detect which ancient vs modern deserts are shared (producing a TRUE/FALSE vector)
  deserts_ancient <- detect_deserts(mean_win_df[, ..ancient_samples], cutoff = 0.05)
  deserts_modern <- detect_deserts(mean_win_df[, ..modern_samples], cutoff = 0.05)

  # detect which reshuffled windows correspond to a shared A vs M desert
  deserts_shared <- deserts_ancient & deserts_modern

  # return the result in a tidy form
  df <- data.table(
    rep_i = rep_i,
    win_i = seq_along(deserts_shared),
    desert_ancient = as.integer(deserts_ancient),
    desert_modern = as.integer(deserts_modern),
    desert_shared = as.integer(deserts_shared)
  )

  return(df)
}

cat("Simulating bootstrap replicates... ")

if (!file.exists("bootstrap_reps.rds")) {
  tstart <- Sys.time()

  # run 100 desert reshuffling bootstrap replicates
  bootstrap_reps <- lapply(1:100, function(rep_i) run_replicate(rep_i, cov_df)) %>%
    do.call(rbind, .)

  tend <- Sys.time()

  cat("Total bootstrap time:", tend - tstart, "\n")

  saveRDS(bootstrap_reps, "bootstrap_reps.rds")

} else {
  bootstrap_reps <- readRDS("bootstrap_reps.rds")
}

cat("done.\n")

# the bootstrap loop produces a data frame with three columns:
#   - rep_i: replicate number
#   - win_i: index of a window
#   - desert_ancient: whether that window is an 'ancient desert' (1) or not (0)
#   - desert_modern: whether that window is an 'modern desert' (1) or not (0)
#   - desert_shared: whether that window is a both ancient and modern desert
bootstrap_reps

# the observed count of the shared deserts
observed_count <- length(new_deserts)

# counts observed in each bootstrapping iteration
bootstrap_counts <- bootstrap_reps[, .(desert_count = sum(desert_shared)), by = rep_i]
bootstrap_counts

###############################################################
# Putting a p-value on the result
###############################################################

# compute the empirical CDF
e <- ecdf(bootstrap_counts$desert_count)
plot(e, xlim = c(0, max(bootstrap_counts$desert_count)))
abline(v = observed_count, col = "red", lty = 2)

# what's the probability of observing a value as extreme (or more extreme)
# than the value we observed?
1 - e(observed_count)

# histogram of the bootstrap counts along with the observed value
ggplot(bootstrap_counts) +
  geom_histogram(aes(desert_count)) +
  geom_vline(xintercept = observed_count, linetype = "dashed", color = "red") +
  coord_cartesian(xlim = c(0, max(bootstrap_counts$desert_count))) +
  labs(x = "number of shared deserts", y = "replicate simulations") +
  theme_minimal()
