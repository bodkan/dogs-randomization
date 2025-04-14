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

###############################################################
# read all input data
###############################################################

# coordinates of windows
win_gr <- fread("cat data/dogs_allchrom_windows_cov_500kb.txt | cut -f1-3 | sort | uniq") %>%
  setNames(c("chrom", "start", "end")) %>%
  mutate(chrom = as.integer(gsub("chr", "", chrom))) %>%
  arrange(chrom, start, end) %>%
  makeGRangesFromDataFrame()
win_gr

# ROH in modern samples
modern_gr <- fread("data/ref-panel_allchrom_sample-snp_filltags_filter_MAF_0.01_all_sites_hom_win_het_1_dogs.hom")
modern_gr <- tibble(modern_gr) %>% select(sample = FID, chrom = CHR, start = POS1, end = POS2)
modern_gr$set <- "modern"
modern_gr

# ROH in ancient samples
ancient_gr <- fread("data/merged_phased_annotated.allchrom_MAF_0.01_recalibrated_INFO_0.8_all_sites_hom_win_het_1_dogs.hom")
ancient_gr <- tibble(ancient_gr) %>% select(sample = FID, chrom = CHR, start = POS1, end = POS2)
ancient_gr$set <- "ancient"
ancient_gr

# combine ROH sets into a single GRanges object
roh_gr <- rbind(modern_gr, ancient_gr)
roh_gr <- makeGRangesFromDataFrame(roh_gr, keep.extra.columns = TRUE)
roh_gr

# ancient sites
ancient_sites <- fread("data/merged_phased_annotated.allchrom_MAF_0.01_recalibrated_INFO_0.8_all_sites_hom_win_het_1_dogs.hom.summary.gz")
ancient_sites <- tibble(ancient_sites) %>% select(chrom = CHR, pos = BP) %>% mutate(ancient = TRUE)
# modern sites
modern_sites <- fread("data/ref-panel_allchrom_sample-snp_filltags_filter_MAF_0.01_all_sites_hom_win_het_1_dogs.hom.summary.gz")
modern_sites <- tibble(modern_sites) %>% select(chrom = CHR, pos = BP) %>% mutate(modern = TRUE)
# combine both sets for simplicity of the downstream code
merged_sites <- left_join(modern_sites, ancient_sites, by = c("chrom", "pos")) %>%
  mutate(
    ancient = if_else(is.na(ancient), FALSE, ancient)
  )
sites_gr <- makeGRangesFromDataFrame(merged_sites, keep.extra.columns = TRUE, start.field = "pos", end.field = "pos")

# extract sample names for later use
ancient_samples <- unique(roh_gr[roh_gr$set == "ancient"]$sample)
ancient_samples
modern_samples <- unique(roh_gr[roh_gr$set == "modern"]$sample)
modern_samples

all_samples <- c(ancient_samples, modern_samples)
all_samples

###############################################################
# detecting sites overlapping ROH in each individual
###############################################################

# assign each site to a window (adding a window index column)
assign_sites <- function(sites_gr, win_gr) {
  # find which window overla each SNP
  hits <- findOverlaps(sites_gr, win_gr)

  # remove SNPs not hitting a window
  sites_gr <- sites_gr[queryHits(hits)]

  # annotate each SNP with its respective window
  sites_gr$win_i <- subjectHits(hits)

  # return the modified SNP table
  return(sites_gr)
}

# Score each site in a given individual as TRUE or FALSE depending on
# whether or not it overlaps with ROH. This code is probably terribly
# inefficient but works as long as its run with a machine with a huge
# amount of RAM (I haven't benchmarked this, but ran it on a 750 Gb RAM
# number cruncher server).
sites_coverage <- function(samples, sites_gr, roh_gr) {
  hits_inds <- mclapply(seq_along(samples), function(i) {
    ind <- samples[i]

    # extract the flag for whether that individual is ancient or modern
    set <- roh_gr[roh_gr$sample == ind]$set[1]

    cat(sprintf("Processing %s (%s) [%s/%s]\n", ind, set, i, length(samples)))

    # take a respective subset of all sites (either ancient set or modern set)
    subset <- mcols(sites_gr)[[set]]

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
# TODO CHECK THAT THE MODERN AND ANCIENT COLUMNS CORRESPOND REALLY TO THE NA VALUES IN ANCIENT
# TODO ALSO IMPLEMENT THE PADDING
  return(hits_df)
}

# assign a window number to each site
sites_gr <- assign_sites(sites_gr, win_gr)

if (!file.exists("cov_df.qs")) {
  # detect TRUE or FALSE for each site in each individual depending on whether or
  # not a given site overlaps an ROH in that individual
  cov_df <- sites_coverage(all_samples, sites_gr, roh_gr)
  qsave(cov_df, "cov_df.qs", preset = "high")
} else {
  cov_df <- qread("cov_df.qs")
}

# 3. just a visual test that the # of ROHs in ancient vs modern match our expectation
# TODO: adapt to data.table in a single merged form
#roh_counts <- lapply(coverages, function(ind) {
#  data.frame(set = ind$set[1], ind = ind$sample[1], roh_coverage = mean(ind$coverage))
#}) %>% do.call(rbind, .)
#ggplot(roh_counts) +
#  geom_jitter(aes(set, roh_coverage, color = set))


###############################################################
# Bootstrap components (1, 2, 3)
###############################################################

# shuffle sites in a given individual (keeping the sites sitting on the window together),
# then return a corresponding shuffled vector of TRUE/FALSE/NA ROH status of each site
# in this individual
shuffle_ind <- function(ind, cov_df, win_list) {
  # shuffle the list of windows/rows to get new row indices to...
  shuffled_rows <- do.call(rbind, sample(win_list))
  # shuffle the sites in this individual sample
  shuffled_cov_df[shuffled_rows$row, ..ind]
  shuffled_cov_df
}

# shuffle windows (and, therefore, sites) in all given individuals
shuffle_samples <- function(samples, cov_df) {
  # get a list of rows of the sites table corresponding to each window
  win_list <- cov_df[, .(win_i, row = 1:.N)] %>% { split(., .$win_i) }

  shuffled_cov_df <- mclapply(
    samples, function(ind) {
      cat(sprintf("Shuffling sites/windows in %s\n", ind))
      shuffle_ind(ind, cov_df, win_list)
    }, mc.cores = detectCores()) %>%
    do.call(cbind, .) %>%
    cbind(cov_df[, .(chrom, pos, win_i)], .)

  shuffled_cov_df
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
# sanity checking the desert inference in the original data
###############################################################

mean_cov_df <- cov_df[, lapply(.SD, mean, na.rm = TRUE), by = win_i, .SDcols = all_samples]

deserts_ancient <- detect_deserts(mean_cov_df[, ..ancient_samples], cutoff = 0.05)
deserts_modern <- detect_deserts(mean_cov_df[, ..modern_samples], cutoff = 0.05)

sum(deserts_ancient)
# 2677
sum(deserts_modern)
# 212

# count which windows represent a shared A vs M desert
shared_deserts <- deserts_ancient & deserts_modern
sum(shared_deserts)
# 159

# For debugging purposes, add frequencies of ROH overlapping each
# desert (in ancient and modern individuals) to the original table of windows
original_win <- win_gr
original_win$cov_ancient <- rowMeans(df[, ancient_samples])
original_win$cov_modern <- rowMeans(df[, modern_samples])
original_win$desert <- shared_deserts
original_win

# just the deserts in ancient individuals
original_win[original_win$cov_ancient < 0.05]

# just the deserts in modern individuals
original_win[original_win$cov_modern < 0.05]

# deserts in both ancient and modern individuals
original_win[original_win$desert]









###############################################################
# bootstrapping begins here
###############################################################

# run a single replicate of the desert inference
run_replicate <- function(windows, rep_i) {
  # 1. reshuffle windows in each individual
  shuffled_win <- shuffle_windows(windows)

  # 2. convert list of GRanges into a single merged data frame
  shuffled_merged <- merge_windows(shuffled_win)

  # 3. detect which ancient vs modern deserts are shared (producing a TRUE/FALSE vector)
  deserts_ancient <- detect_deserts(shuffled_merged[, ancient_samples], cutoff = 0.05)
  deserts_modern <- detect_deserts(shuffled_merged[, modern_samples], cutoff = 0.05)

  # detect which reshuffled windows correspond to a shared A vs M desert
  deserts_shared <- deserts_ancient & deserts_modern

  # return the result in a tidy form
  df <- tibble(
    rep_i = rep_i,
    win_i = seq_along(deserts_shared),
    match = deserts_shared
  )

  return(df)
}

# run 100 bootstrap desert reshuffling simulations
reps <- mclapply(1:100, function(rep_i) run_replicate(windows, rep_i = rep_i),
                 mc.cores = detectCores()) %>%
  do.call(rbind, .) %>%
  mutate(match = as.integer(match))

# the bootstrap loop produces a data frame with three columns:
#   - rep_i: iteration number
#   - win_i: index of a window (1:4-thousand-and0something)
#   - match: for each window, whether that window is a shared desert (1) or not (0)
reps

# Katia's original count of the shared deserts
empirical_count <- 124

# counts observed in each bootstrapping iteration
boot_counts <- group_by(reps, rep_i) %>% summarise(count = sum(match))
boot_counts

# histogram of the bootstrap counts
ggplot(boot_counts) + geom_histogram(aes(count))


###############################################################
# A toy example of how to then use this in a significance test
###############################################################

# get some count values (let's say from bootstrapping)
counts <- rpois(1000, 10)
counts

# let's say this is the observed count
obs_count <- 20

# plot the 'bootstrap distribution' along with the 'observed value'
plot(density(counts))
abline(v = obs_count)

# compute the empirical CDF
e <- ecdf(counts)
plot(e)
abline(v = obs_count)

# what's the probability of observing a value as extreme (or more extreme)
# than the value we observed?
1 - e(obs_count)












snps <- makeGRangesFromDataFrame(
  data.frame(ind = "ind1", chrom = "chr1",
             start = c(1, 6, 10, 20, 100, 1000, 2000),
             end = c(1, 6, 10, 20, 100, 1000, 2000),
             win_i = c(1, 1, 2, 2, 2, 3, 3)),
  keep.extra.columns = TRUE)
snps

sample(snps$win_i)

split(snps, snps$win_i) %>% sample %>% unlist
