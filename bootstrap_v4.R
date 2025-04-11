library(data.table)
library(ggplot2)
library(dplyr)
library(tidyr)
library(tibble)
library(readr)
library(parallel)
library(GenomicRanges)
library(windowscanr)

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
# detecting windows overlapping ROH in each individual
###############################################################

# assign each SNP to a window
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
# whether or not it overlaps with ROH
sites_coverage <- function(samples, sites_gr, roh_gr) {
  # make a copy of site coordinates only
  hits_gr <- granges(sites_gr)
  hits_gr$win_i <- sites_gr$win_i

  for (i in seq_along(samples)) {
    ind <- samples[i]

    # extract ROHs only in the given individual
    ind_roh_gr <- roh_gr[roh_gr$sample == ind]

    # extract the flag for whether that individual is ancient or modern
    set <- ind_roh_gr$set[1]

    cat(sprintf("Processing %s (%s) [%s/%s]\n", ind, set, i, length(samples)))

    # take a respective subset of all sites (either ancient set or modern set)
    subset <- mcols(sites_gr)[[set]]

    # detect which sites in this individual overlap with an ROH -- a site
    # in the overall data set will be either TRUE or FALSE (depending on whether
    # or not it overlaps or ROH) or NA (if its a site only available in the modern
    # set but an ancient individual is being currently processed)
    roh_hits <- rep(FALSE, length(sites_gr[subset]))
    hits <- findOverlaps(sites_gr[subset], ind_roh_gr)
    roh_hits[queryHits(hits)] <- TRUE
    mcols(hits_gr)[[ind]] <- rep(NA, length(hits_gr))
    mcols(hits_gr[subset])[[ind]] <- roh_hits
  }

  return(hits_gr)
}

# assign a window number to each site
sites_gr <- assign_sites(sites_gr, win_gr)

# detect TRUE or FALSE for each site in each individual depending on whether or
# not a given site overlaps an ROH in that individual
coverages <- sites_coverage(all_samples, sites_gr, roh_gr)
saveRDS(coverages, "coverages.rds")

# a couple of sanity checks:

# 1. the number of SNPs in each sample must be the same
sapply(coverages, length) %>% unique

# 2. the number of ancient and modern samples must be matching the original data
sapply(coverages, function(ind) unique(ind$set)) %>% table

# 3. just a visual test that the # of ROHs in ancient vs modern match our expectation
roh_counts <- lapply(coverages, function(ind) {
  data.frame(set = ind$set[1], ind = ind$sample[1], roh_coverage = mean(ind$coverage))
}) %>% do.call(rbind, .)
ggplot(roh_counts) +
  geom_jitter(aes(set, roh_coverage, color = set))


###############################################################
# Bootstrap components (1, 2, 3)
###############################################################

# 1. take a list of sites for each individual as a GRanges object and
# shuffle them in blocks in which they are sitting together in each window
shuffle_windows <- function(sites_gr) {
  shuffled_windows <- lapply(sites_gr, function(ind_sites_gr) {
    split(ind_sites_gr, ind_sites_gr$win_i) %>% sample %>% unlist
  })

  return(shuffled_windows)
}

# 2. convert a list of GRanges windows containing the ROH coverages (one
# GRanges object per individual) into one merged data frame (one column
# of coverage values for a windows per individual)
merge_windows <- function(windows) {
  # from each reshuffled individual, take just the vector of ROH coverages
  df <- lapply(windows, function(ind) {
    ind_df <- as_tibble(mcols(ind)["coverage"])
    names(ind_df) <- ind$sample[1]
    ind_df
  }) %>% do.call(cbind, .) %>%
    as_tibble

  return(df)
}

# 3. detect deserts in a given set of individuals (df will be a subset of
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

# Let's apply step 2. (merging the per-individual window-ROH coverages) and
# 3. (detecting deserts) on the original, non-shuffled data. This should
# produce counts of shared ancient-vs-modern deserts similar to the
# original analysis prior to the review.

# 2. merge the window-ROH hits into a single data frame
df <- merge_windows(coverages)
df

# 3. detect deserts in ancient and modern samples
deserts_ancient <- detect_deserts(df[, ancient_samples], cutoff = 0.05)
deserts_modern <- detect_deserts(df[, modern_samples], cutoff = 0.05)

sum(deserts_ancient)
# 2757
sum(deserts_modern)
# 228

# count which windows represent a shared A vs M desert
shared_deserts <- deserts_ancient & deserts_modern
sum(shared_deserts)
# --> 172

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
