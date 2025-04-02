library(data.table)
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
modern_gr$set <- "M"
modern_gr

# ROH in ancient samples
ancient_gr <- fread("data/merged_phased_annotated.allchrom_MAF_0.01_recalibrated_INFO_0.8_all_sites_hom_win_het_1_dogs.hom")
ancient_gr <- tibble(ancient_gr) %>% select(sample = FID, chrom = CHR, start = POS1, end = POS2)
ancient_gr$set <- "A"
ancient_gr

# combine ROH sets into a single GRanges object
roh_gr <- rbind(modern_gr, ancient_gr)
roh_gr <- makeGRangesFromDataFrame(roh_gr, keep.extra.columns = TRUE)
roh_gr

# extract sample names for later use
ancient_samples <- unique(roh_gr[roh_gr$set == "A"]$sample)
ancient_samples
modern_samples <- unique(roh_gr[roh_gr$set == "M"]$sample)
modern_samples

all_samples <- c(ancient_samples, modern_samples)
all_samples

###############################################################
# detecting windows overlapping ROH in each individual
###############################################################

# for a given individual ind, score each window as TRUE or
# FALSE depending on whether it overlaps and ROH in that sample
flag_windows <- function(ind, roh_gr, win_gr) {
  # extract ROHs only in the given individual
  ind_gr <- roh_gr[roh_gr$sample == ind]

  # create a window set just for that individual (and initialize
  # with additional column)
  ind_win_gr <- win_gr
  ind_win_gr$sample <- ind
  ind_win_gr$set <- ind_gr$set[1]
  ind_win_gr$roh <- FALSE

  # find which windows hit at least one ROH in that individual
  hits <- findOverlaps(ind_win_gr, ind_gr)

  # annotate each window with TRUE if that window hit some ROH
  ind_win_gr[queryHits(hits)]$roh <- TRUE

  # return the ROH-hitting windows in that individual
  return(ind_win_gr)
}

# score each original window in each sample in the original data
# as either hitting an ROH or not
windows <- mclapply(all_samples, function(ind) flag_windows(ind, roh_gr, win_gr),
                    mc.cores = detectCores())

# a couple of sanity checks:

# 1. the number of windows in each sample must be the same
sapply(windows, length) %>% unique

# 2. the number of ancient and modern samples must be matching the original data
sapply(windows, function(ind) unique(ind$set)) %>% table

# 3. just a visual test that the # of ROHs in ancient vs modern match our expectation
roh_counts <- lapply(windows, function(ind) {
  data.frame(set = ind$set[1], ind = ind$sample[1], roh_count = sum(ind$roh))
}) %>% do.call(rbind, .)
ggplot(roh_counts) +
  geom_jitter(aes(set, roh_count, color = set))


###############################################################
# Bootstrap components (1, 2, 3)
###############################################################

# 1. take a list of windows for each individual as a GRanges object and
# shuffle just the `roh` column in each of them, returning back
# GRanges object modified in this way
shuffle_windows <- function(windows) {
  shuffled_windows <- lapply(windows, function(ind) {
    ind$roh <- sample(ind$roh)
    ind
  })

  return(shuffled_windows)
}

# 2. convert a list of GRanges windows (one GRanges object per individual) into
# one merged data frame (one column of TRUE/FALSE for a windows per individual)
merge_windows <- function(windows) {
  # from each reshuffled individual, take just the TRUE/FALSE vector
  # of the windows hitting an ROH
  df <- lapply(windows, function(ind) {
    ind_df <- as_tibble(mcols(ind)["roh"])
    ind_df$roh <- as.integer(ind_df$roh)
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
  freq <- rowMeans(df)

  # flag each window as a desert or not
  cutoff <- 0.05
  deserts <- freq < cutoff

  return(deserts)
}

###############################################################
# sanity checking the desert inference in the original data
###############################################################

# Let's apply step 2. (merging the per-individual window-ROH hits) and
# 3. (detecting deserts) on the original, non-shuffled data. This should
# produce counts of shared ancient-vs-modern deserts similar to the
# original analysis prior to the review.

# 2. merge the window-ROH hits into a single data frame
df <- merge_windows(windows)
df

# 3. detect deserts in ancient and modern samples
deserts_ancient <- detect_deserts(df[, ancient_samples], cutoff = 0.05)
deserts_modern <- detect_deserts(df[, modern_samples], cutoff = 0.05)

# count which windows represent a shared A vs M desert
shared_deserts <- deserts_ancient & deserts_modern
sum(shared_deserts)
# --> 53

# For debugging purposes, add frequencies of ROH overlapping each
# desert (in ancient and modern individuals) to the original table of windows
original_win <- win_gr
original_win$freq_ancient <- rowMeans(df[, ancient_samples])
original_win$freq_modern <- rowMeans(df[, modern_samples])
original_win$desert <- shared_deserts
original_win$desert %>% sum
original_win

# just the deserts in modern individuals
original_win[original_win$freq_modern < 0.05]

# just the deserts in ancient individuals
original_win[original_win$freq_ancient < 0.05]

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
