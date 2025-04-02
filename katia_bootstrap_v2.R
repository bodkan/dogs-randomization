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
win_gr <- fread("cat dogs_allchrom_windows_cov_500kb.txt | cut -f1-3 | sort | uniq") %>%
  setNames(c("chrom", "start", "end")) %>%
  mutate(chrom = as.integer(gsub("chr", "", chrom))) %>%
  arrange(chrom, start, end) %>%
  makeGRangesFromDataFrame()
win_gr

# ROH in modern samples
modern_gr <- fread("ref-panel_allchrom_sample-snp_filltags_filter_MAF_0.01_all_sites_hom_win_het_1_dogs.hom")
modern_gr <- tibble(modern_gr) %>% select(sample = FID, chrom = CHR, start = POS1, end = POS2)
modern_gr$set <- "M"
modern_gr

# ROH in ancient samples
ancient_gr <- fread("merged_phased_annotated.allchrom_MAF_0.01_recalibrated_INFO_0.8_all_sites_hom_win_het_1_dogs.hom")
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
# step 1. -- detecting which windows overlap some ROH
###############################################################

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

# score each original window in each sample in the original data as either hitting an ROH or not
windows <- mclapply(all_samples, function(ind) flag_windows(ind, roh_gr, win_gr), mc.cores = detectCores())

# a couple of sanity checks
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


# 1. take a list of windows for each individual as a GRanges object and shuffle just
# the `roh` column in each of them
shuffle_windows <- function(windows) {
  shuffled_windows <- lapply(windows, function(ind) {
    ind$roh <- sample(ind$roh)
    ind
  })

  return(shuffled_windows)
}

# 2. convert the list of GRanges windows (one element per individual) into a merged
# data frame (one column of TRUE/FALSE for a windows per individual)
merge_windows <- function(windows) {
  # from each reshuffled individual, take just the TRUE/FALSE vector of the windows hitting an ROH
  df <- lapply(windows, function(ind) {
    ind_df <- as_tibble(mcols(ind)["roh"])
    ind_df$roh <- as.integer(ind_df$roh)
    names(ind_df) <- ind$sample[1]
    ind_df
  }) %>% do.call(cbind, .) %>%
    as_tibble

  return(df)
}

# 3. detect deserts in a given set of individuals' windows
detect_deserts <- function(windows, cutoff) {
  freq <- rowMeans(windows)

  # flag each window as a desert or not
  cutoff <- 0.05
  deserts <- freq < cutoff

  return(deserts)
}

###############################################################
# sanity checking the desert inference in the original data
###############################################################

# test that detecting shared deserts works on the original (non-shuffled!) windows
# as we expect
original_merged <- merge_windows(windows)
original_merged

deserts_ancient <- detect_deserts(original_merged[, ancient_samples], cutoff = 0.05)
deserts_modern <- detect_deserts(original_merged[, modern_samples], cutoff = 0.05)

# detect which reshuffled windows correspond to a shared A vs M desert
shared_deserts <- deserts_ancient & deserts_modern
sum(shared_deserts)

original_win <- win_gr
original_win$freq_ancient <- rowMeans(original_merged[, ancient_samples])
original_win$freq_modern <- rowMeans(original_merged[, modern_samples])
original_win$desert <- shared_deserts
original_win$desert %>% sum
original_win

# just the deserts
original_win[original_win$freq_modern < 0.05]

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

reps <- mclapply(1:10, function(rep_i) run_replicate(windows, rep_i = rep_i),
                 mc.cores = detectCores()) %>%
  do.call(rbind, .) %>%
  mutate(match = as.integer(match))

reps

empirical_count <- 124

rep_counts <- group_by(reps, rep_i) %>% summarise(count = sum(match))
rep_counts

ggplot(rep_counts) + geom_histogram(aes(count))

counts <- rpois(1000, 5)
counts

plot(density(counts))
abline(v = 10)

e <- ecdf(counts)
1 - e(10)


counts <- rep(0, 100)
e <- ecdf(counts)
plot(e)
1 - e(10)
