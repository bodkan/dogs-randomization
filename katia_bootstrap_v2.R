library(data.table)
library(ggplot2)
library(dplyr)
library(tidyr)
library(tibble)
library(readr)
library(parallel)
library(GenomicRanges)

win_gr <- fread("cat dogs_allchrom_windows_cov_500kb.txt | cut -f1-3 | sort | uniq") %>%
  setNames(c("chrom", "start", "end")) %>%
  mutate(chrom = as.integer(gsub("chr", "", chrom))) %>%
  arrange(chrom, start, end) %>%
  makeGRangesFromDataFrame()
win_gr

modern_gr <- fread("ref-panel_allchrom_sample-snp_filltags_filter_MAF_0.01_all_sites_hom_win_het_1_dogs.hom")
modern_gr <- tibble(modern_gr) %>% select(sample = FID, chrom = CHR, start = POS1, end = POS2)
modern_gr$set <- "M"
modern_gr

ancient_gr <- fread("merged_phased_annotated.allchrom_MAF_0.01_recalibrated_INFO_0.8_all_sites_hom_win_het_1_dogs.hom")
ancient_gr <- tibble(ancient_gr) %>% select(sample = FID, chrom = CHR, start = POS1, end = POS2)
ancient_gr$set <- "A"
ancient_gr

roh_gr <- rbind(modern_gr, ancient_gr)
roh_gr <- makeGRangesFromDataFrame(roh_gr, keep.extra.columns = TRUE)
roh_gr

ancient_samples <- unique(roh_gr[roh_gr$set == "A"]$sample)
ancient_samples
modern_samples <- unique(roh_gr[roh_gr$set == "M"]$sample)
modern_samples

all_samples <- c(ancient_samples, modern_samples)
all_samples

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

# score each window in each sample in the original data as either hitting an ROH or not
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


simulate_replicate <- function(windows, rep_i) {
  # reshuffle windows in each individual
  shuffled_windows <- lapply(windows, function(ind) {
    ind$roh <- sample(ind$roh)
    ind
  })

  # from each reshuffled individual, take just the TRUE/FALSE vector of the windows hitting an ROH
  shuffled_deserts <- lapply(shuffled_windows, function(ind) {
    ind_df <- as_tibble(mcols(ind)["roh"])
    ind_df$roh <- as.integer(ind_df$roh)
    names(ind_df) <- ind$sample[1]
    ind_df
  }) %>% do.call(cbind, .) %>% as_tibble

  # split the windows into ancient and modern sets
  shuffled_ancient <- shuffled_deserts[, ancient_samples]
  shuffled_modern <- shuffled_deserts[, modern_samples]

  freq_ancient <- rowMeans(shuffled_ancient)
  freq_modern <- rowMeans(shuffled_modern)

  # flag each window as a desert or not
  cutoff <- 0.05
  deserts_ancient <- freq_ancient < cutoff
  deserts_modern <- freq_modern < cutoff

  # detect which reshuffled windows correspond to a shared A vs M desert
  deserts_shared <- deserts_ancient & deserts_modern

  df <- tibble(
    rep_i = rep_i,
    win_i = seq_along(deserts_shared),
    match = deserts_shared
  )

  return(df)
}

reps <- mclapply(1:100, function(rep_i) simulate_replicate(windows, rep_i = rep_i), mc.cores = detectCores()) %>%
  do.call(rbind, .) %>%
  mutate(match = as.integer(match))

reps

empirical_count <- 124

rep_counts <- group_by(reps, rep_i) %>% summarise(count = sum(match))
rep_counts

ggplot(rep_counts) + geom_histogram(aes(count))
