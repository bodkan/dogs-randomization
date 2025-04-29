suppressPackageStartupMessages({
library(data.table)
library(R.utils)
library(bit)
library(qs2)
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
modern_samples <- unique(roh_gr[roh_gr$set == "modern"]$sample)

all_samples <- c(ancient_samples, modern_samples)

# look-up table for detecting whether a given sample is ancient or modern
sample_lookup <- c(rep("ancient", length(ancient_samples)), rep("modern", length(modern_samples)))
names(sample_lookup) <- all_samples

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

# coordinates of windows
windows_gr <- fread("cat data/dogs_allchrom_windows_cov_500kb.txt | cut -f1-3 | sort | uniq") %>%
  setNames(c("chrom", "start", "end")) %>%
  mutate(chrom = as.integer(gsub("chr", "", chrom))) %>%
  arrange(chrom, start, end) %>%
  makeGRangesFromDataFrame()

cat("done.\n")
cat("---\n")

# remove windows which don't have any SNPs (there's only one of those)
hits <- findOverlaps(windows_gr, merged_sites_gr)
to_remove <- setdiff(seq_along(windows_gr), queryHits(hits))

cat("Removing windows without any SNPs:\n")
# 3920
print(as.data.frame(windows_gr)[to_remove, ])
cat("---\n")

windows_gr <- subsetByOverlaps(windows_gr, merged_sites_gr)

###############################################################
# processing site tables
###############################################################

# assign each site to a window (adding a window index column)
assign_sites <- function(sites_gr, windows_gr) {
  # find which window overlaps each SNP
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
                            dummy = TRUE,
                            modern = NA, ancient = NA, win_i = win_gr$win_i[1])
    } else {
      padding_gr <- NULL
    }
    win_gr$dummy <- FALSE
    win_gr <- do.call(c, list(win_gr, padding_gr))

    win_gr
  }) %>% GRangesList(compress = FALSE)

  # turn the list of GRanges object into a single GRanges table again
  unlist(snp_windows_padded)
}


cat("Processing sites... ")
s <- Sys.time()

if (!file.exists("sites_gr.qs")) {
  # assign a window number to each site
  sites_gr <- assign_sites(merged_sites_gr, windows_gr)

  # pad each window with dummy NA sites so that each window is of the same length
  # (and we can reshuffle windows in each individual during bootstrapping)
  sites_gr <- add_padding(sites_gr)

  qs_save(sites_gr, "sites_gr.qs")
} else {
  sites_gr <- qs_read("sites_gr.qs")
}

e <- Sys.time()
cat("done.\n")
print(e - s)
cat("---\n")

# sanity check -- all windows now must be of the same length (real SNPs as well
# as NA 'dummy SNPs')
expect_true(length(unique(table(sites_gr$win_i))) == 1)

###############################################################
# detecting sites overlapping ROH in each individual
###############################################################

detect_overlaps <- function(samples, sites_gr, roh_gr) {
  hits_inds <- mclapply(seq_along(samples), function(i) {
    ind <- samples[i]

    # extract the flag for whether that individual is ancient or modern
    set <- roh_gr[roh_gr$sample == ind]$set[1]

    #cat(sprintf("Processing %s (%s) [%s/%s]\n", ind, set, i, length(samples)))

    # take a respective subset of all sites (either ancient set or modern set),
    # which are not 'dummy sites' used for padding windows to the same length
    subset <- mcols(sites_gr)[[set]]
    subset <- subset & !is.na(subset)

    # detect which sites in this individual overlap with an ROH -- a site
    # in the overall data set will be either TRUE or FALSE (depending on whether
    # or not it overlaps or ROH) or NA (if its a site only available in the modern
    # set but an ancient individual is being currently processed) -- note that
    # the bit-wise representation of TRUE/FALSE here doesn't allow NA values though,
    # and the dummy-or-not mask will be handled via a dedicated data structure downstream
    roh_hits <- bit(length(subset))
    hits <- findOverlaps(sites_gr[subset], roh_gr[roh_gr$sample == ind])
    roh_hits[which(subset)][queryHits(hits)] <- TRUE

    # split the long bit vector into a list of bit vectors (one element for each window)
    roh_hits <- split(roh_hits, sites_gr$win_i) %>% lapply(as.bitwhich)
    roh_hits
  }, mc.cores = detectCores())
  # })
  names(hits_inds) <- samples

  return(hits_inds)
}

cat("Assigning sites to ROH segments... ")
s <- Sys.time()

if (!file.exists("roh_overlaps.qs")) {
  roh_overlaps <- detect_overlaps(all_samples, sites_gr, roh_gr)
  qs_save(roh_overlaps, "roh_overlaps.qs")
} else {
  roh_overlaps <- qs_read("roh_overlaps.qs")
}

e <- Sys.time()
cat("done.\n")
print(e - s)
cat("---\n")

cat("Generating bit masks for each window... ")

# convert dummy/modern/ancient TRUE/FALSE columns from a data frame format (inefficient
# for massive-scale parallel processing during bootstrap, even when stored as data.table)
# to lists of window-based bit-masks -- note that we have to store the dummy-or-not status
# for each site as a TRUE/FALSE bit mask because bit vectors don't allow NA values (whether
# or not is any site in any individual a dummy site would be indicated by an NA value in
# that individual's data.table column)
masks <- as.data.table(sites_gr)[, .(index = 1:.N, win_i, dummy, ancient, modern)] %>%
  split(., .$win_i) %>%
  lapply(function(win_df) {
    with(win_df,
         list(
           dummy = as.bit(dummy),
           ancient = as.bit(ancient),
           modern = as.bit(modern)
          ))
  })

cat("done.\n")
cat("---\n")

###############################################################
# sanity checking after ROH-assignment, and bit-conversion
###############################################################

# bit vectors of non-dummy modern sites have the same number as the original table of modern sites
expect_true(
  sum(sapply(masks, function(win) sum(!win$dummy & win$modern))) ==
    length(sites_gr[!sites_gr$dummy & sites_gr$modern])
)

# bit vectors of non-dummy ancient sites have the same number as the original table of ancient sites
expect_true(
  sum(sapply(masks, function(win) sum(!win$dummy & win$ancient))) ==
    length(sites_gr[!sites_gr$dummy & sites_gr$ancient])
)

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
windows_coverage <- function(samples, sample_lookup, roh_overlaps, masks) {
  mclapply(samples, function(ind) {
    # iterate over each window in a given individual...
    sapply(seq_len(length(masks)), function(win_i) {
      # ... extract the ROH-hit TRUE/FALSE status of each site in the window ...
      roh_sites <- roh_overlaps[[ind]][[win_i]]
      # ... extract the filter of the dummy status of each site in the window...
      dummy <- masks[[win_i]]$dummy
      # ... and the ancient or modern status of each site in the window ...
      set <- masks[[win_i]][[sample_lookup[[ind]]]]

      # the proportion of sites overlapping ROH in this individual's window is the number
      # of sites which hit some ROH divided by the total number of available sites
      # (the number of non-dummy sites which are also among the ancient or modern set)
      sum(roh_sites) / sum(!dummy & set)
    })
  }, mc.cores = detectCores()) %>%
    as.data.table %>%
    setNames(samples)
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

# compute ROH/SNP coverage in each window

cat("Computing ROH coverage of each SNP in each window... ")
s <- Sys.time()

mean_win_df <- windows_coverage(all_samples, sample_lookup, roh_overlaps, masks)

e <- Sys.time()
cat("done.\n")
print(e - s)
cat("---\n")

# assign desert status to each window (TRUE or FALSE)
deserts_ancient <- detect_deserts(mean_win_df[, ..ancient_samples], cutoff = 0.05)
deserts_modern <- detect_deserts(mean_win_df[, ..modern_samples], cutoff = 0.05)

# count the deserts
cat("Number of deserts in ancient individuals:", sum(deserts_ancient), "\n")
# 2677
cat("Number of deserts in modern individuals:", sum(deserts_modern), "\n")
# 212
deserts_shared <- deserts_ancient & deserts_modern
cat("Number of shared deserts:", sum(deserts_shared), "\n")
# 159
cat("---\n")

# for debugging purposes, add ROH frequencies for each window (in ancient and
# modern individuals) to the original table of windows for easier reference
original_win <- windows_gr

original_win$mean_ancient <- rowMeans(mean_win_df[, .SD, .SDcols = ancient_samples])
original_win$mean_modern <- rowMeans(mean_win_df[, .SD, .SDcols = modern_samples])
original_win$desert <- deserts_shared

# just the deserts in ancient individuals
# original_win[original_win$mean_ancient < 0.05]

# just the deserts in modern individuals
# original_win[original_win$mean_modern < 0.05]

# deserts in both ancient and modern individuals
# original_win[original_win$desert]

# save coordinates of shared ancient & modern deserts
deserts <- original_win[original_win$desert]
deserts_df <- as.data.table(deserts)[, .(chrom = seqnames, start, end, mean_ancient, mean_modern, desert)]

cat("Coordinates of shared deserts:\n\n")

print(as.data.frame(deserts_df))

cat("---\n")

fwrite(deserts_df, "deserts.tsv", sep = "\t", row.names = FALSE)

###############################################################
# Get data for the comparison of the pre-review winScan code
###############################################################

# save window information for ancient individuals
roh_wins_anc_df <- as.data.table(original_win)
roh_wins_anc_df <- roh_wins_anc_df[, .(chrom = seqnames, start, end, mean_ancient)]
fwrite(roh_wins_anc_df, "roh_wins_anc_df.tsv", sep = "\t", row.names = FALSE)

# save window information for modern individuals
roh_wins_mod_df <- as.data.table(original_win)
roh_wins_mod_df <- roh_wins_mod_df[, .(chrom = seqnames, start, end, mean_modern)]
fwrite(roh_wins_mod_df, "roh_wins_mod_df.tsv", sep = "\t", row.names = FALSE)

# ancient samples
ancient_states_df <-
  mclapply(ancient_samples, function(ind) {
    set <- sample_lookup[ind]
    lapply(seq_along(masks), function(win_i) {
      available_sites <- which(as.logical(!masks[[win_i]]$dummy & masks[[win_i]][[set]]))
      roh_overlaps[[ind]][[win_i]][available_sites]
    }) %>% unlist
  }) %>% as.data.table %>% setNames(ancient_samples)

ancient_coord_df <- as.data.table(sites_gr)[
  (ancient) & !(dummy), .(chrom = seqnames, pos = start, win_i, modern, ancient, dummy)
]

ancient_states_df <- cbind(ancient_coord_df, ancient_states_df)
qs_save(ancient_states_df, "ancient_states_df.qs")

# modern samples
modern_states_df <-
  mclapply(modern_samples, function(ind) {
    set <- sample_lookup[ind]
    lapply(seq_along(masks), function(win_i) {
      available_sites <- which(as.logical(!masks[[win_i]]$dummy & masks[[win_i]][[set]]))
      roh_overlaps[[ind]][[win_i]][available_sites]
    }) %>% unlist
  }) %>% as.data.table %>% setNames(modern_samples)

modern_coord_df <- as.data.table(sites_gr)[
  (modern) & !(dummy), .(chrom = seqnames, pos = start, win_i, modern, ancient, dummy)
]

modern_states_df <- cbind(modern_coord_df, modern_states_df)
qs_save(modern_states_df, "modern_states_df.qs")

###############################################################
# Comparison of pre-review and post-review desert windows
###############################################################

# read coordinates of new deserts
deserts <-
  fread("deserts.tsv")[(desert)] %>%
  makeGRangesFromDataFrame(keep.extra.columns = TRUE)

# read coordinates of deserts from the pre-review manuscript
old_deserts <-
  fread("data/imputed_modern_window_bed_0.01_INFO_0.8_all_sites_hom_win_het_1_dogs_deserts.bed",
        col.names = c("chrom", "start", "end")) %>%
  makeGRangesFromDataFrame(keep.extra.columns = TRUE)

cat("Number of new deserts:", length(deserts), "\n")
# 159

cat("Number of original deserts:", length(old_deserts), "\n")
# 124

# detect deserts which were present both pre- and post-review
overlapping_deserts <- findOverlaps(deserts, old_deserts, type = "equal")
cat("Number of deserts shared between both sets:", length(overlapping_deserts), "\n")
# 124 -- old desert windows are all within the new desert windows

cat("---\n")

###############################################################
# Bootstrapping
###############################################################

# shuffle windows in a given individual (keeping the sites sitting on the window together)
# -- the windows are stored as lists of bit-encoded vectors, each element of a vector
#    indicating whether or not does a site overlap an ROH or not in this individual
#    (TRUE and FALSE states encoded as sparse bit vectors for maximum memory efficiency)
shuffle_one <- function(windows, shuffled_indices = NULL) {
  if (is.null(shuffled_indices))
    shuffled_indices <- sample(seq_along(windows))
  windows[shuffled_indices]
}

# shuffle windows (and, therefore, sites) in all given individuals
shuffle_samples <- function(samples, roh_overlaps) {
  lapply(roh_overlaps, shuffle_one)
}

cat("Testing that reshuffling works correctly and uses constant amount of memory... ")
s <- Sys.time()

{
  # pick a single individual
  ind <- sample(all_samples, 1)
  windows <- roh_overlaps[[ind]]
  # reshuffle its sites/windows
  shuffled_indices <- sample(seq_along(windows))
  shuffled_windows <- windows[shuffled_indices]

  # check that the pointers to bit-encoded vectors (one for each window) point
  # to the same address in memory (proving that the reshuffling works correctly
  # and that no memory has been copied in the process, only the reference counter
  # has been increased by one
  original_pointers <- sapply(seq_along(windows),
                              function(i) rlang::obj_address(windows[[i]]))
  shuffled_pointers <- sapply(seq_along(windows),
                              function(i) rlang::obj_address(shuffled_windows[[which(shuffled_indices == i)]]))
  expect_equal(original_pointers, shuffled_pointers)
}

cat("done.\n")
e <- Sys.time()
print(e - s)
cat("---\n")

cat("Testing that identity permutation results in identical desert counts... ")
s <- Sys.time()

{
  # "reshuffle" windows in each individual to their original positions
  shuffled_overlaps <- lapply(roh_overlaps, function(windows) {
    original_indices <- seq_along(windows)
    shuffle_one(windows, shuffled_indices = original_indices)
  })
  # compute the mean ROH/site coverage of each window
  mean_win_df <- windows_coverage(all_samples, sample_lookup, shuffled_overlaps, masks)
  # detect windows which are ROH deserts
  deserts_ancient <- detect_deserts(mean_win_df[, ..ancient_samples], cutoff = 0.05)
  deserts_modern <- detect_deserts(mean_win_df[, ..modern_samples], cutoff = 0.05)
  # make sure the identity permutations of windows give the original result
  expect_true(sum(deserts_ancient & deserts_modern) == 159)
}

cat("done.\n")
e <- Sys.time()
print(e - s)
cat("---\n")

# run a single bootstrap replicate of the desert inference
run_replicate <- function(rep_i, roh_overlaps, masks) {
  cat(paste0("Running replicate #", rep_i, "... "))

  # reshuffle windows in each individual
  shuffled_overlaps <- shuffle_samples(all_samples, roh_overlaps)

  # compute SNP-ROH coverage in each window
  mean_win_df <- windows_coverage(all_samples, sample_lookup, shuffled_overlaps, masks)

  # detect which ancient vs modern deserts are shared (producing TRUE/FALSE vectors)
  deserts_ancient <- detect_deserts(mean_win_df[, ..ancient_samples], cutoff = 0.05)
  deserts_modern <- detect_deserts(mean_win_df[, ..modern_samples], cutoff = 0.05)

  # detect which reshuffled windows correspond to a shared A vs M desert
  deserts_shared <- deserts_ancient & deserts_modern

  # return the result in a tidy form -- a data frame with these columns:
  #   - rep_i: replicate number
  #   - win_i: index of a window
  #   - desert_ancient: whether that window is an 'ancient desert' (1) or not (0)
  #   - desert_modern: whether that window is an 'modern desert' (1) or not (0)
  #   - desert_shared: whether that window is a desert in both ancient and modern dogs
  result <- data.table(
    rep_i = rep_i,
    win_i = seq_along(deserts_shared),
    desert_ancient = as.integer(deserts_ancient),
    desert_modern = as.integer(deserts_modern),
    desert_shared = as.integer(deserts_shared)
  )

  cat("done.\n")

  return(result)
}

cat("Starting the bootstrap procedure... \n")
s <- Sys.time()

n_reps <- 100

if (!file.exists("replicates_df.qs")) {
  # run 100 reshuffling bootstrap replicates
  replicates_df <-
    lapply(seq_len(n_reps), function(rep_i) run_replicate(rep_i, roh_overlaps, masks)) %>%
    do.call(rbind, .)

  qs_save(replicates_df, "replicates_df.qs")

} else {
  replicates_df <- qs_read("replicates_df.qs")
}

cat("Bootstrap procedure finished... \n")
e <- Sys.time()
print(e - s)

cat("Time for a single iteration:\n")
print((e - s) / n_reps)
cat("---\n")

###############################################################
# Putting a p-value on the result
###############################################################

# the observed number of shared deserts
observed_count <- length(deserts)

# desert numbers observed in each bootstrap iteration
bootstrap_counts <- replicates_df[, .(desert_count = sum(desert_shared)), by = rep_i]

# compute and plot the empirical CDF
e <- ecdf(bootstrap_counts$desert_count)

pdf("ecdf.pdf", width = 10, height = 7)

plot(e, xlim = c(0, max(bootstrap_counts$desert_count)),
     main = "ECDF of the bootstrapped shared desert counts")

plot(e, xlim = c(0, max(bootstrap_counts$desert_count, observed_count)),
     main = "ECDF of the bootstrapped shared desert counts (along with the observed count)")
abline(v = observed_count, col = "red", lty = 2)
legend(x = 120, y = 0.95, "observed count", fill = "red")

invisible(dev.off())

# what's the probability of observing a value as extreme (or more extreme)
# than the value we observed?
p_value <- 1 - e(observed_count)
cat("Probability of observing the same (or larger) number of shared deserts:", p_value)

# histogram of the bootstrap counts along with the observed value
pdf("bootstrap.pdf", width = 8, height = 5)

ggplot(bootstrap_counts) +
  geom_histogram(aes(desert_count, color = "bootstrap"), binwidth = 1, fill = "darkgray") +
  geom_vline(aes(xintercept = observed_count, color = "observed"), linetype = "dashed") +
  coord_cartesian(xlim = c(0, max(bootstrap_counts$desert_count, observed_count))) +
  guides(color = guide_legend("")) +
  scale_color_manual(values = c("observed" = "red", "bootstrap" = "darkgray")) +
  annotate("label", x = 140, y = 38, label = sprintf("p-value = %s", p_value), color = "red") +
  labs(x = "number of shared deserts", y = "number of bootstrap replicates") +
  theme_minimal()

invisible(dev.off())
