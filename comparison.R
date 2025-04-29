library(data.table)
library(qs2)
library(dplyr)
library(GenomicRanges)

# coverage filter
win_cov_df <- fread("data/dogs_allchrom_windows_cov_500kb.txt")
names(win_cov_df) <- c("chrom", "start", "end", ".", ".", ".", "coverage", ".", ".")
win_cov_df <- win_cov_df[, .(coverage = sum(coverage) / length(ancient_samples)), by = c("chrom", "start", "end")]

sum_cov_df <- win_cov_df[, .(mean = mean(coverage), std = 2 * sd(coverage))]

win_cov_df[, valid := coverage > sum_cov_df$mean - 2*sum_cov_df$std & coverage < sum_cov_df$mean + 2*sum_cov_df$std]
win_cov_df <- win_cov_df[-3920]
win_cov_df

# detecting SNP mismatches

aold <- fread("data/merged_phased_annotated.allchrom_MAF_0.01_recalibrated_INFO_0.8_all_sites_hom_win_het_1_dogs.hom.summary.gz")
names(aold) <- c("chrom", ".", "pos", "AFF", "UNAFF")

anew <- qs_read("ancient_states_df.qs")
ancient_samples <- setdiff(colnames(anew), c("chrom", "pos", "win_i", "modern", "ancient", "dummy"))
anew <- anew[, .(chrom, pos, n = rowSums(.SD)), .SDcols = ancient_samples]
anew[, chrom := as.integer(chrom)]

aold[UNAFF > 0]
anew[n > 0]

matches <- aold[anew, on = .(chrom, pos)]

matches[, mean(UNAFF == n)]

# sites present in the PLINK SNPs but missing in my filtered sites
mismatches <- anew[aold, on = .(chrom, pos)]
mismatches

mismatches[UNAFF > 0]

# histogram of the UNAFF ROH counts
table(mismatches[is.na(n), UNAFF])

# ROH frequencies
table(mismatches[is.na(n), UNAFF] / length(ancient_samples))

# sites above 5% frequency
mismatches[is.na(n) & UNAFF / length(ancient_samples) > 0.05]
# all within 380kb...
mismatches[is.na(n) & UNAFF / length(ancient_samples) > 0.05, .(min(pos), max(pos), max(pos) - min(pos))]
# ... in a single spot on chromosome 27
mismatches[is.na(n) & UNAFF / length(ancient_samples) > 0.05,
           plot(pos, rep(1, length(pos)), xlim = c(0, 50e6))]



# finding the missing SNPs
aold <- makeGRangesFromDataFrame(aold, keep.extra.columns = TRUE, start.field = "pos", end.field = "pos")
aold

anew <- makeGRangesFromDataFrame(anew, keep.extra.columns = TRUE, start.field = "pos", end.field = "pos")
anew

# coordinates of windows
windows_gr <- fread("cat data/dogs_allchrom_windows_cov_500kb.txt | cut -f1-3 | sort | uniq") %>%
  setNames(c("chrom", "start", "end")) %>%
  mutate(chrom = as.integer(gsub("chr", "", chrom))) %>%
  arrange(chrom, start, end) %>%
  makeGRangesFromDataFrame()

subsetByOverlaps(aold, windows_gr)
subsetByOverlaps(anew, windows_gr)


# sanity checking the winScan function on ancient individuals

library(windowscanr)
library(assertthat)

hom_sum <- aold %>% as.data.table() %>% dplyr::rename(CHR = seqnames, BP = start)
hom_sum[, index := 1:nrow(.SD)]

running_roh <- winScan(x = hom_sum,
                       groups = "CHR",
                       position = "BP",
                       values = "UNAFF",
                       win_size = 500000,
                       win_step = 500000,
                       funs = c("mean"),
                       cores = 8)

running_roh <- running_roh[-3920]

running_roh %>%
  mutate(UNAFF_mean = UNAFF_mean/length(ancient_samples)) %>%
  filter(UNAFF_n > 0) -> running_roh_p

running_roh_p$region <- paste(running_roh_p$CHR,":",running_roh_p$win_start,"-",running_roh_p$win_end, sep="")

setDT(running_roh_p)

running_roh_p <- running_roh_p[win_cov_df$valid]

roh_wins_df <- fread("roh_wins_anc_df.tsv")[win_cov_df$valid]

running_roh_p
roh_wins_df

running_roh_p[UNAFF_mean < 0.05]
roh_wins_df[mean_ancient < 0.05]


# sanity checking the winScan function on modern individuals

mold <- fread("data/ref-panel_allchrom_sample-snp_filltags_filter_MAF_0.01_all_sites_hom_win_het_1_dogs.hom.summary.gz")
names(mold) <- c("chrom", ".", "pos", "AFF", "UNAFF")

mold <- makeGRangesFromDataFrame(mold, keep.extra.columns = TRUE, start.field = "pos", end.field = "pos")
mold

hom_sum_modern <- mold %>% as.data.table() %>% dplyr::rename(CHR = seqnames, BP = start)
hom_sum_modern[, index := 1:nrow(.SD)]

running_roh_modern <- winScan(x = hom_sum_modern,
                              groups = "CHR",
                              position = "BP",
                              values = "UNAFF",
                              win_size = 500000,
                              win_step = 500000,
                              funs = c("mean"),
                              cores = 8)

running_roh_modern <- running_roh_modern[-3920]

# take for df:
running_roh_modern$region <- paste(running_roh_modern$CHR,":",running_roh_modern$win_start,"-",running_roh_modern$win_end, sep="")

# remove windows without snps
running_roh_modern %>%
  mutate(UNAFF_mean = UNAFF_mean/length(modern_samples)) %>%
  filter(UNAFF_n > 0) -> running_roh_p_modern

setDT(running_roh_p_modern)

running_roh_p_modern <- running_roh_p_modern[win_cov_df$valid]

roh_wins_mod_df <- fread("roh_wins_mod_df.tsv")[win_cov_df$valid]

running_roh_p_modern
roh_wins_mod_df

plot(running_roh_p_modern$UNAFF_mean, roh_wins_mod_df$mean_modern)

running_roh_p_modern[UNAFF_mean < 0.05]
roh_wins_mod_df[mean_modern < 0.05]



##########
# testing the intersect

old_ancient <- makeGRangesFromDataFrame(running_roh_p, seqnames.field = "CHR",
                                        start.field = "win_start", end.field = "win_end", keep.extra.columns = TRUE)
old_modern <- makeGRangesFromDataFrame(running_roh_p_modern, seqnames.field = "CHR",
                                        start.field = "win_start", end.field = "win_end", keep.extra.columns = TRUE)

(old_ancient$UNAFF_mean < 0.05) %>% sum
(old_modern$UNAFF_mean < 0.05) %>% sum
(old_ancient$UNAFF_mean < 0.05 & old_modern$UNAFF_mean < 0.05) %>% sum



