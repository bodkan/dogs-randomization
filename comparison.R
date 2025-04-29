library(data.table)
library(qs2)
library(dplyr)
library(GenomicRanges)

# coverage filter
win_cov_df <- fread("data/dogs_allchrom_windows_cov_500kb.txt")
names(win_cov_df) <- c("chrom", "start", "end", ".", ".", ".", "coverage", ".", ".")
win_cov_df <- win_cov_df[, .(coverage = sum(coverage) / length(ancient_samples)), by = c("chrom", "start", "end")]
win_cov_df <- win_cov_df[-3920]
win_cov_df

sum_cov_df <- win_cov_df[, .(
  lower_cov = mean(coverage) - 2 * sd(coverage),
  upper_cov = mean(coverage) + 2 * sd(coverage)
)]
sum_cov_df


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

# sanity checking the winScan function on ancient individuals

library(windowscanr)
library(assertthat)

hom_sum <- aold %>% as.data.table() %>% dplyr::rename(CHR = chrom, BP = pos)
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

# running_roh_p <- running_roh_p[win_cov_df$valid]

roh_wins_df <- fread("roh_wins_anc_df.tsv")#[win_cov_df$valid]

running_roh_p
roh_wins_df

running_roh_p[UNAFF_mean < 0.05]
roh_wins_df[mean_ancient < 0.05]


# sanity checking the winScan function on modern individuals

mold <- fread("data/ref-panel_allchrom_sample-snp_filltags_filter_MAF_0.01_all_sites_hom_win_het_1_dogs.hom.summary.gz")
names(mold) <- c("chrom", ".", "pos", "AFF", "UNAFF")

hom_sum_modern <- mold %>% as.data.table() %>% dplyr::rename(CHR = chrom, BP = pos)
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

#running_roh_p_modern <- running_roh_p_modern[win_cov_df$valid]

roh_wins_mod_df <- fread("roh_wins_mod_df.tsv")#[win_cov_df$valid]

running_roh_p_modern
roh_wins_mod_df

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




# coverage stuff ------------------------------------------------------------------------------

library(stringr)

all_cov_window <- read.delim("data/dogs_allchrom_windows_cov_500kb.txt", header=FALSE)

all_cov_window$V1 <- as.integer(str_replace_all(all_cov_window$V1, 'chr', ''))
all_cov_window$region <- paste(all_cov_window$V1,":",all_cov_window$V2,"-",all_cov_window$V3, sep="")

# count number of bam samples used for the window coverage calculation (same number as the imputed dogs or wolves)
num_ind <- length(ancient_samples)

# Group by region and estimate across all samples the average mean depth for each window:
all_cov_window %>%
  group_by(V1, V2, V3) %>%
  summarize(Mean = sum(V7)/num_ind) -> test

#merge the window coverage with the ROH data (either imputed dogs or imputed wolves)
joinn <- left_join(running_roh_p, test, by=c("CHR" = "V1", "win_start" = "V2", "win_end" = "V3"))

#get mean and std
meann_wind <- mean(joinn$Mean)
stdd_wind <- 2*sd(joinn$Mean)
low_wind <- meann_wind-stdd_wind
high_wind <- meann_wind+stdd_wind


