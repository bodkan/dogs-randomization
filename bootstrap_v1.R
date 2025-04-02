# Just a sketch of an idea for individual-based bootstrapping on a toy fake dataset

library(GenomicRanges)
library(boot)
library(dplyr)
library(tidyr)
library(tibble)
library(truncnorm)

ind1_roh <- data.frame(ind = "ind1", chrom = "chr1", start = c(1, 6, 10), end = c(4, 8, 13), set = "M")
ind2_roh <- data.frame(ind = "ind2", chrom = "chr1", start = c(8, 12), end = c(13, 13), set = "M")
ind3_roh <- data.frame(ind = "ind3", chrom = "chr1", start = c(7, 14), end = c(13, 16), set = "M")
ind4_roh <- data.frame(ind = "ind4", chrom = "chr1", start = c(6, 12), end = c(8, 14), set = "A")

original_data <- rbind(ind1_roh, ind2_roh, ind3_roh, ind4_roh)
original_data

# this is what we want to bootstrap over (i.e. individuals)
split(original_data, original_data$ind)

# we can do this and still keep the data as a data frame using nest
data <- nest(original_data, roh = c(chrom, start, end))
data



# summary statistic function
summary_fun <- function(data, ...) {
  unnested_data <- unnest(data, "roh")

  coverage_A <- rtruncnorm(1, a = 0, mean = 2)
  coverage_M <- rtruncnorm(1, a = 0, mean = 80)
  coverage_ratio <- coverage_A / coverage_M
  desert_A <- sample(c(0, 1), 10, replace = TRUE)
  desert_M <- sample(c(0, 1), 10, replace = TRUE)

  return(c(
    coverage_A = coverage_A,   # stat 1
    coverage_M = coverage_M,   # stat 2
    coverage_ratio = coverage_ratio, # stat 3
    desert_sum_A = sum(desert_A),    # stat 4
    desert_sum_M = sum(desert_M),    # stat 5
    desert_sharing = mean(desert_A == 1 & desert_M == 1), # stat 6
    desert_A = desert_A,
    desert_M = desert_M
  ))
}

# function to reshuffle data in each bootstrap iteration
reshuffle_samples <- function(data, ...) {
  anc_data <- data[data$set == "A", ]
  mod_data <- data[data$set == "M", ]

  anc_indices <- sample(1:nrow(anc_data), nrow(anc_data), replace = TRUE)
  mod_indices <- sample(1:nrow(mod_data), nrow(mod_data), replace = TRUE)

  res_anc_data <- anc_data[anc_indices, ]
  res_mod_data <- mod_data[mod_indices, ]

  rbind(res_anc_data, res_mod_data)
}

# run a single fake bootstrap iteration for testing
shuffled_data <- reshuffle_samples(data)
summary_fun(shuffled_data)

# run the full bootstrap
res <- boot(data, summary_fun, R = 100, sim = "parametric", ran.gen = reshuffle_samples)
summary(res)
res

# normal statistics
plot(res, index = 1) # coverage A
plot(res, index = 2) # coverage M
plot(res, index = 3) # coverage ratio
plot(res, index = 4) # desert count in A
plot(res, index = 5) # desert count in M
plot(res, index = 6) # desert sharing proportion

# A vs M "desert" sharing
anc <- res$t[, 7        :  (7 + 10 - 1)]; anc
mod <- res$t[, (7 + 10) : (7 + 10 + 10 - 1)]; mod

# proportion of sharing of deserts
(anc == 1 & mod & 1)
colMeans((anc == 1 & mod & 1))



