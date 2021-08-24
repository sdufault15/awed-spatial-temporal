# Sent to the cluster for efficiency

library(tidyr)
library(dplyr)
library(geodist)
library(furrr)
library(purrr)
library(permute)

load("2020-12-08_work-ts-spat.RData") # Data

########################
# Analysis
########################
source("permutation-function.R")
source("binary-matrix-function_odds.R")

bootstrap_vcd <- work_ts_spat %>% 
  filter(serotype != "unk_serotype")


# Increasing Radius
r_upper <- as.vector(c(50,seq(100, 4000, by = 100)), mode = "list")
r_lower <- as.vector(rep(0,length(r_upper)), mode = "list")
names(r_upper) <- r_upper
names(r_lower) <- r_lower
perms <- as.vector(1:1000, mode = "list")


plan(multisession)
permuted_radii <- future_map_dfr(perms, ~permutation_function(df = bootstrap_vcd,
                                                                 r1 = r_lower,
                                                                 r2 = r_upper,
                                                                 t = 30),
                                     .id = "permutation")

# Concentric circles
r_lower = lapply(r_upper, function(x) {ifelse(x >100, x - 100, x-50)})

plan(multisession)
permuted_concentric <- future_map_dfr(perms, ~permutation_function(df = bootstrap_vcd,
                                                                      r1 = r_lower,
                                                                      r2 = r_upper,
                                                                      t = 30),
                                          .id = "permutation")

output_permuted <- list(permuted_radii = permuted_radii,
               permuted_concentric = permuted_concentric)

save(output_permuted, 
     file = "overall-permutation-output.RData")
