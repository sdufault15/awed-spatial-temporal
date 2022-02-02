# Sent to the cluster for efficiency

library(tidyr)
library(dplyr)
library(geodist)
library(furrr)
library(purrr)
library(permute)

# load("data/2020-12-08_work-ts-spat.RData") # Data
load("2020-12-08_work-ts-spat.RData") # Data

########################
# Analysis
########################
source("permutation-function.R")
source("binary-matrix-function_odds.R")
cluster_wrap_function <- function(df, clusters, perms, r1, r2, t = 30){
  df_cluster <- df[df$cluster == clusters,]
  
  output <- map_dfr(perms, 
                    ~permutation_function(df = df_cluster,
                                          r1 = r1, 
                                          r2 = r2,
                                          t = t),
                    .id = "permutation")
  return(output)
}

perm_vcd <- work_ts_spat %>% 
  filter(serotype != "unk_serotype")

# Increasing Radius
r_upper <- as.vector(c(50,seq(100, 1000, by = 100)), mode = "list")
r_lower <- as.vector(rep(0,length(r_upper)), mode = "list")
names(r_upper) <- r_upper
names(r_lower) <- r_lower
perms <- as.vector(1:1000, mode = "list")


clusters <- as.list(1:24)
names(clusters) <- clusters
plan(multisession)

permuted_radii_t7 <- future_map_dfr(clusters, ~cluster_wrap_function(
  clusters = .x,
  df = perm_vcd,
  perms = perms,
  r1 = r_lower,
  r2 = r_upper,
  t = 7),
  .id = "clusterid"
)

permuted_radii_t14 <- future_map_dfr(clusters, ~cluster_wrap_function(
  clusters = .x,
  df = perm_vcd,
  perms = perms,
  r1 = r_lower,
  r2 = r_upper,
  t = 14),
  .id = "clusterid"
)

# Concentric circles
r_lower = lapply(r_upper, function(x) {ifelse(x > 100, x - 100, x-50)})

plan(multisession)
permuted_concentric_t7 <- map_dfr(clusters, ~cluster_wrap_function(
  clusters = .x,
  df = perm_vcd,
  perms = perms,
  r1 = r_lower,
  r2 = r_upper,
  t = 7),
  .id = "clusterid"
)

permuted_concentric_t14 <- map_dfr(clusters, ~cluster_wrap_function(
  clusters = .x,
  df = perm_vcd,
  perms = perms,
  r1 = r_lower,
  r2 = r_upper,
  t = 14),
  .id = "clusterid"
)

output_cluster_permutations_temp <- list(permuted_radii_t7 = permuted_radii_t7,
                                    permuted_radii_t14 = permuted_radii_t14,
                                    permuted_concentric_t7 = permuted_concentric_t7,
                                    permuted_concentric_t14 = permuted_concentric_t14)

save(output_cluster_permutations_temp, 
     file = "cluster-specific-permutation-output-temporal-sensitivity.RData")
