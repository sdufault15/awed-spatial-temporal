# Sent to the cluster for efficiency

library(tidyverse)
library(geodist)
library(furrr)
library(purrr)

load("2020-12-08_work-ts-spat.RData") # Data

########################
# Analysis
########################
source(here("lib", "bootstrap-function_odds.R"))
source(here("lib", "binary-matrix-function_odds.R"))

bootstrap_vcd <- work_ts_spat %>% 
  filter(serotype != "unk_serotype")


# Increasing Radius
r_upper <- as.vector(c(50,seq(100, 1500, by = 100)), mode = "list")
r_lower <- as.vector(rep(0,length(r_upper)), mode = "list")
names(r_upper) <- r_upper
names(r_lower) <- r_lower
bs <- as.vector(1:1000, mode = "list")

observed <- bootstrap_function_odds(
  df = bootstrap_vcd, 
  r1 = r_lower, 
  r2 = r_upper,
  t = 30,
  # Set replace to FALSE for observed estimates
  replace = FALSE) 

plan(multisession)
bootstrapped_radii <- future_map_dfr(bs, ~bootstap_function_odds(df = bootstrap_vcd,
                                                                 r1 = r_lower,
                                                                 r2 = r_upper,
                                                                 t = 30,
                                                                 replace = TRUE),
                                     .id = "bootstrap")

# Concentric circles
r_lower = lapply(r_upper, function(x) {ifelse(x > 100, x - 100, x-50)})

observed_concentric <- bootstrap_function_odds(
  df = bootstrap_vcd,
  r1 = r_lower,
  r2 = r_upper,
  t = 30,
  # Set replace to FALSE for observed estimates
  replace = FALSE)

plan(multisession)
bootstrapped_concentric <- future_map_dfr(bs, ~bootstrap_function_odds(df = bootstrap_vcd,
                                                                      r1 = r_lower,
                                                                      r2 = r_upper,
                                                                      t = 30,
                                                                      replace = TRUE),
                                          .id = "bootstrap")

output <- list(observed_radii = observed,
               observed_concentric = observed_concentric,
               bootstrapped_radii = bootstrapped_radii,
               bootstrapped_concentric = bootstrapped_concentric)

save(output, 
     file = "overall-bootstrap-output.RData")
