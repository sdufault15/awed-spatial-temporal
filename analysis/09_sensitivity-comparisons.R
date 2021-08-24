library(tidyverse)
library(here)
library(ggpubr)
library(latex2exp)
load(here("data","cluster-specific-permutation-output.RData"))
load(here("data","cluster-specific-permutation-output-temporal-sensitivity.RData"))
load(here("data","cluster-specific-permutation-output-sensitivity-d50.RData"))
load(here("data","2020-12-08_work-ts-spat.RData"))
load(here("data","2021-01-08_distance-all-pts-to-boundary-sensitivity.RData"))

source(here("lib","binary-matrix-function_odds.R"))
source(here("lib","cluster-wrap-function.R"))
source(here("lib","bootstrap-function_odds.R"))
observed_data <- work_ts_spat %>% 
  filter(serotype != "unk_serotype")

tx_list <- work_ts_spat %>% 
  dplyr::select(clusterid = cluster, intervention) %>% 
  distinct() %>% 
  mutate(clusterid = factor(clusterid, levels = 1:24))


# Observed Tau
# Increasing Concentric
r_upper <- as.vector(c(50,seq(100, 1000, by = 100)), mode = "list")
r_lower <- lapply(r_upper, function(x) {ifelse(x > 100, x - 100, x-50)})
names(r_upper) <- r_upper
names(r_lower) <- r_lower

clusters <- as.list(1:24)
names(clusters) <- clusters

observed_concentric <- map_dfr(clusters, ~cluster_wrap_function(
  cluster = .x,
  df = observed_data,
  r1 = r_lower,
  r2 = r_upper,
  t = 30,
  replace = FALSE),
  .id = "clusterid"
)

d50_concentric <- map_dfr(clusters, ~cluster_wrap_function(
  cluster = .x,
  df = sensitivity_data$dist_all_s50,
  r1 = r_lower,
  r2 = r_upper,
  t = 30,
  replace = FALSE),
  .id = "clusterid"
  )

t7_concentric <- map_dfr(clusters, ~cluster_wrap_function(
  cluster = .x,
  df = observed_data,
  r1 = r_lower,
  r2 = r_upper,
  t = 7,
  replace = FALSE),
  .id = "clusterid"
)

t14_concentric <- map_dfr(clusters, ~cluster_wrap_function(
  cluster = .x,
  df = observed_data,
  r1 = r_lower,
  r2 = r_upper,
  t = 14,
  replace = FALSE),
  .id = "clusterid"
)

#######################
# Table
#######################

obs_mean_1 <- observed_concentric %>% 
  dplyr::select(clusterid, tau_obs = tau, r_lower, r_upper) %>% 
  full_join(output_cluster_permutations$permuted_concentric) %>% 
  dplyr::select(-a,-b,-c,-d) %>%
  dplyr::rename("tau_perm" = "tau") %>% 
  full_join(tx_list) %>% 
  mutate(tau_obs = ifelse(is.nan(tau_obs), 0, tau_obs)) %>% 
  mutate(tau_obs_1 = tau_obs + 1) %>% 
  group_by(intervention, r_upper) %>% 
  mutate(geom_mean_1 = exp(mean(log(tau_obs_1)))-1) %>% 
  rowwise() %>% 
  mutate(tau_perm = ifelse(is.nan(tau_perm), 0, tau_perm)) %>% 
  mutate(tau_perm_1 = tau_perm + 1) %>% 
  group_by(intervention, r_upper, permutation) %>% 
  mutate(geom_mean_perm_1 = exp(mean(log(tau_perm_1))) - 1) %>% 
  group_by(intervention, r_upper) %>% 
  mutate(CI.l = quantile(geom_mean_perm_1, probs = 0.025),
         CI.u = quantile(geom_mean_perm_1, probs = 0.975)) %>% 
  ungroup() %>% 
  dplyr::select(intervention, r_lower, r_upper, geom_mean_1, CI.l, CI.u) %>% 
  distinct() %>% 
  mutate(method = "obs")

d50_mean_1 <- d50_concentric %>% 
  dplyr::select(clusterid, tau_obs = tau, r_lower, r_upper) %>% 
  full_join(output_cluster_permutations_d50$permuted_concentric) %>% 
  dplyr::select(-a,-b,-c,-d) %>%
  dplyr::rename("tau_perm" = "tau") %>% 
  full_join(tx_list) %>% 
  mutate(tau_obs = ifelse(is.nan(tau_obs), 0, tau_obs)) %>% 
  mutate(tau_obs_1 = tau_obs + 1) %>% 
  group_by(intervention, r_upper) %>% 
  mutate(geom_mean_1 = exp(mean(log(tau_obs_1)))-1) %>% 
  rowwise() %>% 
  mutate(tau_perm = ifelse(is.nan(tau_perm), 0, tau_perm)) %>% 
  mutate(tau_perm_1 = tau_perm + 1) %>% 
  group_by(intervention, r_upper, permutation) %>% 
  mutate(geom_mean_perm_1 = exp(mean(log(tau_perm_1))) - 1) %>% 
  group_by(intervention, r_upper) %>% 
  mutate(CI.l = quantile(geom_mean_perm_1, probs = 0.025),
         CI.u = quantile(geom_mean_perm_1, probs = 0.975)) %>% 
  ungroup() %>% 
  dplyr::select(intervention, r_lower, r_upper, geom_mean_1, CI.l, CI.u) %>% 
  distinct() %>% 
  mutate(method = "d50")

t7_mean_1 <- t7_concentric %>% 
  dplyr::select(clusterid, tau_obs = tau, r_lower, r_upper) %>% 
  full_join(output_cluster_permutations_temp$permuted_concentric_t7) %>% 
  dplyr::select(-a,-b,-c,-d) %>%
  dplyr::rename("tau_perm" = "tau") %>% 
  full_join(tx_list) %>% 
  mutate(tau_obs = ifelse(is.nan(tau_obs), 0, tau_obs)) %>% 
  mutate(tau_obs_1 = tau_obs + 1) %>% 
  group_by(intervention, r_upper) %>% 
  mutate(geom_mean_1 = exp(mean(log(tau_obs_1)))-1) %>% 
  rowwise() %>% 
  mutate(tau_perm = ifelse(is.nan(tau_perm), 0, tau_perm)) %>% 
  mutate(tau_perm_1 = tau_perm + 1) %>% 
  group_by(intervention, r_upper, permutation) %>% 
  mutate(geom_mean_perm_1 = exp(mean(log(tau_perm_1))) - 1) %>% 
  group_by(intervention, r_upper) %>% 
  mutate(CI.l = quantile(geom_mean_perm_1, probs = 0.025),
         CI.u = quantile(geom_mean_perm_1, probs = 0.975)) %>% 
  ungroup() %>% 
  dplyr::select(intervention, r_lower, r_upper, geom_mean_1 = geom_mean_1, CI.l, CI.u) %>% 
  distinct() %>% 
  mutate(method = "t7")

t14_mean_1 <- t14_concentric %>% 
  dplyr::select(clusterid, tau_obs = tau, r_lower, r_upper) %>% 
  full_join(output_cluster_permutations_temp$permuted_concentric_t14) %>% 
  dplyr::select(-a,-b,-c,-d) %>%
  dplyr::rename("tau_perm" = "tau") %>% 
  full_join(tx_list) %>% 
  mutate(tau_obs = ifelse(is.nan(tau_obs), 0, tau_obs)) %>% 
  mutate(tau_obs_1 = tau_obs + 1) %>% 
  group_by(intervention, r_upper) %>% 
  mutate(geom_mean_1 = exp(mean(log(tau_obs_1)))-1) %>% 
  rowwise() %>% 
  mutate(tau_perm = ifelse(is.nan(tau_perm), 0, tau_perm)) %>% 
  mutate(tau_perm_1 = tau_perm + 1) %>% 
  group_by(intervention, r_upper, permutation) %>% 
  mutate(geom_mean_perm_1 = exp(mean(log(tau_perm_1))) - 1) %>% 
  group_by(intervention, r_upper) %>% 
  mutate(CI.l = quantile(geom_mean_perm_1, probs = 0.025),
         CI.u = quantile(geom_mean_perm_1, probs = 0.975)) %>% 
  ungroup() %>% 
  dplyr::select(intervention, r_lower, r_upper, geom_mean_1, CI.l, CI.u) %>% 
  distinct() %>% 
  mutate(method = "t14")


bind_rows(obs_mean_1, d50_mean_1, t7_mean_1, t14_mean_1) %>% 
  mutate(intervention = ifelse(intervention == 0, "Untreated", "Intervention"),
         method = case_when(method == "d50" ~ "Excluding 50m border",
                            method == "obs" ~ "Full data",
                            method == "t14" ~ "Within 2 weeks",
                            TRUE ~ "Within 1 week")) %>%
  mutate(method = factor(method, levels = c("Full data", 
                                            "Excluding 50m border", "Within 1 week", "Within 2 weeks"))) %>% 
  ggplot(aes(x = r_upper, y = geom_mean_1, linetype = method)) + 
  #scale_color_manual(values = c("#EF7822", "#76C1A8", "#003D58","#0AA4D1")) +
  scale_linetype_manual(values = 1:4) + 
  facet_wrap(~intervention, ncol = 1) + 
  geom_hline(yintercept = 1, 
             lty = 5,
             col = "darkgray") +
  geom_line(#size = 0.6
            ) +
  ylab(TeX("Odds ratio ($\\tau)")) +
  scale_x_continuous(TeX("Distance (m), $d_2"),
                     breaks = seq(100,1000,by = 100),
                     labels = seq(100,1000,by = 100)) +
  theme_pubr() + 
  theme(legend.position = "bottom",
        legend.title = element_blank())

ggsave(filename = here("graphs", 
                       paste0(Sys.Date(),
                              "_sensitivity-comparison-plot.png")),
       device = "png",
       height = 8,
       width = 8,
       units = "in")

# Table 

all_comp <- bind_rows(obs_mean_1, d50_mean_1, t7_mean_1, t14_mean_1) %>% 
  mutate(intervention = ifelse(intervention == 0, "Untreated", "Intervention"),
         method = case_when(method == "d50" ~ "Excluding 50m border",
                            method == "obs" ~ "Full data",
                            method == "t14" ~ "Within 2 weeks",
                            TRUE ~ "Within 1 week")) %>%
  mutate(method = factor(method, levels = c("Full data", 
                                            "Excluding 50m border", 
                                            "Within 1 week", 
                                            "Within 2 weeks"))) %>% 
  mutate(signif = ifelse(geom_mean_1 < CI.l,"below",
                         ifelse(geom_mean_1 > CI.u, "above", "null"))) 

all_comp %>% 
  transmute(intervention, method, signif,
            Interval = paste0(r_lower, " to ", r_upper, "m"),
            estimate = paste0(format(round(geom_mean_1, 2), nsmall = 2),
                                                        " (", format(round(CI.l, 2), nsmall = 2),
                                                        ", ", format(round(CI.u, 2), nsmall = 2),")")) %>% 
    mutate(estimate = cell_spec(estimate, color = ifelse(signif == "below", "red",
                                                         ifelse(signif == "above", "green", "black")))) %>% 
  pivot_wider(id_cols = c(intervention, Interval),
              names_from = method,
              values_from = estimate) %>% 
  kable(escape = F) %>% 
  kable_styling() %>% 
  add_header_above(c(" "=2, "Estimate (95% Null Region)"= 4)) %>% 
  collapse_rows(columns = 1)

## 
p1_int <- all_comp %>% 
  filter(intervention == "Intervention") %>% 
  ggplot(aes(x = r_upper, y = geom_mean_1)) +
  facet_wrap(~method,
             ncol = 4) +
  geom_hline(yintercept = 1,
             lty=2) +
  geom_line() +
  geom_ribbon(aes(ymin = CI.l, ymax = CI.u),
              alpha = 0.3) +
  coord_cartesian(ylim = c(0,5)) +
  ylab(TeX("Odds ratio ($\\tau)")) +
  scale_x_continuous(TeX("Distance (m), $d_2"),
                     breaks = seq(0,1000,by = 200),
                     labels = seq(0,1000,by = 200)) + 
  theme_pubr() +
  ggtitle("Intervention Arm") + 
  theme(axis.title.x = element_blank())

p1_unt <- all_comp %>% 
  filter(intervention == "Untreated") %>% 
  ggplot(aes(x = r_upper, y = geom_mean_1)) +
  facet_wrap(~method,
             ncol = 4) +
  geom_hline(yintercept = 1,
             lty=2) +
  geom_line() +
  geom_ribbon(aes(ymin = CI.l, ymax = CI.u),
              alpha = 0.3) +
  coord_cartesian(ylim = c(0,5)) +
  ylab(TeX("Odds ratio ($\\tau)")) +
  scale_x_continuous(TeX("Distance (m), $d_2"),
                     breaks = seq(0,1000,by = 200),
                     labels = seq(0,1000,by = 200)) + 
  theme_pubr() + 
  ggtitle("Untreated Arm")

ggarrange(p1_int,
          p1_unt,
          ncol = 1,
          align = "v")
ggsave(filename = here("graphs",
                       paste0(Sys.Date(),
                              "_sensitivity-comparison-line-plot.png")),
       device = "png",
       width = 8.5,
       height = 6,
       units = "in")
