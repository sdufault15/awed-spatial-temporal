library(tidyverse)
library(here)
library(ggpubr)
library(latex2exp)

load(here("data","cluster-specific-permutation-output-sensitivity-d50.RData"))
load(here("data","2020-12-08_work-ts-spat.RData"))

tx <- work_ts_spat %>% 
  select(cluster, intervention) %>% 
  distinct() %>% 
  mutate(intervention = ifelse(intervention == 0, "Untreated", "Intervention"))

##################
# Observed tau
##################
source(here("lib","binary-matrix-function_odds.R"))
source(here("lib","cluster-wrap-function.R"))
source(here("lib","bootstrap-function_odds.R"))
observed_data <- work_ts_spat %>% 
  filter(serotype != "unk_serotype")

# Increasing Radius
r_upper <- as.vector(c(50,seq(100, 1000, by = 100)), mode = "list")
r_lower <- as.vector(rep(0,length(r_upper)), mode = "list")
names(r_upper) <- r_upper
names(r_lower) <- r_lower

clusters <- as.list(1:24)
names(clusters) <- clusters
observed_radii <- map_dfr(clusters, ~cluster_wrap_function(
  cluster = .x,
  df = observed_data,
  r1 = r_lower,
  r2 = r_upper,
  t = 30,
  replace = FALSE),
  .id = "clusterid"
)

r_lower <- lapply(r_upper, function(x) {ifelse(x > 100, x - 100, x-50)})
observed_concentric <- map_dfr(clusters, ~cluster_wrap_function(
  cluster = .x,
  df = observed_data,
  r1 = r_lower,
  r2 = r_upper,
  t = 30,
  replace = FALSE),
  .id = "clusterid"
)
##################

##################
## Setting up the data
# Increasing radii
df_radii <- observed_radii %>% 
  full_join(mutate(tx, cluster = paste0(cluster)), by = c("clusterid" = "cluster"))
cis_r <- output_cluster_permutations_d50$permuted_radii %>% 
  group_by(clusterid, r_lower, r_upper) %>% 
  summarise(CI.l = quantile(tau, probs = 0.025, na.rm = TRUE),
            CI.u = quantile(tau, probs = 0.975, na.rm = TRUE))
df_radii <- full_join(df_radii,
                      cis_r, by = c("clusterid", "r_upper"))

plus_one_r <- observed_radii %>% 
  mutate(clusterid = as.numeric(clusterid)) %>% 
  full_join(distinct(dplyr::select(observed_data, cluster, intervention)), by = c("clusterid" = "cluster")) %>% 
  group_by(intervention, r_upper) %>% 
  mutate(tau_obs = ifelse(is.nan(tau), 0, tau)) %>% 
  mutate(tau_1 = tau_obs + 1) %>% 
  summarise(geom_mean_1 = exp(mean(log(tau_1))) - 1) 

plus_one_r_ci <- output_cluster_permutations_d50$permuted_radii %>% 
  mutate(clusterid = as.numeric(clusterid)) %>% 
  full_join(distinct(dplyr::select(observed_data, cluster, intervention)), by = c("clusterid" = "cluster")) %>% 
  mutate(tau_perm = ifelse(is.nan(tau), 0, tau)) %>% 
  mutate(tau_perm_1 = tau_perm + 1) %>% 
  group_by(intervention, r_upper, permutation) %>% 
  summarise(geom_mean_perm_1 = exp(mean(log(tau_perm_1))) - 1) %>% 
  group_by(intervention, r_upper) %>% 
  mutate(CI.l = quantile(geom_mean_perm_1, probs = 0.025),
         CI.u = quantile(geom_mean_perm_1, probs = 0.975)) %>% 
  dplyr::select(intervention, r_upper, CI.l, CI.u) %>% 
  distinct()

plus_one_r <- full_join(plus_one_r, plus_one_r_ci)

# Concentric rings
df_concentric <- observed_concentric %>% 
  full_join(mutate(tx, cluster = paste0(cluster)), by = c("clusterid" = "cluster"))
cis_c <- output_cluster_permutations_d50$permuted_concentric %>% 
  group_by(clusterid, r_upper) %>% 
  summarise(CI.l = quantile(tau, probs = 0.025, na.rm = TRUE),
            CI.u = quantile(tau, probs = 0.975, na.rm = TRUE))
df_concentric <- full_join(df_concentric, 
                           cis_c, by = c("clusterid", "r_upper"))

plus_one_c <- observed_concentric %>% 
  mutate(clusterid = as.numeric(clusterid)) %>% 
  full_join(distinct(dplyr::select(observed_data, cluster, intervention)), by = c("clusterid" = "cluster")) %>% 
  group_by(intervention, r_upper) %>% 
  mutate(tau_obs = ifelse(is.nan(tau), 0, tau)) %>% 
  mutate(tau_1 = tau_obs + 1) %>% 
  summarise(geom_mean_1 = exp(mean(log(tau_1))) - 1) 

plus_one_c_ci <- output_cluster_permutations_d50$permuted_concentric %>% 
  mutate(clusterid = as.numeric(clusterid)) %>% 
  full_join(distinct(dplyr::select(observed_data, cluster, intervention)), by = c("clusterid" = "cluster")) %>% 
  mutate(tau_perm = ifelse(is.nan(tau), 0, tau)) %>% 
  mutate(tau_perm_1 = tau_perm + 1) %>% 
  group_by(intervention, r_upper, permutation) %>% 
  summarise(geom_mean_perm_1 = exp(mean(log(tau_perm_1))) - 1) %>% 
  group_by(intervention, r_upper) %>% 
  mutate(CI.l = quantile(geom_mean_perm_1, probs = 0.025),
         CI.u = quantile(geom_mean_perm_1, probs = 0.975)) %>% 
  dplyr::select(intervention, r_upper, CI.l, CI.u) %>% 
  distinct()

plus_one_c <- full_join(plus_one_c, plus_one_c_ci)

##################

## Forest Plots - Concentric
p_unt_forest_concentric <- df_concentric %>% 
  filter(r_upper <= 300, intervention == "Untreated") %>% 
  mutate(cluster = factor(clusterid, levels = 1:24),
         colorp = ifelse(tau > CI.u, "black", "red"),
         r_upper = ifelse(r_upper > 100, paste0(r_upper - 100, " to ", r_upper,"m"), paste0(r_upper - 50, " to ", r_upper,"m")))  %>% 
  mutate(r_upper = factor(r_upper, 
                          levels = c("0 to 50m",
                                     "50 to 100m",
                                     "100 to 200m",
                                     "200 to 300m"))) %>% 
  ggplot(aes(x = tau, y = cluster)) + 
  facet_wrap(~r_upper, nrow = 1) +
  geom_vline(xintercept = 1,
             lty = 2) + 
  geom_errorbarh(aes(xmin = CI.l, xmax = CI.u),
                 height = 0.5,
                 size = 0.75,
                 col = "darkgray") +
  geom_point(aes(color = colorp)) +
  scale_color_manual(values = c("#0AA4D1","#003D58")) +
  ylab(TeX("Cluster")) + 
  geom_vline(xintercept = 1, 
             lty = 2) +
  scale_x_continuous(name = TeX("Odds ratio ($\\tau)")) +
  theme_pubr() +
  theme(legend.position = "none",
        axis.title.x = element_blank()) +
  ggtitle("Untreated Arm")

p_unt_forest_concentric_overall <- plus_one_c %>% 
  filter(r_upper <= 300, intervention == 0) %>% 
  mutate(colorp = ifelse(geom_mean_1 > CI.u, "black", "red"),
         r_upper = ifelse(r_upper > 100, paste0(r_upper - 100, " to ", r_upper,"m"), paste0(r_upper - 50, " to ", r_upper,"m")))  %>% 
  mutate(r_upper = factor(r_upper, 
                          levels = c("0 to 50m",
                                     "50 to 100m",
                                     "100 to 200m",
                                     "200 to 300m"))) %>% 
  ggplot(aes(x = geom_mean_1, y = 1)) + 
  facet_wrap(~r_upper, nrow = 1) +
  geom_vline(xintercept = 1,
             lty = 2) + 
  geom_errorbarh(aes(xmin = CI.l, xmax = CI.u),
                 height = 0.5,
                 size = 0.75,
                 col = "darkgray") +
  geom_point(aes(color = colorp)) +
  scale_color_manual(values = c("#0AA4D1","#003D58")) +
  scale_y_continuous(name = "",
                     breaks = 1,
                     labels = "Overall") +
  scale_x_continuous(name = TeX("Odds ratio ($\\tau)")) +
  theme_pubr() +
  coord_cartesian(xlim = c(0,15),
                  ylim = c(0.5,1.5)) + 
  theme(axis.title.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.line.y = element_blank(),
        strip.background = element_blank(),
        strip.text.x = element_blank(),
        legend.position = "none") 

p_unt_all_c <- ggarrange(p_unt_forest_concentric, p_unt_forest_concentric_overall,
          align = "v",
          heights = c(4,1),
          ncol = 1)

p_int_forest_concentric <- df_concentric %>% 
  filter(r_upper <= 300, intervention == "Intervention") %>% 
  mutate(cluster = factor(clusterid, levels = 1:24),
         colorp = ifelse(tau > CI.u, "black", "red"),
         r_upper = ifelse(r_upper > 100, paste0(r_upper - 100, " to ", r_upper,"m"), paste0(r_upper - 50, " to ", r_upper,"m")))  %>% 
  mutate(r_upper = factor(r_upper, 
                          levels = c("0 to 50m",
                                     "50 to 100m",
                                     "100 to 200m",
                                     "200 to 300m"))) %>% 
  ggplot(aes(x = tau, y = cluster)) + 
  facet_wrap(~r_upper, nrow = 1) +
  geom_vline(xintercept = 1,
             lty = 2) + 
  geom_errorbarh(aes(xmin = CI.l, xmax = CI.u),
                 height = 0.5,
                 size = 0.75,
                 col = "darkgray") +
  geom_point(aes(color = colorp)) +
  scale_color_manual(values = c("#003D58" ,"#0AA4D1")) +
  ylab(TeX("Cluster")) + 
  geom_vline(xintercept = 1, 
             lty = 2) +
  scale_x_continuous(name = TeX("Odds ratio ($\\tau)")) +
  theme_pubr() +
  theme(legend.position = "none",
        axis.title.x = element_blank()) +
  ggtitle("Intervention Arm")

p_int_forest_concentric_overall <- plus_one_c %>% 
  filter(r_upper <= 300, intervention == 1) %>% 
  mutate(colorp = ifelse(geom_mean_1 > CI.u, "black", "red"),
         r_upper = ifelse(r_upper > 100, paste0(r_upper - 100, " to ", r_upper,"m"), paste0(r_upper - 50, " to ", r_upper,"m")))  %>% 
  mutate(r_upper = factor(r_upper, 
                          levels = c("0 to 50m",
                                     "50 to 100m",
                                     "100 to 200m",
                                     "200 to 300m"))) %>% 
  ggplot(aes(x = geom_mean_1, y = 1)) + 
  facet_wrap(~r_upper, nrow = 1) +
  geom_vline(xintercept = 1,
             lty = 2) + 
  geom_errorbarh(aes(xmin = CI.l, xmax = CI.u),
                 height = 0.5,
                 size = 0.75,
                 col = "darkgray") +
  geom_point(aes(color = colorp)) +
  scale_color_manual(values = c("#003D58", "#0AA4D1")) +
  scale_y_continuous(name = "",
                     breaks = 1,
                     labels = "Overall") +
  scale_x_continuous(name = TeX("Odds ratio ($\\tau)")) +
  theme_pubr() +
  coord_cartesian(xlim = c(0,30),
                  ylim = c(0.5,1.5)) + 
  theme(axis.title.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.line.y = element_blank(),
        strip.background = element_blank(),
        strip.text.x = element_blank(),
        legend.position = "none") 

p_int_all_c <- ggarrange(p_int_forest_concentric, p_int_forest_concentric_overall,
                         align = "v",
                         heights = c(4,1),
                         ncol = 1)


ggarrange(p_int_all_c,
          p_unt_all_c,
          labels = "AUTO",
          ncol = 1)
ggsave(filename = here("graphs",
                       paste0(Sys.Date(), "_forest-plot-spatial-sensitivity-concentric.png")),
       device = "png",
       height = 8,
       width = 9,
       units = "in")

## Forest Plots - Radii
p_unt_forest_radii <- df_radii %>% 
  filter(r_upper <= 300, intervention == "Untreated") %>% 
  mutate(cluster = factor(clusterid, levels = 1:24),
         colorp = ifelse(tau > CI.u, "black", "red"),
         r_upper = ifelse(r_upper > 100, paste0(r_upper - 100, " to ", r_upper,"m"), paste0(r_upper - 50, " to ", r_upper,"m")))  %>% 
  mutate(r_upper = factor(r_upper, 
                          levels = c("0 to 50m",
                                     "50 to 100m",
                                     "100 to 200m",
                                     "200 to 300m"))) %>% 
  ggplot(aes(x = tau, y = cluster)) + 
  facet_wrap(~r_upper, nrow = 1) +
  geom_vline(xintercept = 1,
             lty = 2) + 
  geom_errorbarh(aes(xmin = CI.l, xmax = CI.u),
                 height = 0.5,
                 size = 0.75,
                 col = "darkgray") +
  geom_point(aes(color = colorp)) +
  scale_color_manual(values = c("#0AA4D1","#003D58")) +
  ylab(TeX("Cluster")) + 
  geom_vline(xintercept = 1, 
             lty = 2) +
  scale_x_continuous(name = TeX("Odds ratio ($\\tau)")) +
  theme_pubr() +
  theme(legend.position = "none") +
  ggtitle("Untreated Arm")

p_int_forest_radii <- df_radii %>% 
  filter(r_upper <= 300, intervention == "Intervention") %>% 
  mutate(cluster = factor(clusterid, levels = 1:24),
         colorp = ifelse(tau > CI.u, "black", "red"),
         r_upper = ifelse(r_upper > 100, paste0(r_upper - 100, " to ", r_upper,"m"), paste0(r_upper - 50, " to ", r_upper,"m")))  %>% 
  mutate(r_upper = factor(r_upper, 
                          levels = c("0 to 50m",
                                     "50 to 100m",
                                     "100 to 200m",
                                     "200 to 300m"))) %>% 
  ggplot(aes(x = tau, y = cluster)) + 
  facet_wrap(~r_upper, nrow = 1) +
  geom_vline(xintercept = 1,
             lty = 2) + 
  geom_errorbarh(aes(xmin = CI.l, xmax = CI.u),
                 height = 0.5,
                 size = 0.75,
                 col = "darkgray") +
  geom_point(aes(color = colorp)) +
  scale_color_manual(values = c("#0AA4D1","#003D58")) +
  ylab(TeX("Cluster")) + 
  geom_vline(xintercept = 1, 
             lty = 2) +
  scale_x_continuous(name = TeX("Odds ratio ($\\tau)")) +
  theme_pubr() +
  theme(legend.position = "none") +
  ggtitle("Intervention Arm")

p2 <- df_concentric %>% 
  ggplot(aes(x = r_upper, y = tau, col = clusterid)) + 
  facet_wrap(~intervention, ncol = 1) +
  geom_hline(yintercept = 1,
             lty = 2) + 
  geom_line() +
  geom_smooth(aes(group = 1)) + 
  theme_pubr() + 
  scale_x_continuous("Distance range (m)",
                     breaks = seq(0,1000,by = 200),
                     labels = seq(0,1000,by = 200)) + 
  scale_y_continuous(TeX("$\\tau_{cluster} = \\frac{\\hat{\\theta}(d_1,d_2)}{\\hat{\\theta}(0,\\infty)}")) + 
  theme(legend.position = "none")

ggpubr::ggarrange(p1, p2,
                  ncol = 2,
                  labels = "AUTO",
                  align = "hv")
# ggsave(filename = here("graphs", 
#                        paste0(Sys.Date(), "_cluster-specific-tau-dual-plot-sensitivity50m.png")),
#        height = 5,
#        width = 8,
#        units = "in",
#        device = "png")


###### T = 5
load(here("data","cluster-specific-bootstrap-output-sensitivity-t5.RData"))

## Setting up the data
# Increasing radii
df_radii <- output_cluster_specific_t5$observed_radii %>% 
  full_join(mutate(tx, cluster = paste0(cluster)), by = c("clusterid" = "cluster"))
cis_r <- output_cluster_specific_t5$bootstrapped_radii %>% 
  group_by(clusterid, r_upper) %>% 
  summarise(CI.l = quantile(tau, probs = 0.025, na.rm = TRUE),
            CI.u = quantile(tau, probs = 0.975, na.rm = TRUE))
df_radii <- full_join(df_radii,
                      cis_r, by = c("clusterid", "r_upper"))
# Concentric rings
df_concentric <- output_cluster_specific_t5$observed_concentric %>% 
  full_join(mutate(tx, cluster = paste0(cluster)), by = c("clusterid" = "cluster"))
cis_c <- output_cluster_specific_t5$bootstrapped_concentric %>% 
  group_by(clusterid, r_upper) %>% 
  summarise(CI.l = quantile(tau, probs = 0.025, na.rm = TRUE),
            CI.u = quantile(tau, probs = 0.975, na.rm = TRUE))
df_concentric <- full_join(df_concentric, 
                           cis_c, by = c("clusterid", "r_upper"))

## Overall 
p1 <- df_radii %>% 
  ggplot(aes(x = r_upper, y = tau, col = clusterid)) + 
  facet_wrap(~intervention, ncol = 1) +
  geom_hline(yintercept = 1,
             lty = 2) + 
  geom_line() +
  geom_smooth(aes(group = 1)) + 
  theme_pubr() + 
  scale_x_continuous("Distance range (m)",
                     breaks = seq(0,1000,by = 200),
                     labels = seq(0,1000,by = 200)) + 
  scale_y_continuous(TeX("$\\tau_{cluster} = \\frac{\\hat{\\theta}(0,d_2)}{\\hat{\\theta}(0,\\infty)}")) + 
  theme(legend.position = "none")

p2 <- df_concentric %>% 
  ggplot(aes(x = r_upper, y = tau, col = clusterid)) + 
  facet_wrap(~intervention, ncol = 1) +
  geom_hline(yintercept = 1,
             lty = 2) + 
  geom_line() +
  geom_smooth(aes(group = 1)) + 
  theme_pubr() + 
  scale_x_continuous("Distance range (m)",
                     breaks = seq(0,1000,by = 200),
                     labels = seq(0,1000,by = 200)) + 
  scale_y_continuous(TeX("$\\tau_{cluster} = \\frac{\\hat{\\theta}(d_1,d_2)}{\\hat{\\theta}(0,\\infty)}")) + 
  theme(legend.position = "none")

ggpubr::ggarrange(p1, p2,
                  ncol = 2,
                  labels = "AUTO",
                  align = "hv")
ggsave(filename = here("graphs", 
                       paste0(Sys.Date(), "_cluster-specific-tau-dual-plot-sensitivity-t5.png")),
       height = 5,
       width = 8,
       units = "in",
       device = "png")



###### T = 15
load(here("data","cluster-specific-bootstrap-output-sensitivity-t15.RData"))

## Setting up the data
# Increasing radii
df_radii <- output_cluster_specific_t15$observed_radii %>% 
  full_join(mutate(tx, cluster = paste0(cluster)), by = c("clusterid" = "cluster"))
cis_r <- output_cluster_specific_t15$bootstrapped_radii %>% 
  group_by(clusterid, r_upper) %>% 
  summarise(CI.l = quantile(tau, probs = 0.025, na.rm = TRUE),
            CI.u = quantile(tau, probs = 0.975, na.rm = TRUE))
df_radii <- full_join(df_radii,
                      cis_r, by = c("clusterid", "r_upper"))
# Concentric rings
df_concentric <- output_cluster_specific_t15$observed_concentric %>% 
  full_join(mutate(tx, cluster = paste0(cluster)), by = c("clusterid" = "cluster"))
cis_c <- output_cluster_specific_t15$bootstrapped_concentric %>% 
  group_by(clusterid, r_upper) %>% 
  summarise(CI.l = quantile(tau, probs = 0.025, na.rm = TRUE),
            CI.u = quantile(tau, probs = 0.975, na.rm = TRUE))
df_concentric <- full_join(df_concentric, 
                           cis_c, by = c("clusterid", "r_upper"))

## Overall 
p1 <- df_radii %>% 
  ggplot(aes(x = r_upper, y = tau, col = clusterid)) + 
  facet_wrap(~intervention, ncol = 1) +
  geom_hline(yintercept = 1,
             lty = 2) + 
  geom_line() +
  geom_smooth(aes(group = 1)) + 
  theme_pubr() + 
  scale_x_continuous("Distance range (m)",
                     breaks = seq(0,1000,by = 200),
                     labels = seq(0,1000,by = 200)) + 
  scale_y_continuous(TeX("$\\tau_{cluster} = \\frac{\\hat{\\theta}(0,d_2)}{\\hat{\\theta}(0,\\infty)}")) + 
  theme(legend.position = "none")

p2 <- df_concentric %>% 
  ggplot(aes(x = r_upper, y = tau, col = clusterid)) + 
  facet_wrap(~intervention, ncol = 1) +
  geom_hline(yintercept = 1,
             lty = 2) + 
  geom_line() +
  geom_smooth(aes(group = 1)) + 
  theme_pubr() + 
  scale_x_continuous("Distance range (m)",
                     breaks = seq(0,1000,by = 200),
                     labels = seq(0,1000,by = 200)) + 
  scale_y_continuous(TeX("$\\tau_{cluster} = \\frac{\\hat{\\theta}(d_1,d_2)}{\\hat{\\theta}(0,\\infty)}")) + 
  theme(legend.position = "none")

ggpubr::ggarrange(p1, p2,
                  ncol = 2,
                  labels = "AUTO",
                  align = "hv")
ggsave(filename = here("graphs", 
                       paste0(Sys.Date(), "_cluster-specific-tau-dual-plot-sensitivity-t15.png")),
       height = 5,
       width = 8,
       units = "in",
       device = "png")
