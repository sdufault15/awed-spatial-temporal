library(knitr)
library(here)
library(tidyverse)
library(kableExtra)
library(ggpubr)
library(latex2exp)
library(geodist)

load(here("data","cluster-specific-permutation-output.RData"))
load(here("data","2020-12-08_work-ts-spat.RData"))
source(here("lib","cluster-wrap-function.R"))
source(here("lib","bootstrap-function_odds.R"))
source(here("lib","binary-matrix-function_odds.R"))

# Intervention labels
tx_list <- work_ts_spat %>% 
  dplyr::select(cluster, intervention) %>% 
  mutate(intervention = ifelse(intervention == 0, "Untreated", "Intervention")) %>% 
  distinct() %>% 
  mutate(cluster = factor(cluster, levels = 1:24))

##############
# Observed tau
##############
observed_data <- work_ts_spat %>%
  filter(serotype != "unk_serotype")

# Increasing Radius
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

observed_tau <- observed_concentric %>% 
  dplyr::select(cluster = clusterid, r_lower, r_upper, tau_obs = tau)

# Permutation data
perm_c <- output_cluster_permutations$permuted_concentric %>% 
  dplyr::select(cluster = clusterid,
                tau_perm = tau,
                r_lower, r_upper)

# Merging for graphing purposes
perm_c <- full_join(observed_tau, perm_c)

# Finding 95% null distribution
perm_c <- perm_c %>% 
  group_by(cluster, r_upper) %>% 
  mutate(UL = quantile(tau_perm, probs = 0.975, na.rm = TRUE),
         LL = quantile(tau_perm, probs = 0.025, na.rm = TRUE)) %>% 
  full_join(tx_list) %>% 
  dplyr::select(intervention, cluster, r_lower, r_upper, tau_obs, UL, LL) %>% 
  distinct()  

# Geometric means
plus_one <- observed_tau %>% 
  full_join(tx_list) %>% 
  group_by(intervention, r_upper) %>% 
  mutate(tau_obs = ifelse(is.nan(tau_obs), 0, tau_obs)) %>% 
  mutate(tau_1 = tau_obs + 1) %>% 
  summarise(geom_mean_1 = exp(mean(log(tau_1))) - 1) 

null_bands_perm_1 <- output_cluster_permutations$permuted_concentric %>% 
  full_join(tx_list, by = c("clusterid" = "cluster")) %>% 
  mutate(tau_perm = ifelse(is.nan(tau), 0, tau)) %>% 
  mutate(tau_perm_1 = tau_perm + 1) %>% 
  group_by(intervention, r_upper, permutation) %>% 
  summarise(geom_mean_perm_1 = exp(mean(log(tau_perm_1))) - 1) %>% 
  group_by(intervention, r_upper) %>% 
  mutate(LL = quantile(geom_mean_perm_1, probs = 0.025),
         UL = quantile(geom_mean_perm_1, probs = 0.975)) %>% 
  dplyr::select(intervention, r_upper, UL, LL) %>% 
  distinct()

##############
## Untreated Plot
##############3
fp_obs <- perm_c %>% 
  full_join(tx_list) %>% 
  filter(r_upper <= 300) %>% 
  dplyr::select(cluster, tau_obs)

fp_perm_mean <- plus_one %>% 
  full_join(null_bands_perm_1) %>% 
  filter(r_upper <= 300) 

p_unt <- perm_c %>% 
  dplyr::select(intervention, cluster, r_lower, r_upper, tau_obs, LL, UL) %>% 
  distinct() %>% 
  mutate(cluster = factor(cluster, levels = 1:24),
         colp = ifelse(tau_obs > UL, "red", "black"),
         r_upper = ifelse(r_upper > 100,
                          paste0(r_upper - 100, " to ", r_upper,"m"),
                          paste0(r_upper - 50, " to ", r_upper,"m"))) %>% 
  filter(as.numeric(r_lower) < 300, intervention == "Untreated") %>% 
  mutate(r_upper = factor(r_upper, levels = c("0 to 50m",
                                              "50 to 100m",
                                              "100 to 200m",
                                              "200 to 300m"))) %>% 
  ggplot(aes(x = tau_obs, y = cluster)) + 
  facet_wrap(~r_upper, scales = "free_y", nrow = 1) + 
  geom_errorbarh(aes(xmin = LL, xmax = UL),
                 height = 0.5,
                 size = 0.75,
                 col = "darkgray") + 
  geom_point(aes(col = colp)) + 
  scale_color_manual(values = c("#003D58","#0AA4D1")) +
  xlab(TeX("Odds ratio ($\\tau)")) +
  ylab(TeX("Cluster")) + 
  ggtitle("Untreated Arm") + 
  coord_cartesian(xlim = c(0,15)) + 
  geom_vline(xintercept = 1, 
             lty = 2) +
  #  scale_x_log10() +
  theme_pubr()+
  theme(legend.position = "none")

overall_unt <- fp_perm_mean %>% 
  filter(intervention == "Untreated") %>% 
  mutate(colorp = ifelse(geom_mean_1 > UL, "red", "black"),
         r_upper = ifelse(r_upper > 100,
                          paste0(r_upper - 100, " to ", r_upper,"m"),
                          paste0(r_upper - 50, " to ", r_upper,"m"))) %>% 
  mutate(r_upper = factor(r_upper, levels = c("0 to 50m",
                                              "50 to 100m",
                                              "100 to 200m",
                                              "200 to 300m"))) %>% 
  ggplot(aes(x = geom_mean_1, y = 1)) + 
  facet_wrap(~r_upper,nrow = 1) +
  geom_errorbarh(aes(xmin = LL, xmax = UL),
                 height = 0.5,
                 size = 0.75,
                 col = "darkgray") + 
  geom_point(aes(color = colorp),
             pch = 18,
             size = 3) +
  scale_color_manual(values = c("#0AA4D1","#003D58")) +
  ylab(TeX("Cluster")) + 
  geom_vline(xintercept = 1, 
             lty = 2) +
  scale_x_continuous(name = TeX("Odds ratio ($\\tau)")) +
  scale_y_continuous("", breaks = 1,
                     labels = "Overall") + 
  coord_cartesian(xlim = c(0,15),
                  ylim = c(0.5,1.5)) + 
  theme_pubr() +
  theme(axis.title.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.line.y = element_blank(),
        strip.background = element_blank(),
        strip.text.x = element_blank(),
        legend.position = "none")

p_unt_2 <- 
  p_unt +
  theme(axis.title.x = element_blank())

all_unt <- ggarrange(p_unt_2,
                     overall_unt,
                     ncol = 1,
                     align = "v",
                     heights = c(4,1))

##############
## Intervention Plot
##############
p_int <- perm_c %>% 
  dplyr::select(intervention, cluster, r_lower, r_upper, tau_obs, LL, UL) %>% 
  distinct() %>% 
  mutate(cluster = factor(cluster, levels = 1:24),
         colp = ifelse(tau_obs > UL, "red", "black"),
         r_upper = ifelse(r_upper > 100,
                          paste0(r_upper - 100, " to ", r_upper,"m"),
                          paste0(r_upper - 50, " to ", r_upper,"m"))) %>%
  filter(as.numeric(r_lower) < 300, intervention == "Intervention") %>% 
  mutate(r_upper = factor(r_upper, levels = c("0 to 50m",
                                              "50 to 100m",
                                              "100 to 200m",
                                              "200 to 300m"))) %>% 
  ggplot(aes(x = tau_obs, y = cluster)) + 
  facet_wrap(~r_upper, scales = "free_y", nrow = 1) + 
  geom_errorbarh(aes(xmin = LL, xmax = UL),
                 height = 0.5,
                 size = 0.75,
                 col = "darkgray") + 
  geom_point(aes(col = colp)) + 
  scale_color_manual(values = c("#003D58","#0AA4D1")) +
  xlab(TeX("Odds ratio ($\\tau)")) +
  ylab(TeX("Cluster")) + 
  coord_cartesian(xlim = c(0,27)) +
  geom_vline(xintercept = 1, 
             lty = 2) +
  #scale_x_log10() +
  ggtitle("Intervention Arm") + 
  theme_pubr()+
  theme(legend.position = "none")


overall_int <- fp_perm_mean %>% 
  filter(intervention == "Intervention") %>% 
  mutate(colorp = ifelse(geom_mean_1 > UL, "black", "red"),
         r_upper = ifelse(r_upper > 100,
                          paste0(r_upper - 100, " to ", r_upper,"m"),
                          paste0(r_upper - 50, " to ", r_upper,"m"))) %>% 
  mutate(r_upper = factor(r_upper, levels = c("0 to 50m",
                                              "50 to 100m",
                                              "100 to 200m",
                                              "200 to 300m"))) %>% 
  ggplot(aes(x = geom_mean_1, y = 1)) + 
  facet_wrap(~r_upper, nrow = 1) +
  geom_errorbarh(aes(xmin = LL, xmax = UL),
                 height = 0.5,
                 size = 0.75,
                 col = "darkgray") + 
  geom_point(aes(color = colorp),
             pch = 18,
             size = 3) +
  scale_color_manual(values = c("#003D58","#0AA4D1")) +
  ylab(TeX("Cluster")) + 
  geom_vline(xintercept = 1, 
             lty = 2) +
  scale_x_continuous(name = TeX("Odds ratio ($\\tau)"),
                     breaks = seq(0,30,by = 10),
                     labels = seq(0,30,by = 10)) +
  scale_y_continuous("", breaks = 1,
                     labels = "Overall") + 
  coord_cartesian(xlim = c(0,27),
                  ylim = c(0.5,1.5)) + 
  theme_pubr() +
  theme(axis.title.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.line.y = element_blank(),
        strip.background = element_blank(),
        strip.text.x = element_blank(),
        legend.position = "none")

p_int_2 <- 
  p_int +
  theme(axis.title.x = element_blank())

all_int <- ggarrange(p_int_2,
                     overall_int,
                     ncol = 1,
                     align = "v",
                     heights = c(4,1))

ggarrange(all_int,
          all_unt,
          ncol = 1)
ggsave(filename = here("graphs", 
                       paste0(Sys.Date(), "_forestplot-plus-one.png")),
       device = "png",
       height = 8,
       width = 9,
       units = "in")

### Generating processed data

fp_perm_mean %>% 
  rename(`Geom Mean LL` = LL, `Geom Mean UL` = UL) %>% 
  full_join(perm_c) %>% 
  mutate(cluster = factor(cluster, levels = 1:24)) %>% 
  arrange(intervention, r_upper, cluster) %>% 
  dplyr::select(intervention, r_lower, r_upper, cluster, tau_obs, 
                UL, LL, geom_mean_1, `Geom Mean LL`,`Geom Mean UL`) %>% 
  write.csv(file = here("graphs", paste0("processed-data/", Sys.Date(), "_fig4-data.csv")),
            row.names = FALSE)
  
  
