library(knitr)
library(here)
library(tidyverse)
library(kableExtra)
library(ggpubr)
library(latex2exp)

load(here("data","cluster-specific-permutation-output.RData"))
load(here("data","cluster-specific-bootstrap-output.RData"))
load(here("data","2020-12-08_work-ts-spat.RData"))

# Intervention labels
tx_list <- work_ts_spat %>% 
  dplyr::select(cluster, intervention) %>% 
  mutate(intervention = ifelse(intervention == 0, "Untreated", "Intervention")) %>% 
  distinct() %>% 
  mutate(cluster = factor(cluster, levels = 1:24))

# Observed tau
observed_tau <- output_cluster_specific$observed_concentric %>% 
  dplyr::select(cluster = clusterid, 
                tau_obs = tau, 
                r_lower, r_upper)

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
  full_join(tx_list)

# Classic tau plots by cluster - UNTREATED
p1_u <- perm_c %>% 
  mutate(cluster = factor(cluster, levels = 1:24)) %>% 
  filter(intervention == "Untreated") %>% 
  ggplot(aes(x = r_upper, y = tau_obs)) + 
  facet_wrap(~cluster) + 
  geom_line() +
  geom_ribbon(aes(ymin = UL, ymax = LL),
              alpha = 0.2) + 
  xlab(TeX("$d_2")) + 
  ylab(TeX("Odds ratio ($\\tau)")) + 
  theme_pubr()

# Classic tau plots by cluster - TREATED
p1_t <- perm_c %>% 
  mutate(cluster = factor(cluster, levels = 1:24)) %>% 
  filter(intervention == "Intervention") %>% 
  ggplot(aes(x = r_upper, y = tau_obs)) + 
  facet_wrap(~cluster) + 
  geom_line() +
  geom_ribbon(aes(ymin = UL, ymax = LL),
              alpha = 0.2) + 
  xlab(TeX("$d_2")) + 
  ylab(TeX("Odds ratio ($\\tau)")) + 
  theme_pubr()


ggpubr::ggarrange(p1_u, p1_t, ncol = 1,
                  labels = "AUTO")
# ggsave(filename = here("graphs",
#                        paste0(Sys.Date(), "_cluster-specific-tau-perm-plots.png")),
#        device = "png",
#        width = 7, 
#        height = 7,
#        units = "in")

# Forest plots - UNTREATED
p2_u <- perm_c %>% 
  mutate(cluster = factor(cluster, levels = 1:24)) %>% 
  filter(intervention == "Untreated", r_upper <= 300) %>% 
  dplyr::select(-tau_perm) %>% 
  distinct() %>% 
  ggplot(aes(x = tau_obs, y = cluster)) + 
  facet_wrap(~r_upper) + 
  geom_errorbarh(aes(xmin = UL, xmax = LL),
              lwd = 0.75,
              col = "darkgray",
              height = 0.5) + 
  geom_vline(xintercept = 1,
             lty = 2, 
             col = "gray") +
  geom_point(#pch = 16,
             col = "darkgreen",
             size = 2
             ) +
  ylab("Cluster") + 
  xlab(TeX("Odds ratio ($\\tau)")) + 
  theme_pubr()


# Forest plots - TREATED
p2_t <- perm_c %>% 
  mutate(cluster = factor(cluster, levels = 1:24)) %>% 
  filter(intervention == "Intervention", r_upper <= 300) %>% 
  dplyr::select(-tau_perm) %>% 
  distinct() %>% 
  ggplot(aes(x = tau_obs, y = cluster)) + 
  facet_wrap(~r_upper) + 
  geom_errorbarh(aes(xmin = UL, xmax = LL),
                 lwd = 0.75,
                 col = "darkgray",
                 height = 0.5) + 
  geom_vline(xintercept = 1,
             lty = 2, 
             col = "gray") +
  geom_point(#pch = 16,
    col = "darkgreen",
    size = 2
  ) +
  ylab("Cluster") + 
  xlab(TeX("Odds ratio ($\\tau)")) + 
  theme_pubr()

ggarrange(p2_u, p2_t, ncol = 1,
          labels = "AUTO")
# ggsave(filename = here("graphs",
#                        paste0(Sys.Date(), "_cluster-specific-tau-perm-forest-plots.png")),
#        device = "png",
#        width = 7, 
#        height = 7,
#        units = "in")


# Tile plot
perm_c %>% 
  mutate(cluster = factor(cluster, levels = 1:24)) %>% 
 # filter(intervention == "Untreated") %>% 
  dplyr::select(-tau_perm) %>% 
  distinct() %>% 
  ggplot(aes(x = r_upper, y = cluster)) + 
  facet_wrap(~intervention, scales = "free") + 
  geom_tile(aes(fill = tau_obs)) +
  scale_fill_viridis_c(option = "C") +
  # ylab("Cluster") + 
  # xlab(TeX("Odds ratio ($\\tau)")) + 
  theme_pubr()

