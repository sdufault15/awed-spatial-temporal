library(knitr)
library(here)
library(tidyverse)
library(kableExtra)
library(ggpubr)
library(latex2exp)
library(geodist)


source(here("lib","binary-matrix-function_odds.R"))
source(here("lib","bootstrap-function_odds.R"))
cluster_wrap_function <- function(df, clusters, perms, r1, r2, t = 30){
  df_cluster <- df[df$cluster == clusters,]
  
  output <- bootstrap_function_odds(df = df_cluster,
                                    r1 = r1, 
                                    r2 = r2,
                                    t = t,
                                    replace = FALSE
  )
  return(output)
}
load(here("data","cluster-specific-permutation-output.RData"))
load(here("data","cluster-specific-permutation-output-vcds-only.RData"))
load(here("data","2020-12-08_work-ts-spat.RData"))

# Observed Tau
# Increasing Concentric
r_upper <- as.vector(c(50,seq(100, 1000, by = 100)), mode = "list")
r_lower <- lapply(r_upper, function(x) {ifelse(x > 100, x - 100, x-50)})
names(r_upper) <- r_upper
names(r_lower) <- r_lower

clusters <- as.list(c(1,3:24))
names(clusters) <- clusters

observed_concentric <- map_dfr(clusters, ~cluster_wrap_function(
  cluster = .x,
  df = dplyr::filter(work_ts_spat, dengue == 1, serotype != "unk_serotype"),
  r1 = r_lower,
  r2 = r_upper,
  t = 30),
  .id = "clusterid"
)

observed_concentric <- observed_concentric %>% 
  dplyr::select(cluster = clusterid, tau, r_lower, r_upper) 

# # obs_mean_1 <- 
#   observed_concentric %>% 
#   dplyr::select(cluster, tau_obs = tau, r_lower, r_upper) %>% 
#   full_join(output_cluster_permutations$permuted_concentric) %>% 
#   dplyr::select(-a,-b,-c,-d) %>%
#   dplyr::rename("tau_perm" = "tau") %>% 
#   full_join(tx_list) %>% 
#   mutate(tau_obs = ifelse(is.nan(tau_obs) | is.infinite(tau_obs), 0, tau_obs)) %>% 
#   mutate(tau_obs_1 = tau_obs + 1) %>% 
#   group_by(intervention, r_upper) %>% 
#   mutate(geom_mean_1 = exp(mean(log(tau_obs_1)))-1) %>% 
#   rowwise() %>% 
#   mutate(tau_perm = ifelse(is.nan(tau_perm) | is.infinite(tau_perm), 0, tau_perm)) %>% 
#   mutate(tau_perm_1 = tau_perm + 1) %>% 
#   group_by(intervention, r_upper, permutation) %>% 
#   mutate(geom_mean_perm_1 = exp(mean(log(tau_perm_1))) - 1) %>% 
#   group_by(intervention, r_upper) %>% 
#   mutate(CI.l = quantile(geom_mean_perm_1, probs = 0.025),
#          CI.u = quantile(geom_mean_perm_1, probs = 0.975)) %>% 
#   ungroup() %>% 
#   dplyr::select(intervention, r_lower, r_upper, geom_mean_1, CI.l, CI.u) %>% 
#   distinct() %>% head()
#   mutate(method = "obs")

# Intervention labels
tx_list <- work_ts_spat %>% 
  dplyr::select(cluster, intervention) %>% 
  mutate(intervention = ifelse(intervention == 0, "Untreated", "Intervention")) %>% 
  distinct() %>% 
  mutate(cluster = factor(cluster, levels = 1:24))


# Permutation data
perm_c <- output_cluster_permutations_vcds$permuted_concentric %>% 
  dplyr::select(cluster = clusterid,
                tau_perm = tau,
                r_lower, r_upper)

# Merging for graphing purposes
perm_c <- full_join(observed_concentric, perm_c)

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
  ggplot(aes(x = r_upper, y = tau)) + 
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
  ggplot(aes(x = r_upper, y = tau)) + 
  facet_wrap(~cluster) + 
  geom_line() +
  geom_ribbon(aes(ymin = UL, ymax = LL),
              alpha = 0.2) + 
  xlab(TeX("$d_2")) + 
  ylab(TeX("Odds ratio ($\\tau)")) + 
  theme_pubr()


ggpubr::ggarrange(p1_u, p1_t, ncol = 1,
                  labels = "AUTO")

# Forest plots - UNTREATED
p2_u <- perm_c %>% 
  mutate(cluster = factor(cluster, levels = 1:24)) %>% 
  filter(intervention == "Untreated", r_upper <= 300) %>% 
  dplyr::select(-tau_perm) %>% 
  distinct() %>% 
  ggplot(aes(x = tau, y = cluster)) + 
  facet_wrap(~r_upper, nrow = 1) + 
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
  ggplot(aes(x = tau, y = cluster)) + 
  facet_wrap(~r_upper, nrow = 1) + 
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

ggarrange(p2_t, p2_u, ncol = 1,
          labels = "AUTO")
ggsave(filename = here("graphs",
                       paste0(Sys.Date(), "_cluster-specific-tau-perm-forest-plots-vcds-only.png")),
       device = "png",
       width = 7, 
       height = 7,
       units = "in")


geom_means <- perm_c %>% 
  rowwise() %>% 
  mutate(tau_1 = ifelse(is.na(tau) | is.infinite(tau), 1, tau + 1),
         tau_perm_1 = ifelse(is.na(tau_perm) | is.infinite(tau_perm), 1, tau_perm + 1)) %>% 
  group_by(intervention, r_upper) %>% 
  mutate(geom_mean_1 = exp(mean(log(tau_1)))-1) %>% 
  group_by(cluster) %>% 
  mutate(permutation = row_number()) %>% 
  group_by(permutation) %>% 
  mutate(geom_mean_perm_1 = exp(mean(log(tau_perm_1)))-1) %>% 
  group_by(intervention, r_upper) %>% 
  mutate(CI.l = quantile(geom_mean_perm_1, probs = 0.025),
         CI.u = quantile(geom_mean_perm_1, probs = 0.975)) %>%
  dplyr::select(r_lower, r_upper, intervention, geom_mean_1, CI.l, CI.u) %>% 
  distinct()

geom_means %>% 
  ggplot(aes(x = r_upper, y = geom_mean_1)) + 
  facet_wrap(~intervention, ncol = 1) +
  geom_line() +
  geom_ribbon(aes(ymin = CI.l, ymax = CI.u),
              alpha = 0.3) + 
  geom_hline(yintercept = 1,
             lty = 2) +
  theme_pubr()
