library(knitr)
library(tidyverse)
library(kableExtra)
library(ggpubr)


load(here("data", "cluster-specific-jackknife.RData"))
load(here("data", "cluster-specific-bootstrap-output.RData"))
load(here("data","2020-12-08_work-ts-spat.RData"))

observed <- output_cluster_specific$observed_concentric

######################
## Jackknife estimation
######################
observed.val <- observed %>% 
  dplyr::select(cluster = clusterid,
                observed.tau = tau,
                r_lower, r_upper)
jk.val <- output_jackknife %>% 
  dplyr::select(cluster, 
                jack.tau = tau,
                r_lower, r_upper, left.out)

df <- full_join(observed.val, jk.val)

## Example: Cluster 13
df %>% 
  filter(cluster == 13, r_upper <= 300) %>% 
  ggplot(aes(x = log(jack.tau))) + 
  geom_density() +
  facet_wrap(~r_upper,
             scales = "free",
             ncol = 1) + 
  geom_vline(aes(xintercept = log(observed.tau)),
             col = "red",
             lty = 3) + 
  theme_pubr()

temp <- df %>% 
  filter(cluster == 13, r_upper == 100)

# ALL ON THE LOG SCALE!
df <- df %>% 
  # Count N for each cluster and limit
  group_by(cluster, r_upper) %>% 
  mutate(N = n_distinct(left.out)) %>% 
  ungroup() %>% 
  rowwise() %>% 
  # Generate pseudo-values for tau for variance estimation:
  mutate(pseudo.log.tau = N*log(observed.tau) - (N-1)*log(jack.tau)) %>% 
  # FOR SAKE OF INTEREST
  filter(!is.infinite(pseudo.log.tau)) %>% 
  # Bias-corrected jackknife estimate of tau 
  group_by(cluster, r_upper) %>% 
  mutate(bc.log.tau = mean(pseudo.log.tau, na.rm = TRUE)) %>% 
  ungroup() %>% 
  # Jackknife variance of tau
  mutate(ss = (pseudo.log.tau - bc.log.tau)^2) %>% 
  group_by(cluster, r_upper) %>% 
  mutate(var.log.tau = (1/(N*(N-1)))*sum(ss)) %>% 
  ungroup() %>% 
  mutate(sd.log.tau = sqrt(var.log.tau)) %>% 
  dplyr::select(cluster, r_lower, r_upper, observed.tau, N, bc.log.tau, var.log.tau, sd.log.tau) %>% 
  distinct() %>% 
  # Symmetric CIs?
  mutate(jack.tau.bc = exp(bc.log.tau),
         tau.CI.l = exp(bc.log.tau - 1.96*sd.log.tau),
         tau.CI.u = exp(bc.log.tau + 1.96*sd.log.tau))

p1_100 <- df %>% 
  full_join(mutate(tx_list, cluster = factor(cluster, levels = 1:24)))  %>% 
  mutate(intervention = ifelse(intervention == 0, "Untreated", "Intervention")) %>% 
  filter(r_upper == 100) %>% 
  mutate(cluster = factor(cluster, levels = 1:24)) %>% 
  ggplot(aes(x = jack.tau.bc, y = cluster)) + 
  facet_wrap(~intervention,
             scales = "free_y",
             ncol = 1) +
  geom_point() + 
  geom_point(aes(x = observed.tau),
             shape = 17,
             col = "blue") + 
  geom_errorbarh(aes(xmin = tau.CI.l, xmax = tau.CI.u)) + 
  geom_vline(xintercept = 1, lty = 2, col = "darkgray") + 
  xlab(TeX("Odds ratio ($\\tau)")) +
  scale_x_log10() + 
  theme_pubr()

p1_200 <- df %>% 
  full_join(mutate(tx_list, cluster = factor(cluster, levels = 1:24)))  %>% 
  mutate(intervention = ifelse(intervention == 0, "Untreated", "Intervention")) %>% 
  filter(r_upper == 200) %>% 
  mutate(cluster = factor(cluster, levels = 1:24)) %>% 
  ggplot(aes(x = jack.tau.bc, y = cluster)) + 
  facet_wrap(~intervention,
             scales = "free_y",
             ncol = 1) +
  geom_point() + 
  geom_point(aes(x = observed.tau),
             shape = 17,
             col = "blue") + 
  geom_errorbarh(aes(xmin = tau.CI.l, xmax = tau.CI.u)) + 
  geom_vline(xintercept = 1, lty = 2, col = "darkgray") + 
  scale_x_log10() + 
  xlab(TeX("Odds ratio ($\\tau)")) +
  theme_pubr()

p1_300 <- df %>% 
  full_join(mutate(tx_list, cluster = factor(cluster, levels = 1:24)))  %>% 
  mutate(intervention = ifelse(intervention == 0, "Untreated", "Intervention")) %>% 
  filter(r_upper == 300) %>% 
  mutate(cluster = factor(cluster, levels = 1:24)) %>% 
  ggplot(aes(x = jack.tau.bc, y = cluster)) + 
  facet_wrap(~intervention,
             scales = "free_y",
             ncol = 1) +
  geom_point() + 
  geom_point(aes(x = observed.tau),
             shape = 17,
             col = "blue") + 
  geom_errorbarh(aes(xmin = tau.CI.l, xmax = tau.CI.u)) + 
  geom_vline(xintercept = 1, lty = 2, col = "darkgray") + 
  scale_x_log10() + 
  xlab(TeX("Odds ratio ($\\tau)")) +
  theme_pubr()
  
ggarrange(p1_100, p1_200, p1_300, 
          nrow = 1,
          labels = "AUTO")
ggsave(filename = here('graphs',
                       paste0(Sys.Date(), "_cluster-forest-concentric-jack.png")),
       device = "png",
       width = 9,
       height = 6,
       units = "in")

tx_list <- work_ts_spat %>% 
  dplyr::select(cluster, intervention) %>% 
  distinct()

df_sum_stats <- df %>% 
  dplyr::select(cluster, r_lower, r_upper, observed.tau, bc.tau, var.tau, sd.tau) %>% 
  distinct() %>% 
  mutate(cluster = as.numeric(cluster)) %>% 
  full_join(tx_list)

df_sum_stats %>% 
  mutate(cluster = factor(cluster, levels = 1:24)) %>% 
  ggplot(aes(x = r_upper, y = bc.tau, col = cluster)) + 
  geom_line() + 
  facet_wrap(~intervention) + 
  geom_smooth(aes(group = intervention)) + 
  geom_hline(yintercept = 1, lty = 2) + 
  theme_pubr() + 
  coord_cartesian(ylim = c(0,7.5))

df_sum_stats %>% 
  mutate(CI.l = )


# VISUALIZATION
output_jackknife %>% 
  filter(r_upper <= 500) %>% 
  ggplot(aes(x = tau, col = as.factor(r_upper))) + 
  facet_wrap(~ cluster, scales = "free") +
  geom_density() + 
  theme_pubr()
