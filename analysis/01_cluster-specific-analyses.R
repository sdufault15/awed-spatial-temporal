#############
# Cluster-specific results
#############

library(tidyverse)
library(here)
library(ggpubr)
library(latex2exp)

load(here("data","cluster-specific-bootstrap-output.RData"))
load(here("data","2020-12-08_work-ts-spat.RData"))

tx <- work_ts_spat %>% 
  select(cluster, intervention) %>% 
  distinct() %>% 
  mutate(intervention = ifelse(intervention == 0, "Untreated", "Intervention"))

## Setting up the data
# Increasing radii
df_radii <- output_cluster_specific$observed_radii %>% 
  full_join(mutate(tx, cluster = paste0(cluster)), by = c("clusterid" = "cluster"))
cis_r <- output_cluster_specific$bootstrapped_radii %>% 
  group_by(clusterid, r_upper) %>% 
  summarise(CI.l = quantile(tau, probs = 0.025, na.rm = TRUE),
            CI.u = quantile(tau, probs = 0.975, na.rm = TRUE))
df_radii <- full_join(df_radii,
                      cis_r, by = c("clusterid", "r_upper"))
# Concentric rings
df_concentric <- output_cluster_specific$observed_concentric %>% 
  full_join(mutate(tx, cluster = paste0(cluster)), by = c("clusterid" = "cluster"))
cis_c <- output_cluster_specific$bootstrapped_concentric %>% 
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
  ggplot(aes(x = r_upper, y = tau, shape = clusterid)) + 
  facet_wrap(~intervention, ncol = 1, scales = "free_y") +
  geom_hline(yintercept = 1,
             lty = 2) + 
  #geom_line() +
  scale_shape_manual(values = 1:24) +
  geom_point(position = position_jitter(width = 4, height = 0)) + 
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
                       paste0(Sys.Date(), "_cluster-specific-tau-dual-plot.png")),
       height = 5,
       width = 8,
       units = "in",
       device = "png")

ggsave(plot = p2, 
       filename = here("graphs", 
                       paste0(Sys.Date(), "_cluster-specific-tau-concentric-plot.png")),
       height = 5,
       width = 8,
       units = "in",
       device = "png")

# Distance specific forest plots
p1a <- df_concentric %>% 
  filter(r_upper == 100) %>% 
  mutate(clusterid = factor(clusterid, levels = c(paste0(1:24)))) %>% 
  ggplot(aes(y = clusterid)) +
  facet_wrap(~intervention, ncol = 1, scales = "free_y") +
  geom_point(aes(x = tau)) + 
  geom_errorbarh(aes(xmin = CI.l, xmax = CI.u)) + 
  geom_vline(xintercept = 1, lty = 2) + 
  xlab(TeX("$\\tau_{cluster} = \\frac{\\hat{\\theta}(d_1,d_2)}{\\hat{\\theta}(0,\\infty)}")) + 
  ylab("Cluster ID") + 
  theme_pubr()

p1b <- df_concentric %>% 
  filter(r_upper == 200) %>% 
  mutate(clusterid = factor(clusterid, levels = c(paste0(1:24)))) %>% 
  ggplot(aes(y = clusterid)) +
  facet_wrap(~intervention, ncol = 1, scales = "free_y") +
  geom_point(aes(x = tau)) + 
  geom_errorbarh(aes(xmin = CI.l, xmax = CI.u)) + 
  geom_vline(xintercept = 1, lty = 2) + 
  xlab(TeX("$\\tau_{cluster} = \\frac{\\hat{\\theta}(d_1,d_2)}{\\hat{\\theta}(0,\\infty)}")) + 
  ylab("Cluster ID") + 
  theme_pubr()

p1c <- df_concentric %>% 
  filter(r_upper == 300) %>%
  mutate(clusterid = factor(clusterid, levels = c(paste0(1:24)))) %>% 
  ggplot(aes(y = clusterid)) +
  facet_wrap(~intervention, ncol = 1, scales = "free_y") +
  geom_point(aes(x = tau)) + 
  geom_errorbarh(aes(xmin = CI.l, xmax = CI.u)) + 
  geom_vline(xintercept = 1, lty = 2) + 
  xlab(TeX("$\\tau_{cluster} = \\frac{\\hat{\\theta}(d_1,d_2)}{\\hat{\\theta}(0,\\infty)}")) + 
  ylab("Cluster ID") + 
  theme_pubr()

ggarrange(p1a, p1b, p1c, ncol = 3,
          labels = "AUTO")
ggsave(filename = here('graphs',
                       paste0(Sys.Date(), "_cluster-forest-concentric.png")),
       device = "png",
       width = 8,
       height = 6,
       units = "in")

# Radii
p1a_r <- df_radii %>% 
  filter(r_upper == 100) %>% 
  mutate(clusterid = factor(clusterid, levels = c(paste0(1:24)))) %>% 
  ggplot(aes(y = clusterid)) +
  facet_wrap(~intervention, ncol = 1, scales = "free_y") +
  geom_point(aes(x = tau)) + 
  geom_errorbarh(aes(xmin = CI.l, xmax = CI.u)) + 
  geom_vline(xintercept = 1, lty = 2) + 
  theme_pubr()

p1b_r <- df_radii %>% 
  filter(r_upper == 200) %>% 
  mutate(clusterid = factor(clusterid, levels = c(paste0(1:24)))) %>% 
  ggplot(aes(y = clusterid)) +
  facet_wrap(~intervention, ncol = 1, scales = "free_y") +
  geom_point(aes(x = tau)) + 
  geom_errorbarh(aes(xmin = CI.l, xmax = CI.u)) + 
  geom_vline(xintercept = 1, lty = 2) + 
  theme_pubr()

p1c_r <- df_radii %>% 
  filter(r_upper == 300) %>%
  mutate(clusterid = factor(clusterid, levels = c(paste0(1:24)))) %>% 
  ggplot(aes(y = clusterid)) +
  facet_wrap(~intervention, ncol = 1, scales = "free_y") +
  geom_point(aes(x = tau)) + 
  geom_errorbarh(aes(xmin = CI.l, xmax = CI.u)) + 
  geom_vline(xintercept = 1, lty = 2) + 
  theme_pubr()

ggarrange(p1a_r, 
          p1b_r, 
          p1c_r, ncol = 3,
          labels = "AUTO")
ggsave(filename = here('graphs',
                       paste0(Sys.Date(), "_cluster-forest-radii.png")),
       device = "png",
       width = 8,
       height = 6,
       units = "in")


## Meta-analysis
## Required packages in this paper 
lib2install <- c("metaSEM", "metafor") ## Install them automatically if they are not available in your computer 
for (i in lib2install) { 
  if (!(i %in% rownames(installed.packages()))) install.packages(i)
}
## Load the libraries 
library(metaSEM) 
library(metafor)

names(df_concentric)
# Setting up dataset for meta-analysis
test <- df_concentric %>% 
  #filter(r_upper == 100) %>% 
  mutate(log.tau = log(tau)) 

se.boot <- output_cluster_specific$bootstrapped_concentric %>% 
  filter(tau > 0) %>% 
  group_by(r_upper, clusterid) %>% 
  summarise(var.log.tau = var(log(tau)))

test2 <- full_join(test, se.boot) %>% 
  filter(!is.na(var.log.tau)) 

test_100m <- test2 %>% filter(r_upper == 100)
test_200m <- test2 %>% filter(r_upper == 200)
test_300m <- test2 %>% filter(r_upper == 300)

# 100 m: pooled tau = exp(1.25) = 3.49
summary(meta(y = log.tau, v = var.log.tau, data = test_100m[test_100m$intervention == "Untreated",]))
forest(rma(yi = log.tau, vi = var.log.tau, data = test_100m[test_100m$intervention == "Untreated",], slab = clusterid))

# 200 m: pooled tau = exp(0.90) = 2.46
summary(meta(y = log.tau, v = var.log.tau, data = test_200m[test_200m$intervention == "Untreated",]))
forest(rma(yi = log.tau, vi = var.log.tau, data = test_200m[test_200m$intervention == "Untreated",], slab = clusterid))

# 300 m: didn't converge
summary(meta(y = log.tau, v = var.log.tau, data = test_300m[test_300m$intervention == "Untreated",]))
forest(rma(yi = log.tau, vi = var.log.tau, data = test_300m[test_300m$intervention == "Untreated",], slab = clusterid))
