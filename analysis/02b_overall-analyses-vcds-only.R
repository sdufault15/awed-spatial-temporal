#############
# Overall results
#############

library(tidyverse)
library(here)
library(ggpubr)
library(latex2exp)

load(here("data","overall-bootstrap-output.RData"))
load(here("data","overall-permutation-output.RData"))
load(here("data", "overall-permutation-output-vcds-only.RData"))
load(here("data","2020-12-08_work-ts-spat.RData"))

# tx <- work_ts_spat %>% 
#   select(cluster, intervention) %>% 
#   distinct() %>% 
#   mutate(intervention = ifelse(intervention == 0, "Untreated", "Intervention"))

cis_r <- output$bootstrapped_radii %>% 
  group_by(r_upper) %>% 
  summarise(CI.l = quantile(tau, probs = 0.025, na.rm = TRUE),
            CI.u = quantile(tau, probs = 0.975, na.rm = TRUE))
## Overall 
p1 <- output$observed_radii %>% 
  full_join(cis_r, by = "r_upper") %>% 
  ggplot(aes(x = r_upper, y = tau)) + 
  geom_hline(yintercept = 1,
             lty = 2) + 
  geom_ribbon(aes(ymin = CI.l, ymax = CI.u), alpha = 0.3) + 
  geom_line() +
  theme_pubr() + 
  scale_x_continuous(TeX("Distance (m), d_2"),
                     breaks = seq(0,1000,by = 100),
                     labels = seq(0,1000,by = 100)) + 
  coord_cartesian(xlim = c(50,1000),
                  ylim = c(0,18)) +
  scale_y_continuous(TeX("$\\tau_{overall} = \\frac{\\hat{\\theta}(0,d_2)}{\\hat{\\theta}(0,\\infty)}")) + 
  theme(legend.position = "none")


cis_c <- output$bootstrapped_concentric %>% 
  group_by(r_lower, r_upper) %>% 
  summarise(CI.l = quantile(tau, probs = 0.025, na.rm = TRUE),
            CI.u = quantile(tau, probs = 0.975, na.rm = TRUE))
p2 <- output$observed_concentric %>% 
  full_join(cis_c, by = c("r_lower", "r_upper")) %>% 
  ggplot(aes(x = r_upper, y = tau)) + 
  geom_hline(yintercept = 1,
             lty = 2) + 
  geom_ribbon(aes(ymin = CI.l, ymax = CI.u), alpha = 0.3) + 
  geom_line() +
  theme_pubr() + 
  scale_x_continuous(TeX("Distance (m), d_2"),
                     breaks = seq(0,1000,by = 100),
                     labels = seq(0,1000,by = 100)) + 
  coord_cartesian(xlim = c(50,1000),
                  ylim = c(0,18)) +
  scale_y_continuous(TeX("$\\tau_{overall} = \\frac{\\hat{\\theta}(d_1,d_2)}{\\hat{\\theta}(0,\\infty)}")) + 
  theme(legend.position = "none")


ggpubr::ggarrange(p1, p2,
                  ncol = 1,
                  labels = "AUTO",
                  align = "hv")

p2

p1

### Permutation distribution
perm_c <- output_permuted$permuted_concentric
p3 <- perm_c %>% 
  rename(perm.tau = tau) %>% 
  group_by(r_upper) %>% 
  mutate(LB = quantile(perm.tau, probs = 0.025, na.rm = TRUE),
         UB = quantile(perm.tau, probs = 0.975, na.rm = TRUE)) %>% 
  ungroup() %>% 
  dplyr::select(r_lower, r_upper, LB, UB) %>% 
  distinct() %>% 
  full_join(output$observed_concentric) %>% 
  ggplot(aes(x = r_upper, y = tau)) + 
  geom_line() + 
  geom_ribbon(aes(ymin = LB, ymax = UB),
              alpha = 0.3) + 
  geom_hline(yintercept = 1, lty = 2) + 
  scale_x_continuous(TeX("Distance (m), d_2"),
                     breaks = seq(0,1000,by = 100),
                     labels = seq(0,1000,by = 100)) + 
  coord_cartesian(xlim = c(50,1000),
                  ylim = c(0,10)) +
  scale_y_continuous(TeX("$\\tau_{overall} = \\frac{\\hat{\\theta}(d_1,d_2)}{\\hat{\\theta}(0,\\infty)}")) + 
  theme_pubr()

ggarrange(p2, p3, ncol = 1, labels = "AUTO")
# ggsave(filename = here("graphs",
#                        paste0(Sys.Date(), "_overall-tau-concentric-plot-perm-boot.png")),
#        device = "png",
#        width = 8,
#        height = 6,
#        units = "in")

### Permutation distribution - VCDs only
perm_c_vcd <- output_permuted_vcds$permuted_concentric
# Increasing Concentric
source("lib/bootstrap-function_odds.R")
source("lib/binary-matrix-function_odds.R")
r_upper <- as.vector(c(50,seq(100, 1000, by = 100)), mode = "list")
r_lower <- lapply(r_upper, function(x) {ifelse(x > 100, x - 100, x-50)})
names(r_upper) <- r_upper
names(r_lower) <- r_lower

library(geodist)
observed_concentric <- bootstrap_function_odds(
  df = dplyr::filter(work_ts_spat, dengue == 1, serotype != "unk_serotype"),
  r1 = r_lower,
  r2 = r_upper,
  t = 30,
  replace = FALSE)

p4 <- perm_c_vcd %>% 
  rename(perm.tau = tau) %>% 
  group_by(r_upper) %>% 
  mutate(LB = quantile(perm.tau, probs = 0.025, na.rm = TRUE),
         UB = quantile(perm.tau, probs = 0.975, na.rm = TRUE)) %>% 
  ungroup() %>% 
  dplyr::select(r_lower, r_upper, LB, UB) %>% 
  distinct() %>% 
  full_join(observed_concentric) %>% 
  ggplot(aes(x = r_upper, y = tau)) + 
  geom_line() + 
  geom_ribbon(aes(ymin = LB, ymax = UB),
              alpha = 0.3) + 
  geom_hline(yintercept = 1, lty = 2) + 
  scale_x_continuous(TeX("Distance (m), d_2"),
                     breaks = seq(0,1000,by = 100),
                     labels = seq(0,1000,by = 100)) + 
  coord_cartesian(xlim = c(50,1000),
                  ylim = c(0,10)) +
  scale_y_continuous(TeX("$\\tau_{overall} = \\frac{\\hat{\\theta}(d_1,d_2)}{\\hat{\\theta}(0,\\infty)}")) + 
  theme_pubr()
p3a <- p3 + ggtitle("Including test-negatives")
p4a <- p4 + ggtitle("Excluding test-negatives")
ggarrange(p3a, p4a, ncol = 1, labels = "AUTO")
ggsave(filename = here("graphs",
                       paste0(Sys.Date(), "_overall-tau-concentric-plot-perm-all-v-vds.png")),
       device = "png",
       width = 8,
       height = 6,
       units = "in")
# Tables
library(kableExtra)
library(knitr)
t1 <- output$observed_radii %>% 
  full_join(cis_r, by = "r_upper") %>% 
  filter(r_upper <= 1000) %>% 
  dplyr::select(r_lower, r_upper, tau, CI.l, CI.u) %>% 
  transmute(
    `d<sub>2</sub> (m)` = r_upper, 
    #`Interval (m)` = paste0("(", `r_lower`, ", ", `r_upper`, ")"),
    `Tau<sup>1</sup> (95% CI)` = paste0(format(round(tau, 2), nsmall = 2),
                            " (", 
                            format(round(CI.l, 2), nsmall = 2),
                            ", ",
                            format(round(CI.u, 2), nsmall = 2),
                            ")")) 

t2 <- output$observed_concentric %>% 
  full_join(cis_c, by = c("r_lower", "r_upper")) %>% 
  filter(r_upper <= 1000) %>% 
  dplyr::select(r_lower, r_upper, tau, CI.l, CI.u) %>% 
  transmute(#`Interval (m)` = paste0("(", `r_lower`, ", ", `r_upper`, ")"),
            `Tau<sup>2</sup> (95% CI)` = paste0(format(round(tau, 2), nsmall = 2),
                                    " (", 
                                    format(round(CI.l, 2), nsmall = 2),
                                    ", ",
                                    format(round(CI.u, 2), nsmall = 2),
                                    ")")) 
bind_cols(t1,t2) %>% 
  kable(escape = F) %>% 
  kable_styling() %>% 
  footnote(number = c("Distance intervals (0, d<sub>2</sub>)", "Distance intervals (d<sub>2</sub> - 100, d<sub>2</sub>)"),
           escape = F)
