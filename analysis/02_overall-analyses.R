#############
# Overall results
#############

library(tidyverse)
library(here)
library(ggpubr)
library(latex2exp)

load(here("data","overall-bootstrap-output.RData"))
load(here("data", "overall-permutation-output.RData"))
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
  # geom_ribbon(aes(ymin = CI.l, ymax = CI.u), alpha = 0.5) + 
  geom_errorbar(aes(ymin = CI.l, ymax = CI.u), col = "darkgray") +
  geom_point() +
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
  # geom_ribbon(aes(ymin = CI.l, ymax = CI.u), alpha = 0.5) + 
  geom_errorbar(aes(ymin = CI.l, ymax = CI.u), col = "darkgray") +
  # geom_line() +
  geom_point() + 
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
ggsave(filename = here("graphs",
                       paste0(Sys.Date(), "_overall-tau-dual-plot.pdf")),
       device = "pdf",
       dpi = 300,
       width = 6,
       height = 6,
       units = "in")

p2
ggsave(filename = here("graphs",
                       paste0(Sys.Date(), "_overall-tau-concentric-plot.pdf")),
       device = "pdf",
       dpi = 300,
       width = 6,
       height = 4,
       units = "in")

p1
ggsave(filename = here("graphs",
                       paste0(Sys.Date(), "_overall-tau-radii-plot.pdf")),
       device = "pdf",
       dpi = 300,
       width = 6,
       height = 4,
       units = "in")

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
  # geom_line() + 
  # geom_ribbon(aes(ymin = LB, ymax = UB),
  #             alpha = 0.5) + 
  geom_errorbar(aes(ymin = LB, ymax = UB), col = "darkgray") +
  geom_point() +
  geom_hline(yintercept = 1, lty = 2) + 
  scale_x_continuous(TeX("Distance (m), d_2"),
                     breaks = seq(0,1000,by = 100),
                     labels = seq(0,1000,by = 100)) + 
  coord_cartesian(xlim = c(50,1000),
                  ylim = c(0,18)) +
  scale_y_continuous(TeX("$\\tau_{overall} = \\frac{\\hat{\\theta}(d_1,d_2)}{\\hat{\\theta}(0,\\infty)}")) + 
  theme_pubr()

ggarrange(p2, p3, ncol = 1, labels = "AUTO")
ggsave(filename = here("graphs",
                       paste0(Sys.Date(), "_overall-tau-concentric-plot-perm-boot.pdf")),
       device = "pdf",
       dpi = 300,
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


### Processed data:
perm_c %>% 
  rename(perm.tau = tau) %>% 
  group_by(r_upper) %>% 
  mutate(LB = quantile(perm.tau, probs = 0.025, na.rm = TRUE),
         UB = quantile(perm.tau, probs = 0.975, na.rm = TRUE)) %>% 
  ungroup() %>% 
  dplyr::select(r_lower, r_upper, LB, UB) %>% 
  distinct() %>% 
  full_join(output$observed_concentric) %>%
  dplyr::select(r_lower, r_upper, tau, LB, UB) %>% 
  full_join(cis_c, by = c("r_lower", "r_upper")) %>% 
  filter(r_upper <= 1000) %>% 
  dplyr::select(d1 = r_lower, d2 = r_upper, tau, 
                `Bootstrap CI.l` = CI.l, `Bootstrap CI.u` = CI.u, 
                `Perm Null CI.l` = LB, `Perm Null CI.u` = UB) %>% 
  write.csv(file = here("graphs", paste0("processed-data/", Sys.Date(), "_fig3-data.csv")),
            row.names = FALSE)
