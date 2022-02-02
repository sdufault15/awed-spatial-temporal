#############
# Cluster-specific results
#############

library(tidyverse)
library(here)
library(ggpubr)
library(latex2exp)

load(here("data","overall-arm-specific-bootstrap-output.RData"))
load(here("data","2020-12-08_work-ts-spat.RData"))

# tx <- work_ts_spat %>% 
#   select(cluster, intervention) %>% 
#   distinct() %>% 
#   mutate(intervention = ifelse(intervention == 0, "Untreated", "Intervention"))

cis_r_tx <- output_arm_specific$bootstrapped_radii$tau_treat %>% 
  group_by(r_upper) %>% 
  summarise(CI.l = quantile(tau, probs = 0.025, na.rm = TRUE),
            CI.u = quantile(tau, probs = 0.975, na.rm = TRUE)) %>% 
  mutate(arm = "Intervention")

cis_r_unt <- output_arm_specific$bootstrapped_radii$tau_unt %>% 
  group_by(r_upper) %>% 
  summarise(CI.l = quantile(tau, probs = 0.025, na.rm = TRUE),
            CI.u = quantile(tau, probs = 0.975, na.rm = TRUE)) %>% 
  mutate(arm = "Untreated")

## Overall 
p1 <- bind_rows(output_arm_specific$observed_radii$tau_treat,
  output_arm_specific$observed_radii$tau_unt) %>%
  mutate(arm = c(rep("Intervention", nrow(output_arm_specific$observed_radii$tau_treat)),
                 rep("Untreated", nrow(output_arm_specific$observed_radii$tau_unt)))) %>% 
  full_join(bind_rows(cis_r_tx, cis_r_unt), by = c("r_upper", "arm")) %>% 
  ggplot(aes(x = r_upper, y = tau, col = arm,
             fill = arm)) + 
  geom_hline(yintercept = 1,
             lty = 2) + 
  geom_ribbon(aes(ymin = CI.l, ymax = CI.u), alpha = 0.3) + 
  geom_line() +
  theme_pubr() + 
  scale_x_continuous(TeX("Distance (m), d_2"),
                     breaks = seq(0,1000,by = 100),
                     labels = seq(0,1000,by = 100)) + 
  coord_cartesian(xlim = c(50,1000),
                  ylim = c(0,5)) +
  scale_y_continuous(TeX("$\\tau_{overall} = \\frac{\\hat{\\theta}(0,d_2)}{\\hat{\\theta}(0,\\infty)}")) + 
  theme(legend.position = "bottom")

cis_c_tx <- output_arm_specific$bootstrapped_concentric$tau_treat %>% 
  group_by(r_upper) %>% 
  summarise(CI.l = quantile(tau, probs = 0.025, na.rm = TRUE),
            CI.u = quantile(tau, probs = 0.975, na.rm = TRUE)) %>% 
  mutate(arm = "Intervention")

cis_c_unt <- output_arm_specific$bootstrapped_concentric$tau_unt %>% 
  group_by(r_upper) %>% 
  summarise(CI.l = quantile(tau, probs = 0.025, na.rm = TRUE),
            CI.u = quantile(tau, probs = 0.975, na.rm = TRUE)) %>% 
  mutate(arm = "Untreated")

p2 <- bind_rows(output_arm_specific$observed_concentric$tau_treat,
                output_arm_specific$observed_concentric$tau_unt) %>%
  mutate(arm = c(rep("Intervention", nrow(output_arm_specific$observed_concentric$tau_treat)),
                 rep("Untreated", nrow(output_arm_specific$observed_concentric$tau_unt)))) %>% 
  full_join(bind_rows(cis_c_tx, cis_c_unt), by = c("r_upper", "arm")) %>% 
  ggplot(aes(x = r_upper, y = tau, col = arm, fill = arm)) + 
  geom_hline(yintercept = 1,
             lty = 2) + 
  geom_ribbon(aes(ymin = CI.l, ymax = CI.u), alpha = 0.3) + 
  geom_line() +
  theme_pubr() + 
  scale_x_continuous(TeX("Distance (m), d_2"),
                     breaks = seq(0,1000,by = 100),
                     labels = seq(0,1000,by = 100)) + 
  coord_cartesian(xlim = c(50,1000),
                  ylim = c(0,5)) +
  scale_y_continuous(TeX("$\\tau_{overall} = \\frac{\\hat{\\theta}(d_2 - 100,d_2)}{\\hat{\\theta}(0,\\infty)}")) + 
  theme(legend.position = "bottom")


ggpubr::ggarrange(p1, p2,
                  ncol = 1,
                  labels = "AUTO",
                  align = "hv",
                  common.legend = TRUE,
                  legend = "bottom")

ggsave(filename = here("graphs",
                       paste0(Sys.Date(), "_arm-specific-tau-dual-plot.png")),
       device = "png",
       width = 6,
       height = 6,
       units = "in")

p2
ggsave(filename = here("graphs",
                       paste0(Sys.Date(), "_arm-specific-tau-concentric-plot.png")),
       device = "png",
       width = 6,
       height = 4,
       units = "in")


# Tables
library(kableExtra)
library(knitr)
t1 <- bind_rows(output_arm_specific$observed_radii$tau_treat,
                output_arm_specific$observed_radii$tau_unt) %>%
  mutate(arm = c(rep("Intervention", nrow(output_arm_specific$observed_radii$tau_treat)),
                 rep("Untreated", nrow(output_arm_specific$observed_radii$tau_unt)))) %>% 
  full_join(bind_rows(cis_r_tx, cis_r_unt), by = c("r_upper", "arm")) %>% 
  filter(r_upper <= 1000) %>% 
  dplyr::select(arm, r_lower, r_upper, tau, CI.l, CI.u) %>% 
  transmute(
    Arm = arm,
    `d<sub>2</sub> (m)` = r_upper, 
    #`Interval (m)` = paste0("(", `r_lower`, ", ", `r_upper`, ")"),
    `Tau<sup>1</sup> (95% CI)` = paste0(format(round(tau, 2), nsmall = 2),
                            " (", 
                            format(round(CI.l, 2), nsmall = 2),
                            ", ",
                            format(round(CI.u, 2), nsmall = 2),
                            ")")) 

t2 <- bind_rows(output_arm_specific$observed_concentric$tau_treat,
                output_arm_specific$observed_concentric$tau_unt) %>%
  mutate(arm = c(rep("Intervention", nrow(output_arm_specific$observed_concentric$tau_treat)),
                 rep("Untreated", nrow(output_arm_specific$observed_concentric$tau_unt)))) %>% 
  full_join(bind_rows(cis_c_tx, cis_c_unt), by = c("r_upper", "arm")) %>% 
  dplyr::select(r_lower, r_upper, tau, CI.l, CI.u) %>% 
  transmute(#`Interval (m)` = paste0("(", `r_lower`, ", ", `r_upper`, ")"),
    `Tau<sup>2</sup> (95% CI)` = paste0(format(round(tau, 2), nsmall = 2),
                                        " (", 
                                        format(round(CI.l, 2), nsmall = 2),
                                        ", ",
                                        format(round(CI.u, 2), nsmall = 2),
                                        ")")) 
bind_cols(t1,t2) %>% 
  pivot_wider(names_from = Arm, 
              values_from = c(`Tau<sup>1</sup> (95% CI)`,`Tau<sup>2</sup> (95% CI)`),
              names_repair = "minimal") %>% 
  kable(escape = F) %>% 
  kable_styling() %>% 
  footnote(number = c("Distance intervals (0, d<sub>2</sub>)", "Distance intervals (d<sub>2</sub> - 100, d<sub>2</sub>)"),
           escape = F) 
