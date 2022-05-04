library(tidyverse)
library(here)
library(ggpubr)
library(latex2exp)
library(knitr)
library(kableExtra)

load(here("data","2020-12-08_work-ts-spat.RData"))

load(here("data","2021-01-07_distance-all-pts-to-boundary.RData")) # with distance to boundaries
load(here("data","2021-01-08_distance-all-pts-to-boundary-sensitivity.RData")) # with distance to boundaries

names(dist_all)
t1 <- dist_all %>% 
  mutate(dengue = ifelse(dengue == 0, "Test-negatives", "VCDs"),
         intervention = ifelse(intervention == 0, "Untreated", "Intervention")) %>% 
  group_by(intervention, cluster, dengue) %>% 
  summarise(N = n(),
            avg_distance = mean(distance),
            sd_distance = sd(distance),
            med_distance = median(distance)) %>% 
  pivot_wider(id_cols = c(intervention, cluster),
              names_from = dengue,
              values_from = N:med_distance)

t2 <- sensitivity_data$dist_all_s50 %>% 
  mutate(dengue = ifelse(dengue == 0, "Test-negatives", "VCDs"),
         intervention = ifelse(intervention == 0, "Untreated", "Intervention")) %>% 
  group_by(intervention, cluster, dengue) %>% 
  summarise(N = n(),
            avg_distance = mean(distance),
            sd_distance = sd(distance),
            med_distance = median(distance)) %>% 
  pivot_wider(id_cols = c(intervention, cluster),
              names_from = dengue,
              values_from = N:med_distance)

t2 <- t2 %>% 
  dplyr::select_at(vars(c(intervention, cluster, contains("N_")))) %>% 
  dplyr::rename("N<sup>&dagger;</sup>" = `N_VCDs`,
                "N<sup>&dagger;</sup> tn" = `N_Test-negatives`)

t1_new <- full_join(t1, t2)
t1_new <- t1_new[,c(1,2,4,12,6,8,3,11,5,7)]

names(t1_new) <- c("Intervention", "Cluster", 
                   "N<sup>*</sup>", "N<sup>&dagger;</sup>",
                   "Mean Distance<sup>*</sup> (m)", "SD Distance<sup>*</sup> (m)",
                   "N<sup>*</sup>", "N<sup>&dagger;</sup>",
                   "Mean Distance<sup>*</sup> (m)", "SD Distance<sup>*</sup> (m)")

t1_new %>% 
  kable(digits = 2,
        escape = FALSE,
        align = "c") %>% 
  kable_styling() %>% 
  collapse_rows(columns = 1) %>% 
  add_header_above(header = c(" " = 2, "VCDs" = 4, "Test-negatives" = 4)) %>% 
  footnote(symbol = c("Full data", "Sensitivity data: removing individuals within 50m of the cluster boundary"),
           escape = FALSE)

write_csv(t1_new,
          path = here("data", paste0(Sys.Date(), "_table-distance-sensitivity.csv")))
          
          