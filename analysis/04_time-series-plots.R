library(tidyverse)
library(here)
library(ggpubr)
library(latex2exp)
load(here("data","2020-12-08_work-ts-spat.RData"))

sero_col <- c("#0AA4D1", "#F3C73E", "#EF7822", "#76C1A8", "#003D58")
p1 <- work_ts_spat %>% 
  filter(serotype == "test-negative control") %>% 
  group_by(date = as.Date(cut(illness_onset, "month")), intervention = as.factor(intervention)) %>% 
  summarise(`Test-negative counts (monthly)` = n_distinct(participant_id)) %>% 
  ungroup() %>% 
  ggplot(aes(x = date, 
             y = `Test-negative counts (monthly)`, 
             col = intervention)) + 
  # geom_point() + 
  geom_line(aes(group = intervention, lty = intervention),
            size = 0.75,
            col = "black",
            alpha = 0.8) + 
  # scale_color_manual(values = sero_col[c(1,4)],#wes_palette("Zissou1")[c(5,1)],
  #                    breaks = c("0", "1"),
  #                    labels = c("Untreated", "Intervention")) +
  scale_linetype_manual(values = c(1,2),#wes_palette("Zissou1")[c(5,1)],
                     breaks = c("0", "1"),
                     labels = c("Untreated", "Intervention")) +
  ggpubr::theme_pubr() + 
  #scale_x_discrete(name = "Month of illness onset") +  
  scale_x_date(name = "Month of illness onset",
               date_breaks = "2 months",
               date_labels = "%b %Y ") + 
  theme(axis.text.x = element_text(angle = -47,
                                   hjust = 0,
                                   vjust = 0,
                                   color = "black"),
        axis.text = element_text(size = 12),
        axis.title = element_text(size = 12),
        legend.title = element_blank(),
        legend.text = element_text(size = 10),
        axis.text.y = element_text(color = "black"),
        axis.title.x = element_text(color = "black"),
        plot.margin = unit(c(5.5, 20, 5.5, 5.5), "points"))

p3 <- work_ts_spat %>% 
  filter(serotype != "test-negative control") %>% 
  group_by(date = as.Date(cut(illness_onset, "month")), intervention, serotype) %>% 
  summarise(counts = n_distinct(participant_id)) %>% 
  mutate(counts = replace_na(counts, 0),
         intervention = case_when(intervention == 0 ~ "Untreated",
                                  intervention == 1 ~ "Intervention")) %>% 
  group_by(date, intervention) %>% 
  mutate(totals = sum(counts)) %>% 
  # Missing september 2018, and february
  # bind_rows(data.frame(date = as.Date(c(rep("2018-02-01",2),rep("2018-09-01",2))),
  #                      intervention = rep(c("Intervention", "Untreated"),2),
  #                      serotype = rep(c("unk_serotype","unk_serotype"),2),
  #                      counts = rep(0,4),
  #                      totals = rep(0,4))) %>% 
  ggplot(aes(x = date, y = counts, col = serotype, fill = serotype)) +
  facet_wrap(~intervention,
             ncol = 1) + 
  geom_bar(stat = "identity") +
  scale_color_manual(name = "Serotype",
                     aesthetics = c("color", "fill"),
                     values = sero_col,
                     labels = c(paste0("DENV", 1:4), "Unknown")) + 
  ylab("VCD counts (monthly)") + 
  scale_x_date(name = "Month of illness onset",
               date_breaks = "2 months",
               date_labels = "%b %Y ") + 
  ggpubr::theme_pubr() + 
  theme(axis.text.x = element_text(angle = -47,
                                   hjust = 0,
                                   vjust = 0,
                                   color = "black"),
        axis.text.y = element_text(color = "black"),
        axis.title.x = element_blank(),
        axis.text = element_text(size = 12),
        axis.title = element_text(size = 12),
        legend.text = element_text(size = 10),
        legend.title = element_text(size = 12),
        strip.text = element_text(size = 14))


ggpubr::ggarrange(p1, p3,
                  align = "v",
                  ncol = 1,
                  heights = c(1,1.5),
                  labels = c("A", "B"))
ggsave(filename = here("graphs", paste0(Sys.Date(), "_time-series-plot-serotypes.png")),
       device = "png",
       width = 8.75,
       height = 8,
       units = "in")

### Processed data
work_ts_spat %>% 
 # filter(serotype == "test-negative control") %>% 
  group_by(date = as.Date(cut(illness_onset, "month")), intervention = as.factor(intervention), serotype) %>% 
  arrange(date, intervention, serotype) %>%
  summarise(`Counts (monthly)` = n_distinct(participant_id)) %>% 
  mutate(intervention = ifelse(intervention == 0, "Untreated", "Intervention")) %>% 
  write.csv(file = here("graphs", paste0("processed-data/", Sys.Date(), "_fig1-data.csv")),
            row.names = FALSE)
       