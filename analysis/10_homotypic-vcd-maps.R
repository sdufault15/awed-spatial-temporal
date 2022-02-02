library(here)
library(tidyverse)
library(tmap)
library(geodist)
load(here("data","2020-12-08_work-ts-spat.RData"))

library(rgdal)
library(raster)
library(sp)
Yogya_Adm_1 <- readOGR("~/Box/2019-Research/WMP/CR-TND/awed-spatial-temporal/yogya-shape-files/", "Yk_CaseControl_20160512")
df_col_sero <- data.frame(serotype = c(paste0("DENV", 1:4), "Unknown", "Test-negative"),
                          sero.color = c("#0AA4D1", "#F3C73E", "#EF7822", "#76C1A8", "#003D58", "gray"),
                          sero.shape = c(21:24, 4, 25))
work_ts_spat <- work_ts_spat %>% 
  mutate(serotype = case_when(serotype == "denv1" ~ "DENV1",
                              serotype == "denv2" ~ "DENV2",
                              serotype == "denv3" ~ "DENV3",
                              serotype == "denv4" ~ "DENV4",
                              serotype == "unk_serotype"~ "Unknown",
                              serotype == "test-negative control" ~ "Test-negative"),
         intervention = ifelse(intervention == 0, "Untreated", "Intervention")) %>% 
  full_join(df_col_sero) 

int_assignment <- work_ts_spat %>% 
  dplyr::select(cluster, Intervention = intervention) %>% 
  distinct() %>% 
  mutate(New_Clustr = ifelse(cluster < 10, paste0("0", cluster), as.character(cluster))) %>% 
  dplyr::select(-cluster)
Yogya_Adm_1@data <- left_join(Yogya_Adm_1@data, int_assignment, by = "New_Clustr")

tm_shape(Yogya_Adm_1) + 
  tm_polygons(col = "Intervention",
              alpha = 0.75,
              palette = c("lightgray", "white")) + 
  tm_text(text = "New_Clustr")

# Creating Spatial Points Dataframe
CRTND_dengue_data_SPDF <- SpatialPointsDataFrame(coords = work_ts_spat[,c("longitude", "latitude")],
                                                 data = work_ts_spat[, c("participant_id", "cluster", "intervention", "enrolment_date", "illness_onset", "dengue", "age", "sex", "serotype", "sero.color", "sero.shape")],
                                                 proj4string = CRS("+init=epsg:4326")) # sets the projection to WGS 1984 using lat/long. Optional but good to specify

CRTND_VCD_data_SPDF <- subset(CRTND_dengue_data_SPDF, dengue == 1)
Yogya_Adm_1 <- spTransform(Yogya_Adm_1, crs(CRTND_dengue_data_SPDF))

## Finding distances
work_all_dengue <- work_ts_spat %>% 
  filter(dengue == 1 & serotype != "Unknown")
dist.mat <- geodist(x = cbind(longitude = work_all_dengue$longitude, latitude = work_all_dengue$latitude))
time.pairs <- expand.grid(work_all_dengue$illness_onset, work_all_dengue$illness_onset)
time.dif <- as.Date(time.pairs[,1]) - as.Date(time.pairs[,2])
time.mat <- matrix(time.dif, nrow = nrow(work_all_dengue))
# Serotype
sero.pairs <- expand.grid(work_all_dengue$serotype, work_all_dengue$serotype)
sero.related <- 1*(sero.pairs[,1] == sero.pairs[,2])
sero.mat <- matrix(sero.related, nrow = nrow(work_all_dengue))

diag(sero.mat) <- diag(time.mat) <- diag(dist.mat) <- NA

rownames(sero.mat) <- rownames(dist.mat) <- rownames(time.mat) <- work_all_dengue$participant_id
colnames(sero.mat) <- colnames(dist.mat) <- colnames(time.mat) <- work_all_dengue$participant_id

####################
# 50m
####################
# Exposed VCDs
t <- 30
r <- 50
# Time matrix
vcds_time.bin = 1*(time.mat <= t & time.mat > 0)
# Distance matrix
vcds_dist.bin = 1*(dist.mat < r)

sum(vcds_dist.bin, na.rm = TRUE)
sum(vcds_time.bin, na.rm = TRUE)

sum(vcds_dist.bin*vcds_time.bin, na.rm = TRUE)

ones_row <- apply(vcds_dist.bin*vcds_time.bin*sero.mat, 1, function(x){max(x,na.rm = TRUE)})
ones_col <- apply(vcds_dist.bin*vcds_time.bin*sero.mat, 2, function(x){max(x,na.rm = TRUE)})
temp <- names(ones_row)[which(ones_row == 1)]
temp_2 <- names(ones_col[which(ones_col == 1)])
m50 <- tm_shape(Yogya_Adm_1) + 
  tm_polygons(col = "Intervention",
              alpha = 0.75,
              palette = c("lightgray", "white")) + 
  tm_text(text = "New_Clustr") +
  tm_shape(subset(CRTND_VCD_data_SPDF, participant_id %in% unique(c(temp, temp_2)))) +
  tm_symbols(col = "serotype",
             shape = "serotype",
             legend.col.show = FALSE,
             legend.shape.show = FALSE,
             palette = df_col_sero$sero.color,
             shapes = df_col_sero$sero.shape,
             size = 0.2) + 
  tm_layout(title = "A") +
  tm_add_legend("symbol", 
                col = df_col_sero$sero.color[1:4], 
                shape = df_col_sero$sero.shape[1:4],
                labels = df_col_sero$serotype[1:4]) +
  tm_scale_bar(breaks = c(0,.5,1), text.size = 0.9, position = c("left", "top"))

####################
# 100m
####################
# Exposed VCDs
t <- 30
r <- 100
# Time matrix
vcds_time.bin = 1*(time.mat <= t & time.mat > 0)
# Distance matrix
vcds_dist.bin = 1*(dist.mat < r)

sum(vcds_dist.bin, na.rm = TRUE)
sum(vcds_time.bin, na.rm = TRUE)

sum(vcds_dist.bin*vcds_time.bin, na.rm = TRUE)

ones_row <- apply(vcds_dist.bin*vcds_time.bin*sero.mat, 1, function(x){max(x,na.rm = TRUE)})
ones_col <- apply(vcds_dist.bin*vcds_time.bin*sero.mat, 2, function(x){max(x,na.rm = TRUE)})
temp <- names(ones_row)[which(ones_row == 1)]
temp_2 <- names(ones_col[which(ones_col == 1)])
m100 <- tm_shape(Yogya_Adm_1) + 
  tm_polygons(col = "Intervention",
              alpha = 0.75,
              palette = c("lightgray", "white")) + 
  tm_text(text = "New_Clustr") +
  tm_shape(subset(CRTND_VCD_data_SPDF, participant_id %in% unique(c(temp, temp_2)))) +
  tm_symbols(col = "serotype",
             shape = "serotype",
             legend.col.show = FALSE,
             legend.shape.show = FALSE,
             palette = df_col_sero$sero.color,
             shapes = df_col_sero$sero.shape,
             size = 0.2) + 
  tm_add_legend("symbol", 
                col = df_col_sero$sero.color[1:4], 
                shape = df_col_sero$sero.shape[1:4],
                labels = df_col_sero$serotype[1:4]) +
  tm_layout(title = "B") +
  tm_scale_bar(breaks = c(0,.5,1), text.size = 0.9, position = c("left", "top"))

both <- tmap_arrange(m50, m100, nrow = 1)
# tmap_save(tm = both,
#           filename = here("graphs", paste0(Sys.Date(), "_homotypic-vcd-graphs.png")),
#           width = 8.5,
#           height = 7.5,
#           units = "in")

####################
# 300m
####################
# Exposed VCDs
t <- 30
r <- 300
# Time matrix
vcds_time.bin = 1*(time.mat <= t & time.mat > 0)
# Distance matrix
vcds_dist.bin = 1*(dist.mat < r)

sum(vcds_dist.bin, na.rm = TRUE)
sum(vcds_time.bin, na.rm = TRUE)

sum(vcds_dist.bin*vcds_time.bin, na.rm = TRUE)

ones_row <- apply(vcds_dist.bin*vcds_time.bin*sero.mat, 1, function(x){max(x,na.rm = TRUE)})
ones_col <- apply(vcds_dist.bin*vcds_time.bin*sero.mat, 2, function(x){max(x,na.rm = TRUE)})
temp <- names(ones_row)[which(ones_row == 1)]
temp_2 <- names(ones_col[which(ones_col == 1)])
m300 <- tm_shape(Yogya_Adm_1) + 
  tm_polygons(col = "Intervention",
              alpha = 0.75,
              palette = c("lightgray", "white")) + 
  tm_text(text = "New_Clustr") +
  tm_shape(subset(CRTND_VCD_data_SPDF, participant_id %in% unique(c(temp, temp_2)))) +
  tm_symbols(col = "serotype",
             shape = "serotype",
             legend.col.show = FALSE,
             legend.shape.show = FALSE,
             palette = df_col_sero$sero.color,
             shapes = df_col_sero$sero.shape,
             size = 0.2) + 
  tm_add_legend("symbol",
                col = df_col_sero$sero.color[1:4],
                shape = df_col_sero$sero.shape[1:4],
                labels = df_col_sero$serotype[1:4],
                size = 0.8) +
  tm_scale_bar(breaks = c(0,.5,1), text.size = 0.9, position = c("left", "top")) +
  tm_layout(legend.outside = TRUE,
            legend.outside.position = "right",
            legend.outside.size = 0.25,
            legend.text.size = 1)

tmap_save(tm = m300,
          filename = here("graphs", paste0(Sys.Date(), "_homotypic-vcd-graphs-300.png")),
          width = 7,
          height = 8.85,
          units = "in")
