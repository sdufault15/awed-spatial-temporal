library(tidyverse)
library(here)
library(spatstat)
library(tmap)
library(rgdal)
library(leaflet)
library(raster)

load(here("data","2020-12-08_work-ts-spat.RData"))
df_col_sero <- data.frame(serotype = c(paste0("DENV", 1:4), "Unknown", "Test-negative"),
                          sero.color = c("#0AA4D1", "#F3C73E", "#EF7822", "#76C1A8", "#003D58", "gray"),
                          sero.shape = c(21:24, 4, 25))
# Tidy Up for Mapping
work_ts_spat <- work_ts_spat %>% 
  mutate(serotype = case_when(serotype == "denv1" ~ "DENV1",
                              serotype == "denv2" ~ "DENV2",
                              serotype == "denv3" ~ "DENV3",
                              serotype == "denv4" ~ "DENV4",
                              serotype == "unk_serotype"~ "Unknown",
                              serotype == "test-negative control" ~ "Test-negative"),
         intervention = ifelse(intervention == 0, "Untreated", "Intervention")) %>% 
  full_join(df_col_sero) 

Yogya_Adm_1 <- readOGR("~/Box/2019-Research/WMP/CR-TND/awed-spatial-temporal/yogya-shape-files/", "Yk_CaseControl_20160512")
int_assignment <- work_ts_spat %>% 
  dplyr::select(cluster, Intervention = intervention) %>% 
  distinct() %>% 
  mutate(New_Clustr = ifelse(cluster < 10, paste0("0", cluster), as.character(cluster))) %>% 
  dplyr::select(-cluster)

tpf <- work_ts_spat %>% 
  group_by(cluster) %>% 
  mutate(`Test-positive fraction` = sum(dengue)/(sum(dengue) + sum(1-dengue))) %>% 
  dplyr::select(cluster, `Test-positive fraction`) %>% 
  ungroup() %>% 
  distinct()
work_ts_spat <- full_join(work_ts_spat, tpf)

tpf <- tpf %>% 
  mutate(New_Clustr = ifelse(cluster < 10, paste0("0", cluster), as.character(cluster))) %>% 
  dplyr::select(-cluster)

Yogya_Adm_1@data <- left_join(Yogya_Adm_1@data, int_assignment, by = "New_Clustr")
Yogya_Adm_1@data <- left_join(Yogya_Adm_1@data, tpf, by = "New_Clustr")

map_yogya <- 
  tm_shape(Yogya_Adm_1) +
  tm_borders()

# Creating Spatial Points Dataframe
CRTND_dengue_data_SPDF <- SpatialPointsDataFrame(coords = work_ts_spat[,c("longitude", "latitude")],
                                                 data = work_ts_spat[, c("cluster", "intervention", "enrolment_date", "illness_onset", "dengue", "age", "sex", "serotype", "sero.color", "sero.shape", "Test-positive fraction")],
                                                 proj4string = CRS("+init=epsg:4326")) # sets the projection to WGS 1984 using lat/long. Optional but good to specify

CRTND_VCD_data_SPDF <- subset(CRTND_dengue_data_SPDF, dengue == 1)

# Dengue point plot
map_dots <- map_yogya +
  tm_polygons(col = "Intervention",
              alpha = 0.5,
              palette = c("lightgray", "white")) +
  tm_shape(CRTND_VCD_data_SPDF) + 
  tm_symbols(col = "serotype",
          shape = "serotype",
          legend.col.show = FALSE,
          legend.shape.show = FALSE,
          palette = df_col_sero$sero.color,
          shapes = df_col_sero$sero.shape,
          size = 0.15) + 
  tm_add_legend("symbol", 
                col = df_col_sero$sero.color[1:5], 
                shape = df_col_sero$sero.shape[1:5],
                labels = df_col_sero$serotype[1:5]) +
  tm_layout(legend.text.size = 0.9,
            legend.title.size = 1,
            legend.title.color = "white",
            inner.margins = 0.1,
            title = "A",
            title.fontface = "bold") +
  tm_scale_bar(breaks = c(0,.50,1.00), text.size = 0.9, position = c("left", "top"))

# Dengue point plot
map_tpf <- map_yogya +
  tm_polygons(col = "Test-positive fraction") + #+
  # tm_dots(col = "serotype",
  #         palette = df_col_sero$sero.color,
  #         size = 0.15) + 
  tm_layout(legend.text.size = 0.9,
            legend.title.size = 1,
            # legend.title.color = "white",
            inner.margins = 0.1,
            title = "B",
            title.fontface = "bold") +
  tm_scale_bar(breaks = c(0,.50,1.00), text.size = 0.9, position = c("left", "top")) +
  tm_text(text = "New_Clustr",
          size = 0.8)


### Smoother
Yogya_Adm_1 <- spTransform(Yogya_Adm_1, crs(CRTND_dengue_data_SPDF))
yogya_Owin <- owin(xrange=range(work_ts_spat$longitude),yrange=range(work_ts_spat$latitude))
CaseControl_ppp <- ppp(work_ts_spat$longitude, work_ts_spat$latitude, 
                       window = yogya_Owin, 
                       marks=as.factor(work_ts_spat$dengue))

# Default smooth
risk_est <-  relrisk(CaseControl_ppp, relative = FALSE) 

risk_raster <- raster(risk_est, crs = crs(Yogya_Adm_1))
risk_raster_2 <- mask(risk_raster, Yogya_Adm_1)
p_rr <- tm_shape(Yogya_Adm_1) +
  tm_borders() +
  tm_shape(risk_raster_2) +
  tm_raster(title = "Test-positive fraction",
            style = "fixed",
            breaks = c(0,0.05,0.1,0.15,0.2)) +
  tm_shape(CRTND_VCD_data_SPDF) +
  tm_dots(#shape = "intervention",
    col = "black",
    shape = 3,
    alpha = 0.7,
    size = 0.05) +
  tm_shape(Yogya_Adm_1) +
  tm_text("Intervention",
          size = 0.7) +
  tm_borders() +
  tm_layout(#legend.outside = TRUE,
            # legend.outside.position = "right",
            # legend.outside.size = 0.25
            legend.text.size = 0.9,
            legend.title.size = 1,
            # legend.title.color = "white",
            inner.margins = 0.1,
            title = "C",
            title.fontface = "bold") +
  tm_scale_bar(breaks = c(0,.50,1.00), text.size = 0.9, position = c("left", "top")) 

all <- tmap_arrange(map_dots, map_tpf, p_rr)

tmap_save(tm = all,
          filename = here("graphs", paste0(Sys.Date(), "_overall-denv-tmap.png")),
          width = 17.5,
          height = 7.5,
          units = "in")
