library(tidyverse)
library(terra)
library(exactextractr)
library(tidyterra)
library(patchwork)
library(metR)

#load data
ET <- rast("evapotrans/evapotrans/")
ET_annual <- ET*12

miss <- vect("Mississippi_7030047060.geojson")
yang <- vect("Yangtze_4030009880.geojson")
euph <- vect("TigrisEuphrates_2030073570.geojson")

basins <- rbind(miss,yang,euph)

#extract mean, min, max values across each basin
vals <- extract(ET_annual,basins,fun=c("mean","min","max"))
vals$ID <- c("Mississippi","Yangtze","Tigris-Euphrates")

# > vals
# ID                    mean_evapotrans min_evapotrans max_evapotrans
# 1      Mississippi        556.2609       148.3644       962.9556
# 2          Yangtze        680.9409       108.4476      1001.2320
# 3 Tigris-Euphrates        225.4540        78.8496       395.7804


#plot for each basin (claude generated)
names(ET_annual) <- "ET"
basins$name <- c("Mississippi", "Yangtze", "Tigris-Euphrates")

plot_basin <- function(i) {
  b <- basins[i, ]
  df <- mask(crop(ET_annual, b, snap = "out"), b) |>
    as.data.frame(xy = TRUE, na.rm = TRUE)
  bdf <- as.data.frame(geom(b))   # polygon outline for ggplot
  
  ggplot() +
    geom_raster(data = df, aes(x, y, fill = ET)) +
    geom_polygon(data = bdf, aes(x, y, group = part),
                 fill = NA, color = "grey20", linewidth = 0.8) +
    scale_fill_viridis_c(option = "mako", direction = -1, name = "ET\n(mm/yr)") +
    coord_quickmap() +
    labs(title = b$name, x = NULL, y = NULL) +
    theme_minimal()
}

p <- lapply(seq_len(nrow(basins)), plot_basin)
p[[1]]
wrap_plots(p, nrow = 1)          # all three side by side


