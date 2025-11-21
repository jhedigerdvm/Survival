
#and let's load all the libraries we need
library(terra)
library(ggplot2)
library(rnaturalearth)
library(tidyterra)
library(tidyverse)
library(viridis)
library(mapview)
library(ggspatial)
library(cowplot) #for insets
library(sf)
library(maptiles) #allows to bring in esri basemap
library(here)

#read in the shapefile with vect()
faith_shp <- vect("./raw/faith_shp/Pastures-FR NAD83TXSPSC.shp")
tx <- vect("./raw/shp_files/State.shp")
counties <- vect("./raw/tx_counties/County_Boundaries.shp")
dmp <- vect("./raw/faith_shp/DMP Deer-FR.shp")
usa <-  vect("./raw/shp_files/tl_2012_us_state.shp")
faithperim <- vect("./raw/faith_shp/Faith_Boundary.shp")


#get coord system
tx <- project(tx, "EPSG:6588")
counties <- project(counties, "EPSG:6588")
usa <- project(usa, "EPSG:6588")
# 
# ext(faith_shp)
# extfaith <- ext(1569918.3390678, 1687450.12404436, 1335888.2509673, 13298730.1967523)
# box_faith <- vect(extfaith, crs="EPSG:6588") #and make this extent a spatVector
# 
# ext_sotex <-ext(1490000, 1750000 , 13159888.2509673, 13358730.1967523)
# box_sotex <- vect(ext_sotex, crs="EPSG:6588") #
# mapview(box_sotex)

# unique(dmp$PASTURE)
dmp <- dmp[dmp$PASTURE %in% "W Yana", ]
usa <- usa[!usa$NAME %in% c("Alaska", "Hawaii", "American Samoa",
                                 "Commonwealth of the Northern Mariana Islands",
                                 "Guam", "United States Virgin Islands",
                                 "Puerto Rico")]
counties <- counties[counties$CNTY_NM %in% c("Dimmit", "Webb")]
mapview(counties)

yanas <- faith_shp[faith_shp$PASTURE %in% c("East Yana Pasture", "West Yana Pasture")]

mapview(yanas)


# Convert SpatVector to sf so we can create ESRI basemap
yanas <- st_as_sf(yanas) #convert to dataframe
dmp <- st_as_sf(dmp) #convert to dataframe
dmp <- dmp[,-4] #remove unnec column so we can rbind
dmp$PASTURE <- 'dmp' #rename pasture to dmp
dmp <- st_transform(dmp, st_crs(yanas)) #need to be on same CRS to bind
yanadmp <- rbind(yanas, dmp) #bind yana polygons with dmp polygons
tx <- st_as_sf(tx) 
faithperim <- st_as_sf(faithperim)

# crop <- crop(basemap, yanadmp)



#let's plot this SpatVector
usaplot<-ggplot() +
  geom_spatvector(data = usa_cont, size = 0.5, color = "black", fill = "transparent") +
  theme(
    panel.background = element_rect(fill = NA, color = NA),
    plot.background  = element_rect(fill = NA, color = NA),
    panel.grid = element_blank(),
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    axis.title = element_blank()
  )
theme_bw()

usaplot
ggsave('./figures/usaplot.png', usaplot, width = 10, height = 5)

###############33
#################
#################

tx_plot<-
  ggplot() + 
  geom_spatvector(data = tx, size = 0.5, color = "black", fill = "white") +
  geom_spatvector(data = counties, size = 0.5, color = "black", fill = "red") + #just dimmitt and webb
  # geom_spatvector(data = box_sotex, size = 0.5, color = "red", fill = "transparent")+
  #geom_spatvector(data = faithperim, size = 0.5, color = "red", fill = "red") + 
  
  theme(
    panel.background = element_rect(fill = NA, color = NA),
    plot.background  = element_rect(fill = NA, color = NA),
    panel.grid = element_blank(),
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    axis.title = element_blank()
  )

tx_plot

ggsave('./figures/txplot.png', tx_plot, width = 10, height = 10)

################
##############
###############


county_plot<-
  ggplot() + 
  geom_spatvector(data = counties, size = 0.5, color = "black", fill = "white") + #just dimmitt and webb
  geom_spatvector(data = faithperim, size = 0.5, color = "red", fill = "red") +
  # geom_spatvector(data = box_sotex, size = 0.5, color = "red", fill = "transparent")+
  #geom_spatvector(data = faithperim, size = 0.5, color = "red", fill = "red") + 
  
  theme(
    panel.background = element_rect(fill = NA, color = NA),
    plot.background  = element_rect(fill = NA, color = NA),
    panel.grid = element_blank(),
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    axis.title = element_blank()
  )

county_plot

ggsave('./figures/txplot.png', tx_plot, width = 10, height = 10)

################
##############
###############

# Get basemap (ESRI satellite)
basemap <- get_tiles(faithperim, provider = "Esri.WorldImagery", zoom = 15) #adjust zoom to 17 for high res


# Convert to raster for plotting
faith<-ggplot() +
  layer_spatial(basemap) +
  geom_sf(data = faithperim, color = 'red', linewidth = 1, fill = NA, size = 6) +
  geom_sf(data = yanadmp, fill = 'red' , size = 6) +
  # scale_color_viridis_d(option = 'plasma',
  #                       labels = c("DMP", "East Yana", "West Yana")
  # )+
  coord_sf(xlim = c(380737.484727737, 412310.857041923), #st_bbox(yanadmp)[c("xmin", "xmax")],
           ylim = c(3115000.51580854, 3130120.60288207), #st_bbox(yanadmp)[c("ymin", "ymax")],
           expand = FALSE) +
  theme(
    panel.background = element_rect(fill = NA, color = NA),
    plot.background  = element_rect(fill = NA, color = NA),
    legend.background = element_rect(fill = NA, color = NA),
    panel.grid = element_blank(),
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    axis.title = element_blank(),
    legend.title = element_blank()
    
  )+
  annotation_scale(style="bar",  #alternative style is "ticks"
                   height=unit(.5,"cm"),
                   pad_x = unit(15, "cm"),
                   pad_y = unit(1, "cm"))+
  #Let's add the north arrow - adjust height and padding for your individual screen
  annotation_north_arrow(
    height=unit(1.5, "cm"),
    width=unit(1.2, "cm"),
    pad_x = unit(10, "cm"),
    pad_y = unit(1, "cm"))

faith


inset<-
  ggdraw() +
  draw_plot(faith) +
  draw_plot(tx_plot,
            height = 0.5,
            x = -.35, #you will have to play with these values a bit to get them right!
            y = 0.5)+
  draw_plot(county_plot,
            height = 0.5,
            x = -.1, #you will have to play with these values a bit to get them right!
            y = 0.5)


inset

ggsave('./figures/faithinset.png', inset, width = 10, height = 5)





# faith_plot<-
#   ggplot() + 
#   geom_spatvector(data = counties, size = 0.5, color = "black", fill = "transparent") +
#   geom_spatvector(data = faithperim, size = 0.5, color = "black", fill = "red") + 
#   
#   # geom_spatvector(data = yanas, size = 0.5, color = 'black', fill = 'transparent') + 
#   # geom_spatvector(data = dmp, size = 0.5, color = 'black', fill = 'transparent') +
#   # geom_spatvector(data = feeders, size = 0.5, color = "black", fill = "black") + 
# 
#   # coord_sf(
#   #   xlim = c(ext_sotex$xmin, ext_sotex$xmax),
#   #   ylim = c(ext_sotex$ymin, ext_sotex$ymax),
#   #   expand = FALSE
#   # )+
#   theme(
#     panel.background = element_rect(fill = NA, color = NA),
#     plot.background  = element_rect(fill = NA, color = NA),
#     panel.grid = element_blank(),
#     axis.text = element_blank(),
#     axis.ticks = element_blank(),
#     axis.title = element_blank(),
#     
#   )
# 
# faith_plot

inset<-ggdraw() +
  draw_plot(faith_plot) +
  draw_plot(tx_plot,
            height = 0.4,
            x = 0.15, #you will have to play with these values a bit to get them right!
            y = 0.55)


ggsave('./figures/inset.png', inset, width = 5, height = 10)


################3
#################

# Convert to raster for plotting
basemap <- get_tiles(yanadmp, provider = "Esri.WorldImagery", zoom = 15) #adjust zoom to 17 for high res

yana<-ggplot() +
  layer_spatial(basemap) +
  geom_sf(data = yanadmp, aes(color = PASTURE), linewidth = 1, fill = NA, size = 6) +
  scale_color_viridis_d(option = 'plasma',
                        labels = c("DMP", "East Yana", "West Yana")
  )+
  coord_sf(xlim = c(1658000.33357406, 1675913.06965823), #st_bbox(yanadmp)[c("xmin", "xmax")],
           ylim = c(13260966.6762228, 13270751.0856088), #st_bbox(yanadmp)[c("ymin", "ymax")],
           expand = FALSE) +
  theme(
    panel.background = element_rect(fill = NA, color = NA),
    plot.background  = element_rect(fill = NA, color = NA),
    legend.background = element_rect(fill = NA, color = NA),
    panel.grid = element_blank(),
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    axis.title = element_blank(),
    legend.title = element_blank()
    
  )+
  annotation_scale(style="bar",  #alternative style is "ticks"
                   height=unit(.5,"cm"),
                   pad_x = unit(4, "cm"),
                   pad_y = unit(1, "cm"))+
  #Let's add the north arrow - adjust height and padding for your individual screen
  annotation_north_arrow(
    height=unit(1.5, "cm"),
    width=unit(1.2, "cm"),
    pad_x = unit(1, "cm"),
    pad_y = unit(1, "cm"))


yana

ggsave('./figures/yana.png', yana, width = 10, height = 5)
