
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
feeders <- vect("./raw/faith_shp/Feeder Protein-FR.shp")
water <- vect("./raw/faith_shp/Water Well-FR.shp")
waterline <- vect("./raw/faith_shp/Water Line-FR.shp")
#get coord system
tx <- project(tx, "EPSG:6588")
counties <- project(counties, "EPSG:6588")
usa <- project(usa, "EPSG:6588")
# 

mapview(faith_shp)
faith_shp$PASTURE


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
counties <- counties[counties$CNTY_NM %in% c("Dimmit","Webb")]#, 
mapview(counties)

yanas <- faith_shp[faith_shp$PASTURE %in% c("East Yana Pasture", "West Yana Pasture")]

mapview(yanas)
mapview(waterline)

# Convert SpatVector to sf so we can create ESRI basemap
yanas <- st_as_sf(yanas) #convert to dataframe
dmp <- st_as_sf(dmp) #convert to dataframe
dmp <- dmp[,-4] #remove unnec column so we can rbind
dmp$PASTURE <- 'dmp' #rename pasture to dmp
dmp <- st_transform(dmp, st_crs(yanas)) #need to be on same CRS to bind
yanadmp <- rbind(yanas, dmp) #bind yana polygons with dmp polygons
tx <- st_as_sf(tx) 
faithperim <- st_as_sf(faithperim)
feeders <- st_as_sf(feeders)


feeders <- feeders %>% filter(PASTURE %in% c("E Yana", "W Yana", "DMP Pen 1",  "DMP Pen 2", "DMP Pen 3","DMP Pen 4"))
unique(feeders$PASTURE)

feeders1 <- feeders  %>% filter(PASTURE %in% c("E Yana", "W Yana", "DMP Pen 3","DMP Pen 4"))


#let's plot this SpatVector
usaplot<-ggplot() +
  geom_spatvector(data = usa, size = 0.5, color = "black", fill = "white") +
  geom_spatvector(
    data = usa[usa$NAME == "Texas", ],
    fill = "darkgrey",
    color = "black"
  )+
  theme(
    panel.background = element_rect(fill = NA, color = NA),
    plot.background  = element_rect(fill = NA, color = NA),
    panel.grid = element_blank(),
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    axis.title = element_blank()
  )

usaplot
# ggsave('./figures/usaplot.png', usaplot, width = 10, height = 5)

###############33
#################
#################

tx_plot<-
  ggplot() + 
  geom_spatvector(data = tx, size = 0.5, color = "black", fill = "white") +
  geom_spatvector(data = counties[counties$CNTY_NM == "Dimmit"], size = 0.5, color = "black", fill = "darkgrey") + #just dimmitt and webb
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

# ggsave('./figures/txplot.png', tx_plot, width = 10, height = 10)

################
##############
###############

faith_shp

county_plot<-
  ggplot() + 
  geom_spatvector(data = counties[counties$CNTY_NM == "Dimmit"], size = 0.5, color = "black", fill = "white") + #just dimmitt and webb
  geom_spatvector(data =yanadmp, size = 0.5, color = "black", fill = "darkgrey") +
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

# ggsave('./figures/txplot.png', tx_plot, width = 10, height = 10)

################
##############
###############
# 
# # Get basemap (ESRI satellite)
# basemap <- get_tiles(faithperim, provider = "Esri.WorldImagery", zoom = 12) #adjust zoom to 17 for high res
# basemap_masked <- mask(basemap, vect(faithperim))
# 
# 
# # Convert to raster for plotting
# faith<-ggplot() +
#   layer_spatial(basemap_masked) +
#   # geom_spatvector(data = faithperim, size = 0.5, color = "darkgrey", fill = "transparent") +
#   geom_sf(data = faithperim, color = 'white', linewidth = 1, fill = NA, size = 6) +
#   geom_sf(data = yanadmp , aes(color = PASTURE), fill = NA, linewidth = 1,  size = 6) +
#   scale_color_viridis_d(option = 'plasma',
#                         labels = c("DMP", "East Yana", "West Yana"))+
#   # scale_color_viridis_d(option = 'plasma',
#   #                       labels = c("DMP", "East Yana", "West Yana")
#   # )+
#   # coord_sf(xlim = c(380737.484727737, 412310.857041923), #st_bbox(yanadmp)[c("xmin", "xmax")],
#   #          ylim = c(3115000.51580854, 3130120.60288207), #st_bbox(yanadmp)[c("ymin", "ymax")],
#   #          expand = FALSE) +
#   theme(
#     panel.background = element_blank(),# element_rect(fill = NA, color = NA),
#     plot.background  = element_blank(),#, element_rect(fill = NA, color = NA),
#     legend.background = element_blank(), # element_rect(fill = NA, color = NA),
#     panel.grid = element_blank(),
#     axis.text = element_blank(),
#     axis.ticks = element_blank(),
#     axis.title = element_blank(),
#     legend.position = "none"
#     # legend.title = element_blank(),
#     # legend.text = element_blank(),
#     # legend.box = element_blank(),
#     # legend.key = element_blank(),
#     # legend
#     
#   )#+
#   # annotation_scale(style="bar",  #alternative style is "ticks"
#   #                  height=unit(.5,"cm"),
#   #                  pad_x = unit(15, "cm"),
#   #                  pad_y = unit(1, "cm"))+
#   # #Let's add the north arrow - adjust height and padding for your individual screen
#   # annotation_north_arrow(
#   #   height=unit(1.5, "cm"),
#   #   width=unit(1.2, "cm"),
#   #   pad_x = unit(10, "cm"),
#   #   pad_y = unit(1, "cm"))
# 
# faith
# 
# 
# inset<-
#   ggdraw() +
#   draw_plot(faith) +
#   draw_plot(tx_plot,
#             height = 0.5,
#             x = -.35, #you will have to play with these values a bit to get them right!
#             y = 0.5)+
#   draw_plot(county_plot,
#             height = 0.5,
#             x = -.1, #you will have to play with these values a bit to get them right!
#             y = 0.5)
# 
# 
# inset
# 
# ggsave('./figures/faithinset.png', inset, width = 10, height = 5)
# 
# 
# 
# 
# 
# # faith_plot<-
# #   ggplot() + 
# #   geom_spatvector(data = counties, size = 0.5, color = "black", fill = "transparent") +
# #   geom_spatvector(data = faithperim, size = 0.5, color = "black", fill = "red") + 
# #   
# #   # geom_spatvector(data = yanas, size = 0.5, color = 'black', fill = 'transparent') + 
# #   # geom_spatvector(data = dmp, size = 0.5, color = 'black', fill = 'transparent') +
# #   # geom_spatvector(data = feeders, size = 0.5, color = "black", fill = "black") + 
# # 
# #   # coord_sf(
# #   #   xlim = c(ext_sotex$xmin, ext_sotex$xmax),
# #   #   ylim = c(ext_sotex$ymin, ext_sotex$ymax),
# #   #   expand = FALSE
# #   # )+
# #   theme(
# #     panel.background = element_rect(fill = NA, color = NA),
# #     plot.background  = element_rect(fill = NA, color = NA),
# #     panel.grid = element_blank(),
# #     axis.text = element_blank(),
# #     axis.ticks = element_blank(),
# #     axis.title = element_blank(),
# #     
# #   )
# # 
# # faith_plot
# 
# inset<-ggdraw() +
#   draw_plot(faith_plot) +
#   draw_plot(tx_plot,
#             height = 0.4,
#             x = 0.15, #you will have to play with these values a bit to get them right!
#             y = 0.55)
# 
# 
# ggsave('./figures/inset.png', inset, width = 5, height = 10)


################3
#################

# Convert to raster for plotting
basemap <- get_tiles(yanadmp, provider = "Esri.WorldImagery", zoom = 15) #adjust zoom to 17 for high res
basemap_masked <- mask(basemap, vect(yanadmp))

yana<-ggplot() +
  layer_spatial(basemap_masked) +
  geom_sf(data = yanadmp, aes(color = PASTURE), linewidth = 1.5, fill = NA, size = 6) +
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
    legend.title = element_blank(),
    legend.position = "inside", 
    legend.position.inside =  c(.32, .23),
    legend.direction = 'horizontal',
    legend.text = element_text(size = 14),
    plot.margin = margin(t = 2.5, r = 0, b = 0, l = 0, unit = "cm")


    
    
  )+
  annotation_scale(style="bar",  #alternative style is "ticks"
                   height=unit(.5,"cm"),
                   pad_x = unit(4, "cm"),
                   pad_y = unit(1, "cm"),
                   text_cex = 1)+
  #Let's add the north arrow - adjust height and padding for your individual screen
  annotation_north_arrow(
    height=unit(1.5, "cm"),
    width=unit(1.2, "cm"),
    pad_x = unit(1, "cm"),
    pad_y = unit(1, "cm"))


yana


inset<-ggdraw() +
  draw_plot(yanafeeders) +
  draw_plot(usaplot,
            height = 0.25,
            x = -.2, #you will have to play with these values a bit to get them right!
            y = 0.75)+
  draw_plot(tx_plot,
            height = 0.25,
            x = 0, #you will have to play with these values a bit to get them right!
            y = 0.75)+
  draw_plot(county_plot,
            height = 0.25,
            x = .2, #you will have to play with these values a bit to get them right!
            y = 0.75)

inset


ggsave('./figures/insetfaith.jpg', inset, width = 12, height = 10)


# ggsave('./figures/yana.png', yana, width = 10, height = 5)


#faith ranch camera/water locations

# Convert to raster for plotting
basemap <- get_tiles(yanadmp, provider = "Esri.WorldImagery", zoom = 15) #adjust zoom to 17 for high res
basemap_masked <- mask(basemap, vect(yanadmp))

feeders1 <- feeders  %>% filter(PASTURE %in% c("E Yana", "W Yana", "DMP Pen 2", "DMP Pen 3"))


yanafeeders<-ggplot() +
  layer_spatial(basemap_masked) +
  geom_sf(data = yanadmp, aes(color = PASTURE), linewidth = 1.5, fill = NA, size = 6) +
  geom_spatvector(data = feeders1, shape = 21, size = 2, stroke = 1,  color = "black", fill = "white")+
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
    legend.title = element_blank(),
    legend.position = "inside",
    legend.position.inside =  c(.32, .23),
    legend.direction = 'horizontal',
    legend.text = element_text(size = 14),
    plot.margin = margin(t = 2.5, r = 0, b = 0, l = 0, unit = "cm")

  )+
  annotation_scale(style="bar",  #alternative style is "ticks"
                   height=unit(.5,"cm"),
                   pad_x = unit(4, "cm"),
                   pad_y = unit(1, "cm"),
                   text_cex = 1)+
  #Let's add the north arrow - adjust height and padding for your individual screen
  annotation_north_arrow(
    height=unit(1.5, "cm"),
    width=unit(1.2, "cm"),
    pad_x = unit(1, "cm"),
    pad_y = unit(1, "cm"))


yanafeeders

