# ----HEADER------------------------------------------------------------------------------------------------------------
# Author: JF
# Date: 09/05/2018
# Purpose: Run custom function to create U5M maps
# source("/homes/jfrostad/_code/lbd/u5m_mapping/plot.R", echo=T)
#***********************************************************************************************************************

# ----CONFIG------------------------------------------------------------------------------------------------------------
# clear memory
rm(list=ls())

# runtime configuration
if (Sys.info()["sysname"] == "Linux") {
  j_root <- "/home/j/"
  h_root <- "/homes/jfrostad/"
  arg <- commandArgs()[-(1:3)] # First args are for unix use only
  
  if (length(arg)==0) {
    # arg <- c("IND", #current project iteration
    #          "8", #output version
    #          1) #number of cores provided to multicore functions
  }
  

  # necessary to set this option in order to read in a non-english character shapefile on a linux system (cluster)
  Sys.setlocale(category = "LC_ALL", locale = "C")
  
} else {
  j_root <- "J:"
  h_root <- "H:"
  # arg <- c("IND", #current project iteration
  #          "4", #output version
  #          1) #number of cores provided to multicore functions
}

#use cairo to render instead of quartz (quartz causes big slowdowns with geom_sf)
if(!identical(getOption("bitmapType"), "cairo") && isTRUE(capabilities()[["cairo"]])){
  options(bitmapType = "cairo")
}

## Set core_repo location and indicator group
user            <- Sys.info()['user']
core_repo       <- '/homes/jfrostad/_code/lbd/lbd_core/'
commondir       <- paste(core_repo, 'mbg_central/share_scripts/common_inputs', sep = '/')
package_list    <- c(t(read.csv(paste(commondir, 'package_list.csv', sep = '/'), header = FALSE)))

#load packages
package_lib    <- sprintf('%s_code/_lib/pkg',h_root)
## Load libraries and  MBG project functions.
.libPaths(package_lib)
pacman::p_load(data.table, RMySQL, feather, ggplot2, gridExtra, maptools, 
               parallel, raster, RColorBrewer, rgdal, rgeos, scales, survival, sf, tictoc, tidyverse, viridis) 

# Use setup.R functions to load common LBD packages and mbg_central "function" scripts
message('Loading in required R packages and MBG functions')
source(paste0(core_repo, '/mbg_central/setup.R'))
mbg_setup(package_list = package_list, repos = core_repo)

#capture date
today <- Sys.Date() %>% gsub("-", "_", .)

#options
run_date <- '2018_08_14_14_43_37'
indicator_group <- 'u5m'
indicator <- 'died_under5'
type <- 'mean'
raked <- T
start_year <- 2000
end_year <- 2017
cores <- 10
#***********************************************************************************************************************

# ----IN/OUT------------------------------------------------------------------------------------------------------------
###Input###
#raw data
data.dir <- file.path('/share/geospatial/mbg/u5m/died_under5/output/', run_date)
u5m.paths <- data.table(raster=file.path(data.dir, 'died_under5_mean_raked_2000_2017.tif'),
                        admin1=file.path(data.dir, 'died_under5_mean_raked_ad1.csv'),
                        admin2=file.path(data.dir, 'died_under5_mean_raked_ad2.csv'))

###Output###
out.dir  <- file.path(j_root, 'temp/jfrostad/u5m_mapping/')
#***********************************************************************************************************************

# ---FUNCTIONS----------------------------------------------------------------------------------------------------------
##function lib##
#map functions#
hap.function.dir <- file.path(h_root, '_code/lbd/u5m_mapping/_lib')
#this pulls hap collapse helper functions
file.path(hap.function.dir, '/map_fx.R') %>% source
#***********************************************************************************************************************

# ---PREP DATA----------------------------------------------------------------------------------------------------------
#read in the proper annotations (borders, lakes, mask)
annotations <- load_map_annotations()

#read in input data and prepare it for mapping
tic('total')
tic('loading data')
data <- load_map_results(indicator, indicator_group, run_date, type, raked, start_year, end_year, single_year=2017,
                         custom_path = u5m.paths,
                         use.sf = F,
                         #geo_levels='admin2',
                         cores=cores)
toc()

#define extent of map
zoom.afr <- data.table(x1=-10, x2=50, y1=-20, y2=40)
zoom.global <- data.table(x1=-120, x2=150, y1=-40, y2=55)
#***********************************************************************************************************************

# ---CREATE MAPS--------------------------------------------------------------------------------------------------------
#render and save maps
colors <- c('#ffffe0','#ffe4ac','#ffc879','#ffa84c','#ff8725','#ff5c03','#f12861','#cc117d','#a60383','#800080')
colors <- magma(10, direction=-1)
color_values <- c(seq(0, .05, length.out = 4), seq(.05, .10, length.out = 4), seq(.10, .25, length.out = 4)) %>%
  unique %>%
  rescale

tic('plotting')
gg <- plot_map(data$admin1, annotations, limits=c(0, .25), title='', 
               legend_colors=colors, legend_color_values=color_values,
               legend_breaks=seq(0, .25, .05), legend_labels=c(seq(0, .20, .05), ".25+"),
               legend_title='under-5 mortality', custom_scale=T,
               pop.mask=F, lake.mask=T, stage3.mask=T, borders=T,
               zoom=zoom.global)
toc()

tic('ggsaving')
ggsave(filename=file.path(out.dir, 'test.png'), plot=gg, 
       width=10, height=6, units='in', dpi=300)
toc()

tic('ggsaving600')
ggsave(filename=file.path(out.dir, 'test600_nosf.png'), plot=gg, 
       width=10, height=6, units='in', dpi=600)
toc()

tic('png saving')
png(filename=file.path(out.dir, 'testreg_nosf.png'), 
    units='px', width=2700, height=1500)
print(gg)
dev.off()
toc()

tic('bmp saving')
png(filename=file.path(out.dir, 'testreg_nosf.bmp'), 
    units='px', width=2700, height=1500)
print(gg)
dev.off()
toc()

toc()
#***********************************************************************************************************************

# ---SCRAPS-------------------------------------------------------------------------------------------------------------