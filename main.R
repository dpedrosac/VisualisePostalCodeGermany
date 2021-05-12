# This is code to plot data on a german map;
# Code developed by David Pedrosa, adapted from:
# https://arilamstein.com/blog/2016/05/16/use-case-study-mapping-german-zip-codes/

# General code
wdir <- "/media/storage/skripte/visPLZ/"
setwd(wdir)

# subfunctions 
ipak <- function(pkg){ # taken from https://gist.github.com/stevenworthington/3178163
new.pkg <- pkg[!(pkg %in% installed.packages()[, "Package"])]
if (length(new.pkg)) 
    install.packages(new.pkg, dependencies = TRUE)
sapply(pkg, require, character.only = TRUE)
}

## First specify the packages of interest and install them if necessary
packages <- c("maptools", "dplyr", "Rcpp", "RColorBrewer", "osmar", "XML", "sf", "ggmap", "tidyverse", "viridis", "stringr", 
"spdep", "proj4", "choroplethr", "ggplot2", "R6", "readr", "gpclib", "rgdal", "rgeos")
ipak(packages)

ger_plz <- readOGR(dsn = ".", layer = "plz-gebiete")
gpclibPermit()

#convert the raw data to a data.frame as ggplot works on data.frames
ger_plz@data$id <- rownames(ger_plz@data)
ger_plz.point <- fortify(ger_plz, region="id")
ger_plz.df <- inner_join(ger_plz.point,ger_plz@data, by="id")

#data file
df <- read.table(paste0(wdir, "de_plz_einwohner.csv"), sep =';', header=TRUE, colClasses = c("character", "integer"))

# variable name ‘region’ is needed for choroplethr
ger_plz.df$region <- ger_plz.df$plz
head(ger_plz.df)

# subclass choroplethr to make a class for your my need [sic]
GERPLZChoropleth <- R6Class("GERPLZChoropleth",
    inherit = choroplethr:::Choropleth,
    public = list(
        initialize = function(user.df) 
		{
            super$initialize(ger_plz.df, user.df)
        }
    )
)

#choropleth needs these two columnames - 'region' and 'value'
colnames(df) <- c("region", "value")

#instantiate new class with data
c <- GERPLZChoropleth$new(df)

# Plot data
c$ggplot_polygon = geom_polygon(aes(fill = value), color = NA)
c$legend <- "Number of Inhabitants"
c$set_num_colors(1)
fig <- c$render()
fig + scale_fill_gradientn(colours=brewer.pal(7,"Blues"), na.value = "transparent")

# 
hosp.raw <- read.table(paste0(wdir, "liste_kh.csv"), sep =';', header=TRUE, 
			colClasses = c("factor", "factor", "factor", "character", "character", "character", "factor", rep("integer", 51)))
hosp.raw[ hosp.raw == "-" ] <- 0
idx1 <- match(hosp.raw$PLZ, df$region)
idx1 <- subset(idx1, !is.na(idx1))
idx2 <- unique(idx1)
dhosp <- df
dhosp[,2] <- rep(0,dim(df)[1])

for (i in 1:length(idx2)){
idx_temp <- which(idx1 == idx2[i])
dhosp[idx2[i],2] <- sum(hosp.raw$Betten_Ins[idx_temp])
}


#choropleth needs these two columnames - 'region' and 'value'
colnames(dhosp) = c("region", "value")

#instantiate new class with data

c <- GERPLZChoropleth$new(dhosp)

#plot the data
c$ggplot_polygon = geom_polygon(aes(fill = value), color = NA)
c$legend= "Number of Inhabitants"
c$set_num_colors(1)
fig = c$render()
fig + scale_fill_gradientn(colours=brewer.pal(7,"Blues"), na.value = "transparent")


## Prepare data
# Define directories and files to read shapes from
shapefile 	<- "post_pl.shp"
my_plz 	<- "DE.txt"
datfile 		<- "results_plz.csv"
datfile_pop <- "bevoelkerung.csv"
de_plz 		<- paste0(wdir, my_plz)
dat_plz		<- paste0(wdir, datfile)
dat_pop	<- paste0(wdir, datfile_pop)

# Read the german postal codes and change the table according to the information  
plz_df <- read_tsv(file = de_plz, col_names = FALSE)  # read german postcodes as tab separated data
plz_df <- plz_df %>% rename(country_code = X1,
         postal_code = X2,
         place_name = X3,
         state = X4,
         state_code = X5,
         county = X6,
         county_code = X7,
         community = X8,
         community_code = X9,
         lat = X10,
         long = X11,
         accuracy = X12)  # accuracy of lat/lng from 1=estimated to 6=centroid

glimpse(plz_df) # sanity check, if one wants to see a part of the information
write.csv(plz_df, file = paste0(wdir, "data.csv"))

# Read the shape information into the workspace
shape_dest <- paste0(wdir,shapefile) # define the name of the file
de_shape <- sf::st_read(shape_dest) # read shapes
de_shape$PLZORT99 <- as.character(de_shape$PLZORT99) # in order to avoid conflict with special german characters ("Umlaute"), the next few lines are needed
Encoding(de_shape$PLZORT99) <- "latin1"
slice(de_shape$PLZORT99, 90:100)

# Create information of interest, i.e. add data to shapefile
de_shape$val <- 0
mydata = read.csv(dat_plz, head = FALSE, sep=";", colClasses='character') # number of patients according to postal code

for (i in 1:dim(mydata)[1]){
	idx <- match(mydata[i,1], as.factor(de_shape$PLZ99))
	de_shape$val[idx] = mydata[i,2]
  }

# Smoothing
#shapetot  <- readShapePoly(paste0(wdir, shapefile))
shapetot  <- rgdal::readOGR(paste0(wdir, shapefile))
proj4string(shapetot) <- CRS("+proj=longlat")
shapetemp <-spTransform(shapetot, CRS("+init=epsg:3857"))
#shapesmooth <- as(shapetemp, "data.frame")
#shapesmooth[is.na(shapesmooth)] <-0
coords <- coordinates(shapetemp)

IDs<-row.names(as(shapetemp, "data.frame"))
nsmooth <- 4
knn50 <- knn2nb(knearneigh(coords, k = nsmooth), row.names = IDs)
knn50 <- include.self(knn50)
  
# Creating the localG statistic for each of counties, with a k-nearest neighbor value of 5, and round this to 3 decimal places
smoothval <- localG(x = as.numeric(de_shape$val), listw = nb2listw(knn50, style = "B"), zero.policy = TRUE)
smoothval <- round(smoothval,3)
de_shape$sval <- log10(abs(smoothval))

## Start with the plotting of the data
# first plot (Germany)
nsteps <- 9
my_col <- sample(1:nsteps, length(de_shape), replace = T) # draws a random number between 1:7 for the length of (de_shape)
plot(de_shape[c("val", "geometry")],
     col = brewer.pal(nsteps, "Greens")[my_col], border = F)
	
# Second plot (regions)
# include a string which makes the selection of different states possible
#plz_reg <- plz_df %>% filter(state == "Hessen" | state == "Nordrhein-Westfalen" | state == "Rheinland-Pfalz")
#plz_reg <- plz_df %>% filter(state == "Hessen")
plz_reg <- plz_df
plz_reg_vec <- plz_reg$postal_code 						# returns the postal codes of interest for the selected states
my_rows <- de_shape$PLZ99 %in% plz_reg_vec				# returns a list of rows that correspond to the data of interest
my_col <- sample(1:nsteps, length(de_shape), replace = T)
has_data <- de_shape[my_rows, c("PLZ99", "val", "sval", "geometry")]
plot(has_data,
     col = brewer.pal(nsteps, "Blues")[my_col], border = F)
has_data$PLZ6 <- str_extract(has_data$PLZ99, "\\d\\d")

ggplot(data = has_data) +
  geom_sf(aes(fill = sval)) +
	scale_fill_gradient(low = "#FFFFFF", high = "#D0442D",
  space = "Lab", na.value = "grey50", guide = "colourbar",
  aesthetics = "fill") + #scale_fill_gradientn(colours = cividis(256)) + # scale_fill_viridis_c() +
  guides(fill = FALSE) -> p_region
p_region

data.all = read.csv(paste0(wdir, "data_all.csv"), head = TRUE, sep=";", colClasses='character')
# Create information of interest, i.e. add data to shapefile
de_shape$pop <- 0

for (i in 1:dim(data.all)[1]){
	idx <- match(as.factor(data.all$postal_code)[i], as.factor(de_shape$PLZ99))
	de_shape$pop[idx] = as.numeric(data.all$Var14[i])
  }

# Second plot (regions)
# include a string which makes the selection of different states possible
#plz_reg <- plz_df %>% filter(state == "Hessen" | state == "Nordrhein-Westfalen" | state == "Rheinland-Pfalz")

#plz_reg <- plz_df %>% filter(state == "Hessen")
plz_reg <- plz_df
plz_reg_vec <- plz_reg$postal_code 						# returns the postal codes of interest for the selected states
my_rows <- de_shape$PLZ99 %in% plz_reg_vec				# returns a list of rows that correspond to the data of interes
my_col <- sample(1:nsteps, length(de_shape), replace = T)
has_data <- de_shape[my_rows, c("PLZ99", "val", "sval", "geometry", "pop")]
has_data$PLZ6 <- str_extract(has_data$PLZ99, "\\d\\d")
has_data$pop = as.numeric(has_data$pop)

ggplot(data = has_data) +
  geom_sf(aes(fill = pop)) +
	scale_fill_gradient(low = "#FFFFFF", high = "#D0442D",
  space = "Lab", na.value = "#FFFFFF", guide = "colourbar",
  aesthetics = "fill") + #scale_fill_gradientn(colours = cividis(256)) + # scale_fill_viridis_c() +
  guides(fill = FALSE) -> p_region
p_region
