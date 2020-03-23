---
title: "cbp v. salty coastlines"
output: 
  html_document:
    keep_md: true
---

### To get a raster of distance to coastline
from an [answer](https://stackoverflow.com/questions/35555709/global-raster-of-geographic-distances "Global Raster of Geographic Distances") on StackOverflow.


```r
library(raster)
```

```
## Warning: package 'raster' was built under R version 3.5.2
```

```
## Loading required package: sp
```

```
## Warning: package 'sp' was built under R version 3.5.2
```

```r
library(maptools)
```

```
## Warning: package 'maptools' was built under R version 3.5.2
```

```
## Checking rgeos availability: TRUE
```

```r
data(wrld_simpl)

# Create a raster template for rasterizing the polygons. 
# (set the desired grid resolution with res)
r <- raster(xmn=-180, xmx=180, ymn=-90, ymx=90, res=1)

# Rasterize and set land pixels to NA
r2 <- rasterize(wrld_simpl, r, 1)
r3 <- mask(is.na(r2), r2, maskvalue=1, updatevalue=NA)

# Calculate distance to nearest non-NA pixel
d <- distance(r3)

plot(d)
```

![](cbp_occurance_salty_coasts_files/figure-html/unnamed-chunk-1-1.png)<!-- -->

I got the soil salinity data here, from the [ISRIC Data Hub.](https://data.isric.org/geonetwork/srv/eng/catalog.search#/metadata/ca880bd4-cff8-11e9-8046-0cc47adaa92c "WoSIS snapshot: September 2019")

The _Capsella bursa-pastoris_ occurences were downloaded from iNaturalist.

### iNaturalist data and global soil chemical profiles

```r
# read inaturalist data
obs <- read.csv("~/Documents/Josephs_lab/worldwide_salt_analysis/cbp-observations-81425.csv")

# read soil profile data
profiles <- read.delim("~/Documents/Josephs_lab/worldwide_salt_analysis/wosis_201909_profiles.tsv")

# read chemical data
#chem <- read.delim("~/Documents/Josephs_lab/worldwide_salt_analysis/wosis_201909_layers_chemical.tsv")

# just want to see what all the columns are
#colnames(chem)

# only lat and long of the profiles atm
profiles_drop <- profiles[ , c("profile_id", "latitude", "longitude")]
```

#### My brain needs pictures


```r
library("ggplot2")
```

```
## Warning: package 'ggplot2' was built under R version 3.5.2
```

```r
theme_set(theme_bw())
library("sf")
```

```
## Warning: package 'sf' was built under R version 3.5.2
```

```
## Linking to GEOS 3.7.2, GDAL 2.4.2, PROJ 5.2.0
```

```r
library("rnaturalearth")
library("rnaturalearthdata")

world <- ne_countries(scale = "medium", returnclass = "sf")
class(world)
```

```
## [1] "sf"         "data.frame"
```

```r
#plot Capsella bursa-pastoris observations with salt data
ggplot(data = world) +
    geom_sf() + geom_point(data = obs, aes(x = longitude, y = latitude), size = 1.5, 
        shape = 23, fill = "darkred") 
```

```
## Warning: Removed 117 rows containing missing values (geom_point).
```

![](cbp_occurance_salty_coasts_files/figure-html/unnamed-chunk-3-1.png)<!-- -->

The same plot but with the soil points on top because I wanted to see the levelf of overlap between the two.

```r
wsalt <- ggplot(data = world) +
    geom_sf() + geom_point(data = obs, aes(x = longitude, y = latitude), size = 1.5, 
        shape = 23, fill = "darkred") + geom_point(data = profiles_drop, aes(x = longitude, y = latitude), size = 1, alpha = 0.3, fill = "skyblue4")

wsalt
```

```
## Warning: Removed 117 rows containing missing values (geom_point).
```

![](cbp_occurance_salty_coasts_files/figure-html/unnamed-chunk-4-1.png)<!-- -->

I abandonded this idea for the time being and decided I would try maxent.

### MAXENT

```r
# I think I need all of these packages?
library(raster)
library(rgdal)
```

```
## Warning: package 'rgdal' was built under R version 3.5.2
```

```
## rgdal: version: 1.4-7, (SVN revision 845)
##  Geospatial Data Abstraction Library extensions to R successfully loaded
##  Loaded GDAL runtime: GDAL 2.4.2, released 2019/06/28
##  Path to GDAL shared files: /Library/Frameworks/R.framework/Versions/3.5/Resources/library/sf/gdal
##  GDAL binary built with GEOS: FALSE 
##  Loaded PROJ.4 runtime: Rel. 5.2.0, September 15th, 2018, [PJ_VERSION: 520]
##  Path to PROJ.4 shared files: /Library/Frameworks/R.framework/Versions/3.5/Resources/library/sf/proj
##  Linking to sp version: 1.3-1
```

```r
library(maps)
library(mapdata)
library(dismo)  
library(lubridate)
```

```
## 
## Attaching package: 'lubridate'
```

```
## The following object is masked from 'package:base':
## 
##     date
```

```r
library(ggmap)
```

```
## Warning: package 'ggmap' was built under R version 3.5.2
```

```
## Google's Terms of Service: https://cloud.google.com/maps-platform/terms/.
```

```
## Please cite ggmap if you use it! See citation("ggmap") for details.
```

```
## 
## Attaching package: 'ggmap'
```

```
## The following object is masked from 'package:dismo':
## 
##     geocode
```

```r
library(ggplot2)
library(maptools)
library(dplyr)
```

```
## Warning: package 'dplyr' was built under R version 3.5.2
```

```
## 
## Attaching package: 'dplyr'
```

```
## The following objects are masked from 'package:lubridate':
## 
##     intersect, setdiff, union
```

```
## The following objects are masked from 'package:raster':
## 
##     intersect, select, union
```

```
## The following objects are masked from 'package:stats':
## 
##     filter, lag
```

```
## The following objects are masked from 'package:base':
## 
##     intersect, setdiff, setequal, union
```

```r
library(sf)
#library(rJava)
```

The shapefile of coastlines comes from [here.](https://www.naturalearthdata.com/downloads/10m-physical-vectors/10m-coastline/) (Also, I just realized i could probably get it directly from the rnaturalearth package.)

Loosely following [this tutorial.](https://cmerow.github.io/RDataScience/3_6_Teaching_Ecoinformatics.html) which is pretty much the same as the R documentation for MAXENT.


```r
# get rid of inaturalist observations without location information
cbp=subset(obs, !is.na(longitude) & !is.na(latitude))

# find and eliminate duplicate locations
cbpdups=duplicated(cbp[, c("longitude", "latitude")])
cbp_nodups <-cbp[!cbpdups, ]

# withold 20% of the data for testing the model
occr=cbind.data.frame(cbp_nodups$longitude,cbp_nodups$latitude)
fold <- kfold(occr, k=5)
cbptest <- occr[fold == 1, ]
cbptrain <- occr[fold != 1, ]
```

I end up not using this shapefile of coastlines (in the next code chunk) because the shapes are not polygons but linestrings and I was not sure what to do with it.

```r
# load coastline shape file
coasts <- st_read("~/Downloads/ne_10m_coastline/ne_10m_coastline.shp")

# view metadata
coasts

# plot the shape datafile as a sanity check.
ggplot() + 
  geom_sf(data = coasts, size = 3, color = "black", fill = "cyan1") +
  ggtitle("Coast Shape Data") + 
  coord_sf()
```
Yup, looks like coastlines.



```r
# fit MAXENT model to training set
# d is the raster of distance from coast created earlier
cbp.me <- maxent(d, cbptrain)
```
Above code chunk causes this error:

Error in match.names(clabs, names(xi)) : 
  names do not match previous names

There is an issue with making a MaxEnt model with only one predictor (distance to coastline, in this case.)
Workaround from StackOverflow [answer](https://stackoverflow.com/questions/47599192/maxent-error-in-match-namesclabs-namesxi-names-do-not-match-previous-nam).


```r
# notice the dropped matrix
extract(d[[1]], cbptrain[1:2,])

# create a raster stack with the single layer
prd <- stack(d[[1]])
cbp.me <- maxent(x=prd, p =cbptrain)
```

...and error.

### Trying something else.
Found [a tutorial](https://dominicroye.github.io/en/2019/calculating-the-distance-to-the-sea-in-r/ "Calculating Distance to the sea in R") on making a raster file of distance from coastlines. The tutorial uses a "fishnet" of points, and works really well for small countries but not all land masses. Instead of a fishnet of points, maybe I can use the observation points of _C. b-p_


```r
# Load necessary packages
library(rnaturalearth)
library(sf)
library(lwgeom)
```

```
## Warning: package 'lwgeom' was built under R version 3.5.2
```

```
## Linking to liblwgeom 3.0.0beta1 r16016, GEOS 3.7.2, PROJ 5.2.0
```

```r
# Pull 50 km resolution map of all land masses, as a simple feature file
#land50 <- ne_download(scale = 50, type = "land", category = "physical", returnclass = "sf")
world <- ne_countries(scale = 50, returnclass = "sf")

#WGS84 equivalent EPSG 4326

# Transform to UTM so it is easier to calculate distances
#land50 <- st_transform(land50, 4326)
world <- st_transform(world, 4326)

#create the points
points <- occr

#convert lat/long point to Simple Feature
points_sf = st_as_sf(points, coords = c("cbp_nodups$longitude", "cbp_nodups$latitude"), crs = 4326)

#only extract the points in the limits of land masses?
cbpmap <- st_intersection(points_sf, world) 
```

```
## although coordinates are longitude/latitude, st_intersection assumes that they are planar
```

```
## Warning: attribute variables are assumed to be spatially constant
## throughout all geometries
```

```r
#transform world map from polygon shape to line
world <- st_cast(world, "MULTILINESTRING")

#calculation of the distance between the coast and our points
#step takes a bit of time
dist <- st_distance(world, cbpmap)

#distance with unit in meters
head(dist)
```

```
## Units: [m]
## [1]  5075021 15318229  8154173  5746752 12210195 13096919
```
#### Plotting distance from coast

I haven't gotten this chunk to plot properly yet.

```r
library(RColorBrewer)

#create a data.frame with the distance and the coordinates of the points
df <- data.frame(dist = as.vector(dist)/1000,
                    st_coordinates(cbpmap))

# check structure
str(df)

# save it just in case
write.table(df, file = "~/Documents/Josephs_lab/worldwide_salt_analysis/dist_coords_cbp_fromCoasts.csv", sep = ",", eol = "\n", quote = FALSE, row.names = FALSE)

# pick cartographically pleasing color palette
col_dist <- brewer.pal(9, "GnBu")

# plot
ggplot(df, aes(X, Y, fill = dist))+ #variables
         geom_tile()+ #geometry
           scale_fill_gradientn(colours = rev(col_dist))+ #colors for plotting the distance
             labs(fill = "Distance (km)")+ #legend name
             theme_void()+ #map theme
              theme(legend.position = "bottom") #legend position
```



### The Java Thing
I can load the rJava package just fine.

```r
install.packages("rJava")
```

But, if I try to load the library to use it...

```r
library(rJava)
```
I get this Error:

Error: package or namespace load failed for ‘rJava’:
 .onLoad failed in loadNamespace() for 'rJava', details:
  call: dyn.load(file, DLLpath = DLLpath, ...)
  error: unable to load shared object '/Library/Frameworks/R.framework/Versions/3.5/Resources/library/rJava/libs/rJava.so':
  dlopen(/Library/Frameworks/R.framework/Versions/3.5/Resources/library/rJava/libs/rJava.so, 6): Library not loaded: /Library/Java/JavaVirtualMachines/jdk-11.0.1.jdk/Contents/Home/lib/server/libjvm.dylib
  Referenced from: /Library/Frameworks/R.framework/Versions/3.5/Resources/library/rJava/libs/rJava.so
  Reason: no suitable image found.  Did find:
	/Library/Java/JavaVirtualMachines/1.6.0.jdk/Contents/Libraries/libjvm.dylib: mach-o, but wrong architecture
	/Library/Java/JavaVirtualMachines/1.6.0.jdk/Contents/Libraries/libclient.dylib: mach-o, but wrong architecture
In addition: Warning message:
package ‘rJava’ was built under R version 3.5.2 


### In case I can't figure out the Java issue...
The soil dataset does not explicitly list salinity so I am going to use electrical conductivity which [may be okay.](https://www.nrcs.usda.gov/Internet/FSE_DOCUMENTS/nrcs142p2_053280.pdf)

[Appendix A](https://www.earth-syst-sci-data.net/12/299/2020/ "Appendix A") gives the definitions of the values.

On first pass, I want a dataframe that consists of the latitude, longitude, profile_ID (for merging the datasets), and the electrical conductivity measure.

```r
chem_drop <- chem[, c("profile_id", "elco20_value", "elco20_value_avg", "elco20_method", "elco25_value", "elco25_value_avg","elco25_method","elco50_value", "elco50_value_avg", "elco50_method", "elcosp_value", "elcosp_value_avg", "elcosp_method")]

#essentially c-binding the latitude and longitude from "profiles_drop" to the conductivity data, based on profile ID but keeping all of the rows in chem_drop..?
## May or may not be a good idea because I really just want a location and a conductivity...

test_merged <- merge.data.frame(chem_drop, profiles_drop, by="profile_id", all.x = TRUE)
```

