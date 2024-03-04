
setwd('D:/Dropbox/ForestScaling/HARVdata/')
list.files()
#reading in the data
hf.plot <- read.csv('FieldData/HemlockRemoval/hf126-01-plot.csv')
  #  Easting spm. Massachusetts State Plane Meters; Datum = NAD1983 (unit: meter)
  # Longitude wgs1984, Datum = WGS 1984 (unit: degree
hf.tree <- read.csv('FieldData/HemlockRemoval/hf126-02-tree.csv')
hf.tree <- hf.tree[!is.na(hf.tree[,'xcoord']),]
harv <- read.csv('RS/Weinstein_2020_eLife/HARV_2019.csv')
  # UTM Geographic
head(hf.plot)
head(hf.tree)
head(harv)

#lets turn the Lat/Lon columns in hf.plot into UTM
plot.xy <- hf.plot[,c('longitude.wgs1984','latitude.wgs1984')]
plot.xy
#install.packages('sf')
library(sf)
?st_as_sf
#create shapefile with sf package, setting crs
plot.xy.sf <- st_as_sf(data.frame(plot.xy),coords=c(1,2),crs="+init=epsg:4326")
st_crs(plot.xy.sf)
plot.xy.sf
#transform from lon/lat to utm
plot.xy.utm <- st_transform(plot.xy.sf,crs="+proj=utm +zone=18 +datum=WGS84")
plot.xy.utm
#extract utm coordinates
plot.utm <- st_coordinates(plot.xy.utm)


#install.packages('vec2dtransf')
library(vec2dtransf)

U.plots <- unique(hf.plot[,'plot'])
tree.coords.utm <- matrix(nrow=dim(hf.tree)[[1]],ncol=2,
                          dimnames=list(NULL,c('x.utm','y.utm')))
  for(p in 1:length(U.plots)){
    f.rows <- hf.plot[,'plot'] == U.plots[p]
    control.pts <- data.frame(hf.plot[f.rows,c('plot.x','plot.y')],
                              plot.utm[f.rows,])
    transformation <- AffineTransformation(controlPoints=control.pts)
    calculateParameters(transformation)
    f.rows <- hf.tree[,'plot'] == U.plots[p]
    tree.original <- SpatialPoints(data.frame(hf.tree[f.rows,c('xcoord','ycoord')]))
    tree.transformed <- applyTransformation(transformation,tree.original)
    tree.coords.utm[f.rows,]<- coordinates(tree.transformed)
  }  
hf.tree.output <- (data.frame(hf.tree,tree.coords.utm))
write.csv(file = 'hf126_tree_utm.csv', hf.tree.output,row.names=FALSE)









##############
#
# JUNK BELOW
#
############
#If I understand correctly,
#if we subtract the plot.x and plot.y from the utm coords,
#this should give us the origins of the plot coordinate system
origin.x <- plot.utm[,1] - hf.plot[,'plot.x']
origin.y <- plot.utm[,2] - hf.plot[,'plot.y']
plot.origins <- data.frame(plot=hf.plot[,'plot'],origin.x,origin.y)
#I assume the differences of origin location are sampling error
#but I could be wrong if I'm misunderstanding
plot.origins
data.frame(plot.origins,hf.plot)
#so, lets average across each plot
#get unique plot 
U.plots <- unique(plot.origins[,'plot'])
plot.origins.sd <- plot.origins.mean <- matrix(nrow=length(U.plots),ncol=2,dimnames=list(U.plots,c('x.utm','y.utm')))
for(p in 1:length(U.plots)){
  f.rows <- plot.origins[,'plot'] == U.plots[p]
    plot.origins.mean[p,] <- apply(plot.origins[f.rows,c('origin.x','origin.y')],2,mean)
    plot.origins.sd[p,] <- apply(plot.origins[f.rows,c('origin.x','origin.y')],2,sd)
}
plot.origins.mean
#note some slop, like plot 4 & 6 the SD of y-coord of origin is 7m
plot.origins.sd

#so, if we've got the origins right, all we need to do is
#add the x/y coords from hf.tree to the origin
head(hf.tree)
tree.coords.utm <- matrix(nrow=dim(hf.tree)[[1]],ncol=2,dimnames=list(NULL,c('x.utm','y.utm')))
for(p in 1:length(U.plots)){
  f.rows <- hf.tree[,'plot'] == U.plots[p]
  tree.coords.utm[f.rows,'x.utm'] <- hf.tree[f.rows,'xcoord'] +   plot.origins.mean[p,'x.utm']
  tree.coords.utm[f.rows,'y.utm'] <- hf.tree[f.rows,'ycoord'] +   plot.origins.mean[p,'y.utm']
}

