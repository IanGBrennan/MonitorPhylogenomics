library(rase)
library(rgeos)
library(sp)
library(raster)


# here's a function for extracting the node distributions
# and turning it into Spatial Data that we can add to our existing records
process.rase <- function(mcmc.object, distribution, new.directory=NULL,
                         remove.extralimital=NULL, range.shape, point.width=1) {
  dist.frame <- NULL
  res.points <- as.data.frame(mcmc.object[,c(1:(ncol(mcmc.object)-2))])
  for (p in 1:(ncol(res.points)/2)) {
    int <- matrix(nrow=nrow(res.points), ncol=3)
    int <- as.data.frame(int)
    node <- test <- strsplit(colnames(res.points)[p], "_")
    int[,1] <- node[[1]][1]
    int[,2] <- res.points[,p]
    int[,3] <- res.points[,(p+(ncol(res.points)/2))]
    dist.frame <- rbind(dist.frame, int)
  }
  colnames(dist.frame) <- c("Name_in_Tree", "Latitude", "Longitude")
  
  if (remove.extralimital==TRUE){
    spatial.frame <- SpatialPointsDataFrame(coords = dist.frame[, c("Longitude", "Latitude")], 
                                            data = dist.frame)
    input.shp <- shapefile(range.shape)
    projection(input.shp) <- projection(spatial.frame) # this line of code written by Carlos Pavón
    inside.data <- spatial.frame[!is.na(over(spatial.frame, as(input.shp, "SpatialPolygons"))), ]
    dist.frame <- as.data.frame(inside.data[,1:3])
  }
  
  # create an ouput folder to catch the plots
  if(!is.null(new.directory)) {
    curr <- getwd()
    dir.path <- paste0(curr,"/", new.directory)
    call <- paste("mkdir", dir.path)
    system(call)
  }
  
  all.taxa <- unique(dist.frame$Name_in_Tree)
  all.sp <- list()
  all.hull <- list()
  all.owin <- list()
  
  if(!is.null(new.directory)) {
    sbbox <- make_bbox(lon = distribution$Longitude, lat = distribution$Latitude, f = .1)
    sq_map <- get_map(location = sbbox, maptype = "terrain", source = "google")
  }
  
  for (p in 1:length(unique(dist.frame$Name_in_Tree))) {
    cat("iteration", p, "of", length(unique(dist.frame$Name_in_Tree)), "\n") #keep track of what tree/loop# we're on
    current.taxon <- all.taxa[[p]]
    current.data <- subset(dist.frame, dist.frame$Name_in_Tree == current.taxon)
    if (nrow(current.data) > 1000) {
      current.data <- current.data[sample(nrow(current.data), 1000, replace=F), ]
    } 
    points <- current.data[,c(2,3)]
    all.sp[[p]] <- distribution.sp <- SpatialPoints(points)
    all.hull[[p]] <- distribution.hull <- gBuffer(distribution.sp, width=point.width)
    all.owin[[p]] <- as.owin(all.hull[[p]])
    
    # only plot the maps below if you need to check/clean the data
    if(!is.null(new.directory)) {
      fortified.data <- fortify(distribution.hull)
      pdf(paste(dir.path,"/", current.taxon, ".pdf", sep=""))
      print(ggmap(sq_map) 
            + geom_polygon(data = fortified.data, aes(x = lat, y = long, group=group, fill="tomato"), alpha=0.3)
            + geom_point(data = current.data, mapping = aes(x=Longitude, y=Latitude), color="tomato3"))
      dev.off()
    }
  }
  names(all.sp) <- unique(dist.frame$Name_in_Tree)
  names(all.hull) <- unique(dist.frame$Name_in_Tree)
  
  if(!is.null(new.directory)) {
    plot.location.call <- paste("plots printed to directory:", dir.path)
    print(plot.location.call)
  }
  
  return(list(DistData=dist.frame, SpatialPoints=all.sp, ConvexHulls=all.hull, OWin=all.owin))
}
# mcmc.object: the output from 'rase' (after burnin via 'mcmc')
# distribution: the lat/long data frame used as an input
# new directory: if you want plots of the nodes, give a name for the new directory
# remove.extralimital: would you like to drop points outside of a boundary?
# range.shape: if 'yes' above, provide the path to a .shx shape file (e.g."/Users/Ian/Desktop/Australia.shx")
### the one thing I'm thinking about now is that this function gets the points for each node, and then
### when there's a cladogenetic event, we want to compare a node's distribution to a branch, but we 
### end up comparing a node to another node that may not be at the same time point. 
### we could use the 'rase.slice' function from rase, but it could be a massive headache to try and
### and figure out the ranges on branches each time there's a speciation event. Use tree.slice maybe? 
### On second thought,
### this would be WAAAAY to much work, considering how many branches there would be towards the tips
### I'll leave this here for when I inevitably think of this again later.