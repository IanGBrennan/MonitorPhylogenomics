require(sp)
require(adehabitatHR)
require(rgeos)
require(rworldmap); require(ggmap)
require(phangorn)
require(R.utils)
require(dplyr)
require(phytools)

## Right now the function determines sympatry using gOverlaps which is binary
## It would be nice to incorporate the amount of overlap, which we could then try
## to correlate with the strength of character displacement. E.g. plot pairwise trait
## distance (felsen) against S*overlap
## We could try to do this using gIntersection and gArea. 
## Here's an example of how it can be done:
          #   v.ac <- tips$ConvexHulls$Varanus.primordius
          #   v.ba <- tips$ConvexHulls$Varanus.gouldii  
          #   gOverlaps(v.ac, v.ba)
          #   v.acba <- gIntersection(v.ac, v.ba)
          #   v.ac.overlap <- gArea(v.acba)/gArea(v.ac)
          #   v.ba.overlap <- gArea(v.acba)/gArea(v.ba)

#data(Anolis.data)
##Create a geography.object with a modified edge matrix
##First, specify which region each branch belonged to:
#Anolis.regions<-c(rep("cuba",14),rep("hispaniola",17),"puerto_rico")
#map<-cbind(Anolis.data$phylo$edge,Anolis.regions)
#  # or  #
#map <- read.csv("/Users/Ian/Google.Drive/ANU Herp Work/Adaptive Radiation/Distribution_Data/ANOLIS_TEST_DATA.csv", header=T)
#    map <- map[,c("Name_in_Tree", "Latitude", "Longitude")]
##CreateGeoObject(Anolis.data$phylo,map=Anolis.map)
#phylo <- Anolis.data$phylo
#
#distribution <- read.csv("/Users/Ian/Google.Drive/ANU Herp Work/Adaptive Radiation/Distribution_Data/ANOLIS_TEST_DATA.csv", header=T)
#    distribution <- distribution[,c("Name_in_Tree", "Latitude", "Longitude")]
#        distribution <- distribution[complete.cases(distribution),] 
#
CreateGeoObject_SP <- function (phylo, map, point.width=0.5) {
  if (any(grepl("___", phylo$tip.label))) {
    stop("script will not work with '___' in tip labels; remove extra underscores")
  }
  beginning <- Sys.time()
  
  TWOLETTERS <- paste(rep(LETTERS, each = 26), LETTERS, sep = "")
  THREELETTERS <- paste(rep(TWOLETTERS, each = 26), LETTERS, 
                        sep = "")
  nodeDist <- vector(mode = "numeric", length = phylo$Nnode) # vector of the internal nodes
  totlen <- length(phylo$tip.label) # the number of tips
  root <- totlen + 1 # number of tips +1 for the root
  heights <- nodeHeights(phylo)
  for (i in 1:dim(phylo$edge)[1]) {
    nodeDist[[phylo$edge[i, 1] - totlen]] <- heights[i] # heights above root of each internal node
  }
  nodeDist <- c(nodeDist, max(heights)) # add the total height on the end of the list of nodeheights
  nodeDiff <- diff(nodeDist) # get the distances between the internal nodes
  flag = 0
  if (sum(nodeDiff < 0) > 0) {
    node.order <- match(rank(heights[, 1], ties.method = "min"), 
                        seq(1, by = 2, len = phylo$Nnode))
    node.order <- node.order + totlen
    old.edge <- phylo$edge
    old.phylo <- phylo
    phylo$edge[, 1] <- node.order
    for (j in 1:length(phylo$edge[, 2])) {
      if (phylo$edge[j, 2] > totlen) {
        phylo$edge[j, 2] <- phylo$edge[, 1][match(phylo$edge[j, 
                                                             2], old.edge[, 1])]
      }
    }
    nodeDist <- vector()
    for (i in 1:dim(phylo$edge)[1]) {
      nodeDist[[phylo$edge[i, 1] - totlen]] <- heights[i]
    }
    nodeDist <- c(nodeDist, max(heights))
    nodeDiff <- diff(nodeDist)
    flag = 1
  }
  mat <- matrix(nrow = 0, ncol = 3)
  counter_three_letters <- 0
  
  # I can copy the names for MRCAs here and apply them to nodes on the tree 
  for (i in 1:phylo$Nnode) {
    other <- phylo$edge[phylo$edge[, 1] == i + totlen, 2]
    for (b in other) {
      int <- matrix(ncol = 3)
      int[1] <- i + totlen
      if (b > totlen) {
        counter_three_letters <- counter_three_letters + 
          1
        int[2] <- paste(".", THREELETTERS[counter_three_letters], 
                        sep = "")
        int[3] <- b
      }
      else {
        int[2] <- phylo$tip.label[[b]]
        int[3] <- 0
      }
      mat <- rbind(mat, int)
    }
  }
  #######
  ####### if (any(class(map) == "phylo")) {
  #  for (i in 1:length(mat[, 1])) {
  #    if (mat[i, 3] == 0) {
  #      mat[i, 3] <- as.character(match(mat[i, 2], phylo$tip.label))
  #    }
  #  }
  #  maps.object <- map$maps
  #  brchange <- vapply(maps.object, function(x) any(length(names(x)) > 
  #                                                    1), 1)
  #  if (sum(brchange) != 0) {
  #    newtimes <- vector()
  #    for (i in 1:length(brchange)) {
  #      if (brchange[i] != 0) {
  #        intlen <- length(maps.object[[i]])
  #        inttime <- vector()
  #        for (j in 1:(intlen - 1)) {
  #          inttime <- c(inttime, maps.object[[i]][j])
  #          nt <- as.numeric(heights[i, 1] + sum(inttime))
  #          newtimes <- c(newtimes, nt)
  #        }
  #      }
  #    }
  #    old.Dist <- nodeDist
  #    old.Diff <- nodeDiff
  #    nodeDist <- sort(c(nodeDist, newtimes))
  #    nodeDiff <- diff(nodeDist)
  #  }
  #  nat <- list()
  #  nodecount = 1
  #  for (i in 1:length(nodeDiff)) {
  #    if (i == 1) {
  #      hold.m <- mat[as.numeric(mat[, 1]) <= (totlen + 
  #                                               nodecount), c(2, 3)]
  #      int <- dim(hold.m)[1]
  #      for (m in 1:int) {
  #        hold.m[m, 2] <- names(maps.object[[which(phylo$edge[, 
  #                                                            2] == as.numeric(hold.m[m, 2]))]])[1]
  #      }
  #      nat[[i]] <- hold.m
  #    }
  #    else {
  #      if (nodeDist[i] %in% old.Dist) {
  #        P <- mat[as.numeric(mat[, 1]) <= (totlen + 
  #                                            nodecount), c(2, 3)]
  #        hold.m <- rbind(P[as.numeric(P[, 2]) <= totlen, 
  #                          ], P[as.numeric(P[, 2]) > (totlen + nodecount), 
  #                               ])
  #        int <- dim(hold.m)[1]
  #        for (m in 1:int) {
  #          iden <- which(phylo$edge[, 2] == as.numeric(hold.m[m, 
  #                                                             2]))
  #          if (brchange[iden] == 0 || !(hold.m[m, 1] %in% 
  #                                       nat[[i - 1]][, 1])) {
  #            hold.m[m, 2] <- names(maps.object[[iden]])[1]
  #          }
  #          else {
  #            num = 1
  #            while (round(nodeDist[i + 1] - old.Dist[phylo$edge[iden, 
  #                                                               1] - totlen], 6) > round(sum(maps.object[[iden]][1:num]), 
  #                                                                                        6)) {
  #              num = num + 1
  #            }
  #            hold.m[m, 2] <- names(maps.object[[iden]])[num]
  #          }
  #        }
  #        nat[[i]] <- hold.m
  #      }
  #      else {
  #        P <- mat[as.numeric(mat[, 1]) <= (totlen + 
  #                                            nodecount), c(2, 3)]
  #        hold.m <- rbind(P[as.numeric(P[, 2]) <= totlen, 
  #                          ], P[as.numeric(P[, 2]) > (totlen + nodecount), 
  #                               ])
  #        int <- dim(hold.m)[1]
  #        for (m in 1:int) {
  #          iden <- which(phylo$edge[, 2] == as.numeric(hold.m[m, 
  #                                                             2]))
  #          if (brchange[iden] == 0 | !(hold.m[m, 1] %in% 
  #                                      nat[[i - 1]][, 1])) {
  #            hold.m[m, 2] <- names(maps.object[[iden]])[1]
  #          }
  #          else {
  #            num = 1
  #            while (round(nodeDist[i + 1] - old.Dist[phylo$edge[iden, 
  #                                                               1] - totlen], 6) > round(sum(maps.object[[iden]][1:num]), 
  #                                                                                        6)) {
  #              num = num + 1
  #            }
  #            hold.m[m, 2] <- names(maps.object[[iden]])[num]
  #          }
  #        }
  #        nat[[i]] <- hold.m
  #      }
  #    }
  #    if (nodeDist[i + 1] %in% old.Dist) {
  #      nodecount = nodecount + 1
  #    }
  #  }
  #  geography.matrix <- list()
  #  for (i in 1:length(nodeDiff)) {
  #    timei <- unlist(nat[[i]])
  #    var.list <- timei[, 1]
  #    len = length(var.list)
  #    int.mat <- matrix(nrow = len, ncol = len)
  #    rownames(int.mat) <- var.list
  #    colnames(int.mat) <- var.list
  #    diag(int.mat) <- 1
  #    for (j in 1:length(var.list)) {
  #      for (k in 1:length(var.list)) {
  #        if ((lower.tri(int.mat, diag = TRUE)[j, k] == 
  #             TRUE)) {
  #          int.mat[j, k] <- ifelse(timei[match(var.list[j], 
  #                                              timei[, 1]), 2] == timei[match(var.list[k], 
  #                                                                             timei[, 1]), 2], 1, 0)
  #        }
  #      }
  #    }
  #    int.mat[upper.tri(int.mat) == TRUE] <- t(int.mat)[upper.tri(t(int.mat)) == 
  #                                                        TRUE]
  #    geography.matrix[[i]] <- int.mat
  #  }
  #}
  #######
  if (is.data.frame(map)) {
    nat <- list()
    for (i in 1:length(nodeDiff)) {
      if (i == 1) {
        nat[[i]] <- list(mat[mat[, 1] == (totlen + i), 
                             2])
      }
      else {
        IN <- vector()
        P <- mat[as.numeric(mat[, 1]) <= (totlen + i), 
                 c(2, 3)]
        IN <- c(IN, P[P[, 2] == "0", 1], P[as.numeric(P[, 
                                                        2]) > (totlen + i), 1])
        nat[[i]] <- list(IN)
      }
    }
    
    # I think this is the gap to fit my code, it would make a matrix for each list in the geo_object list
  #####
   # for (i in 1:length(mat[, 1])) {
   #   if (mat[i, 3] == 0) {
   #     mat[i, 3] <- as.character(match(mat[i, 2], phylo$tip.label))
   #   }
   # }
   # for (i in 1:length(mat[, 1])) {
   #   if (mat[i, 3] == 0) {
   #     mat[i, 3] <- as.character(match(mat[i, 2], phylo$tip.label))
   #   }
   # }
   # if (flag == 0) {
   #   map.v <- map[match(mat[, 3], map[, 2]), 3]
   # }
   # ## I think I can insert the overlaps here, or try to call them in the third column
   # if (flag == 1) {
   #   map2 <- cbind(phylo$edge, map[, 3])
   #   map.v <- map2[match(mat[, 3], map2[, 2]), 3]
   # }
   # mat <- data.frame(mat, map.v)
   # geography.matrix <- list()
   # for (i in 1:phylo$Nnode) {
   #   var.list <- unlist(nat[[i]])
   #   len = length(var.list)
   #   int.mat <- matrix(nrow = len, ncol = len)
   #   rownames(int.mat) <- var.list
   #   colnames(int.mat) <- var.list
   #   diag(int.mat) <- 1
   #   for (j in 1:length(var.list)) {
   #     for (k in 1:length(var.list)) {
   #       if ((lower.tri(int.mat, diag = TRUE)[j, k] == 
   #            TRUE)) {
   #         int.mat[j, k] <- ifelse(mat[match(var.list[j], 
   #                                           mat[, 2]), 4] == mat[match(var.list[k], 
   #                                                                      mat[, 2]), 4], 1, 0)
   #       }
   #     }
   #   }
   #   int.mat[upper.tri(int.mat) == TRUE] <- t(int.mat)[upper.tri(t(int.mat)) == 
   #                                                       TRUE]
   #   geography.matrix[[1]] <- int.mat
   # }
  ##### 
    int.names <- rownames(vcvPhylo(phylo)) # this does not include a crown taxon representing all tips
    counter_three_letters <- 1
    for (p in (Ntip(phylo)+1):length(int.names)) {
      int.names[p] <- paste(".", THREELETTERS[counter_three_letters], 
                            sep = "")
      counter_three_letters <- counter_three_letters+1
    } # this does not g
    
    all.taxa <- unique(map[,1])
    all.sp <- list()
    all.hull <- list()

    for (p in 1:length(unique(map[,1]))) {
      #cat("iteration", p, "of", length(unique(map$Name_in_Tree)), "\n") #keep track of what tree/loop# we're on
      current.taxon <- all.taxa[[p]]
      current.data <- subset(map, map[,1] == current.taxon)
      if (nrow(current.data) > 1000) {
        current.data <- current.data[sample(nrow(current.data), 1000, replace=F), ]
      } 
      points <- current.data[,c(2,3)]
      all.sp[[p]] <- distribution.sp <- SpatialPoints(points)
      all.hull[[p]] <- distribution.hull <- gBuffer(distribution.sp, width=point.width)
    }
    names(all.sp) <- unique(map[,1])
    names(all.hull) <- unique(map[,1])
    
    # make a loop that creates distributions at nodes by combining daughter distributions
    for (kk in (length(int.names)+1):(length(phylo$tip.label)+2)){
      descendantz <- Descendants(phylo, kk, type="children")
      combined.range <- gUnion(all.hull[[descendantz[1]]], all.hull[[descendantz[2]]])
      all.hull[[kk]] <- combined.range
    }
    all.hull[(Ntip(phylo)+1):(length(int.names))] <- all.hull[(Ntip(phylo)+2):(length(int.names)+1)]
    all.hull <- all.hull[1:length(int.names)]
    names(all.hull)[(Ntip(phylo)+1) : (length(all.hull))] <- int.names[(Ntip(phylo)+1):(length(int.names))]
    
    unique.combo <-  combn(int.names, m = 2)
    range_overlap <- rep(NA, ncol(unique.combo)) # make a ma
    arrange <- data.frame(species1 = unique.combo[1,], 
                          species2 = unique.combo[2,], 
                          range_overlap, stringsAsFactors=F)
    for (t in 1:length(arrange[,1])) {
      overlap <- gOverlaps(all.hull[[arrange[t,1]]], all.hull[[arrange[t,2]]])
      arrange[t,3] <- overlap
    }
    arrange[,3][arrange[,3] == "TRUE"] <- 1
    
    master.matrix <- matrix(nrow=length(int.names), ncol=length(int.names))
    master.matrix <- as.data.frame(master.matrix)
    rownames(master.matrix) <- int.names
    colnames(master.matrix) <- int.names
    #diag(master.matrix) <- 999
    for (j in 1:length(int.names)) {
      for (k in 1:length(int.names)) {
        if (j == k) {
          master.matrix[j,k] <- 1
        } else {
          master.matrix[j,k] <- filter(arrange, (species1==int.names[j] & species2==int.names[k]) |
                                   (species2==int.names[j]) & species1==int.names[k])$range_overlap
        }
      }
    }
    
    geography.matrix <- list()
    for (i in 1:phylo$Nnode) {
      var.list <- unlist(nat[[i]])
      int.mat <- master.matrix[which(rownames(master.matrix)%in%var.list),]
      int.mat <- int.mat[,which(colnames(int.mat)%in%var.list)]
      int.mat <- as.matrix(int.mat) 
      geography.matrix[[i]] <- int.mat
    }
    
   # geography.matrix <- list()
   # for (i in 1:phylo$Nnode) {
   #   var.list <- unlist(nat[[i]])
   #   len = length(var.list)
   #   int.mat <- matrix(nrow=len, ncol=len)
   #   rownames(int.mat) <- var.list
   #   colnames(int.mat) <- var.list
   #   diag(int.mat) <- 1
   #   for (j in 1:length(var.list)) {
   #     for (k in 1:length(var.list)) {
   #       if (j == k) {
   #         int.mat[j,k] <- 1
   #       } else {
   #         int.mat[j,k] <- filter(arrange, (species1==var.list[j] & species2==var.list[k]) |
   #                          (species2==var.list[j]) & species1==var.list[k])$range_overlap
   #         #int.mat[j,k] <- filter(all.matches, var.list[1]==Species1)$overlap
   #       }
   #     }
   #   }
   #   geography.matrix[[i]] <- int.mat
   # }
  }
  end <- Sys.time()
  duration <- format(end-beginning)
  print(paste("Computation time :", duration))
  
  return(list(geography.object = geography.matrix, times = nodeDist, 
              spans = nodeDiff, name.matrix = mat))
}
#test <- CreateGeoObject_SP(Anolis.data$phylo, distribution)





# I need to create a GeoObject that accounts for pairwise overlap between taxa of two trees
## this could be more complicated than I think, but should be able to do it using the basic 
### format of the single tree approach below
## things to consider: ancestral nodes in the two trees need to be named differently.

### I've done this already, but now I need to be able to take in a RASE output, which
## already includes distributions for ancestral nodes. 
### RASE output
CreateCoEvoGeoObject_SP <- function (phy1, phy2, map, rase.obj1, rase.obj2, point.width=0.25) {
  print("Creating pairwise comparison of distributional overlap between all taxa in trees 1 and 2")
  beginning <- Sys.time()
  
  if (any(grepl("___", phy1$tip.label))) {
    stop("script will not work with '___' in tip labels; remove extra underscores")
  }
  TWOLETTERS <- paste(rep(LETTERS, each = 26), LETTERS, sep = "") # create two letter codes
  THREELETTERS <- paste(rep(TWOLETTERS, each = 26), LETTERS, sep = "") # create three letter codes
  
  tree.height.diff <- abs(max(nodeHeights(phy1)) - max(nodeHeights(phy2))) # get the difference in tree heights
  
# sort the phylogenies to make phylo1 the older of the two, phylo2 the younger, unless they're the same age
  if (max(nodeHeights(phy1))>max(nodeHeights(phy2))) {
    phylo1 <- phy1; phylo2 <- phy2
    rase1 <- rase.obj1; rase2 <- rase.obj2
    print(paste("Note: input tree1 is", tree.height.diff, "units older than input tree2", sep=" "))
  } else if (max(nodeHeights(phy2))>max(nodeHeights(phy1))) {
    phylo1 <- phy2; phylo2 <- phy1
    rase1 <- rase.obj2; rase2 <- rase.obj1
    print(paste("Note: input tree2 is", tree.height.diff, "units older than input tree1", sep=" "))
  } else {
    phylo1 <- phy1; phylo2 <- phy2
    rase1 <- rase.obj1; rase2 <- rase.obj2
    print(paste("Note: trees are of equal depth,", tree.height.diff, "units", sep=" "))
  }

# pull together information about each tree's branching processes
  nodeDist1 <- vector(mode = "numeric", length = phylo1$Nnode) # vector of the internal nodes
    nodeDist2 <- vector(mode = "numeric", length = phylo2$Nnode)
  totlen1 <- length(phylo1$tip.label) # the number of tips
    totlen2 <- length(phylo2$tip.label)
      totlen <- max(c(totlen1, totlen2))
      max.depth <- max(c(max(nodeHeights(phylo1)), max(nodeHeights(phylo2)))) # max height of the older tree
  root1 <- totlen1 + 1 # number of tips +1 for the root
    root2 <- totlen2 + 1
  heights1 <- nodeHeights(phylo1)
    heights2 <- nodeHeights(phylo2)

# get the occurrences of all branching times (relative to the max height of both trees - max.depth)
# this depends on if trees are of equal height or not!
  if (max(nodeHeights(phylo1))>max(nodeHeights(phylo2))) {
    nh1 <- round(heights1[,2], digits=6)
    nh2 <- round(max.depth - c(max(nodeHeights(phylo2)),(max(nodeHeights(phylo2)) - heights2[,2])), digits=6)
    allnodes <- unique(sort(c(0, nh1, nh2)))
  } else if (max(nodeHeights(phylo1))==max(nodeHeights(phylo2))){
    nh1 <- heights1[,2]
    nh2 <- heights2[,2]
    allnodes <- unique(sort(c(0, nh1, nh2)))
  } else {
    stop("your trees are all jacked up, stop and have a look to see what's wrong")
  }
  
  for (i in 1:dim(phylo1$edge)[1]) {
    nodeDist1[[phylo1$edge[i, 1] - totlen1]] <- heights1[i] # heights above root of each internal node
  }
    for (i in 1:dim(phylo2$edge)[1]) {
      nodeDist2[[phylo2$edge[i, 1] - totlen2]] <- heights2[i] # heights above root of each internal node
    }
  nodeDist1 <- c(nodeDist1, max(heights1)) # add the total height on the end of the list of nodeheights (end of the tree)
    nodeDist2 <- c(nodeDist2, max(heights2))
      nodeDist <- unique(sort(c(nodeDist1, nodeDist2)))
  nodeDiff1 <- diff(nodeDist1) # get the distances between the internal nodes
    nodeDiff2 <- diff(nodeDist2)
      nodeDiff <- diff(nodeDist)
  flag = 0
# redo the branching order if the tree branches are mixed up
  if (sum(nodeDiff1 < 0) > 0) {
    node.order <- match(rank(heights1[, 1], ties.method = "min"), 
                        seq(1, by = 2, len = phylo1$Nnode))
    node.order <- node.order + totlen1
    old.edge <- phylo1$edge
    old.phylo <- phylo1
    phylo1$edge[, 1] <- node.order
    for (j in 1:length(phylo1$edge[, 2])) {
      if (phylo1$edge[j, 2] > totlen1) {
        phylo1$edge[j, 2] <- phylo1$edge[, 1][match(phylo1$edge[j, 
                                                             2], old.edge[, 1])]
      }
    }
    nodeDist1 <- vector()
    for (i in 1:dim(phylo1$edge)[1]) {
      nodeDist1[[phylo1$edge[i, 1] - totlen1]] <- heights1[i]
    }
    nodeDist1 <- c(nodeDist1, max(heights1))
    flag = 1
  }
  if (sum(nodeDiff2 < 0) > 0) {
    node.order <- match(rank(heights2[, 1], ties.method = "min"), 
                        seq(1, by = 2, len = phylo2$Nnode))
    node.order <- node.order + totlen2
    old.edge <- phylo2$edge
    old.phylo <- phylo2
    phylo2$edge[, 1] <- node.order
    for (j in 1:length(phylo2$edge[, 2])) {
      if (phylo2$edge[j, 2] > totlen2) {
        phylo2$edge[j, 2] <- phylo2$edge[, 1][match(phylo2$edge[j, 
                                                                2], old.edge[, 1])]
      }
    }
    nodeDist2 <- vector()
    for (i in 1:dim(phylo2$edge)[1]) {
      nodeDist2[[phylo2$edge[i, 1] - totlen2]] <- heights2[i]
    }
    nodeDist2 <- c(nodeDist2, max(heights2))
    flag = 1
  }
  nodeDist <- unique(sort(c(nodeDist1, nodeDist2))) # I don't think I use this for anything later
  
# create a matrix for each tree which includes the edge [,4] starting and [,5] ending times,
# [,1] parent and [,6] daughter node numbers, [,2] names of nodes/tips, and [,3] adjusted daughter node index (0=tip)
  mat1 <- matrix(nrow = 0, ncol = 7)
  counter_three_letters <- 0
  for (i in 1:phylo1$Nnode) {
    other <- phylo1$edge[phylo1$edge[, 1] == i + totlen1, 2]
    for (b in other) {
      int <- matrix(ncol = 6)
      int[1] <- i + totlen1
      if (b > totlen1) {
        counter_three_letters <- counter_three_letters + 
          1
        int[2] <- paste(".", THREELETTERS[counter_three_letters], 
                        sep = "")
        int[3] <- b
      }
      else {
        int[2] <- phylo1$tip.label[[b]]
        int[3] <- 0
      }
      int[4] <- nodeheight(phylo1, (phylo1$edge[phylo1$edge[, 1] == i + totlen1, 1][1]))
      int[5] <- nodeheight(phylo1, b)
      int[6] <- b
      int[7] <- "tree1"
      mat1 <- rbind(mat1, int)
    }
  }
  mat1 <- mat1[match(phylo1$edge[,2],mat1[,6]),] # sort the matrix so it matches the branching order
      rownames(mat1) <- NULL; mat1 <- as.data.frame(mat1)
####################################################################################################
# NOW THERE ARE TWO OPTIONS:
    # if trees are of equal height, it's easy, and we follow the first section below
    # if trees are of unequal height, it's more complicated, and we follow the second section
  if (max(nodeHeights(phylo2))==max.depth){
    mat2 <- matrix(nrow = 0, ncol = 7)
    # I can copy the names for MRCAs here and apply them to nodes on the tree 
    for (i in 1:phylo2$Nnode) {
      other <- phylo2$edge[phylo2$edge[, 1] == i + totlen2, 2]
      for (b in other) {
        int <- matrix(ncol = 6)
        int[1] <- i + totlen2
        if (b > totlen2) {
          counter_three_letters <- counter_three_letters + 
            1
          int[2] <- paste(".", THREELETTERS[counter_three_letters], 
                          sep = "")
          int[3] <- b
        }
        else {
          int[2] <- phylo2$tip.label[[b]]
          int[3] <- 0
        }
        int[4] <- tree.height.diff + nodeheight(phylo2, (phylo2$edge[phylo2$edge[, 1] == i + totlen2, 1][1]))
        int[5] <- tree.height.diff + nodeheight(phylo2, b)
        int[6] <- b
        int[7] <- "tree2"
        mat2 <- rbind(mat2, int)
      }
    }
    mat2 <- mat2[match(phylo2$edge[,2],mat2[,6]),] # sort the matrix so it matches the branching order
    
    # [,1] is the parent node number
    # [,2] are the daughter nodes or tips by name
    # [,3] are the daughter nodes (#>0) or tips (0)
    # [,4] are the starting ages of the branches (above the root)
    # [,5] are the maximum ages of the branches (above the root)
    # [,6] are the original daughter nodes or tips
    mat <- as.data.frame(rbind(mat1, mat2))
        mat[,4] <- as.numeric(as.character(mat[,4])) # treat values as numbers
            mat[,5] <- as.numeric(as.character(mat[,5])) # treat values as numbers

# this is a HUGE step, make a list of vectors of taxa existing in a given period
# it needs to be in the same order as the branching events, so that our geo.object
# matches MatrixA later on
    if (is.data.frame(map)) {
      nat <- list()
      for (i in 1:(length(allnodes)-1)) {
        if (i == 1) {
          nat[[i]] <- list(mat[mat[,4]==0,2])
        } else {
          starts <- filter(mat, mat[,4]>=allnodes[1] & mat[,4]<=allnodes[i])
          ends <- starts[starts[,5]>allnodes[i],2]
          nat[[i]] <- list(ends)
        }
      }
    }

# collect the names of tips and name the internal nodes (this doesn't name the crown!)
    int.names1 <- rownames(vcvPhylo(phylo1)) # this does not include a crown taxon representing all tips
    counter_three_letters <- 1
    for (p in (Ntip(phylo1)+1):length(int.names1)) {
      int.names1[p] <- paste(".", THREELETTERS[counter_three_letters], 
                             sep = "")
      counter_three_letters <- counter_three_letters+1
    } 
    int.names <- int.names1
    # CREATE ANOTHER LOOP FOR THE SECOND TREE
    int.names2 <- rownames(vcvPhylo(phylo2))
    for (p in (Ntip(phylo2)+1):length(int.names2)) {
      int.names2[p] <- paste(".", THREELETTERS[counter_three_letters], 
                             sep = "")
      counter_three_letters <- counter_three_letters+1
    } 
    # THEN APPEND THE NAMES TO THE int.names OBJECT
    int.names <- append(int.names, int.names2)
    
    all.sp1 <- list(); all.sp2 <- list()
    all.hull1 <- list(); all.hull2 <- list()
    
# make a loop that creates distributions of all tips FOR THE FIRST TREE
    for (p in 1:Ntip(phylo1)) {
      #cat("iteration", p, "of", length(unique(map$Name_in_Tree)), "\n") #keep track of what tree/loop# we're on
      current.taxon <- phylo1$tip.label[[p]]
      current.data <- subset(map, map[,1] == current.taxon)
      if (nrow(current.data) > 1000) {
        current.data <- current.data[sample(nrow(current.data), 1000, replace=F), ]
      } 
      points <- current.data[,c(2,3)]
      all.sp1[[p]] <- distribution.sp <- SpatialPoints(points)
      all.hull1[[p]] <- distribution.hull <- gBuffer(distribution.sp, width=point.width)
    }
    names(all.sp1) <- phylo1$tip.label
    names(all.hull1) <- phylo1$tip.label
    
    # then a loop that creates distributions at nodes by combining daughter distributions
    for (kk in (length(int.names1)+1):(length(phylo1$tip.label)+2)){
      descendantz <- Descendants(phylo1, kk, type="children")
      combined.range <- gUnion(all.hull1[[descendantz[1]]], all.hull1[[descendantz[2]]])
      all.hull1[[kk]] <- combined.range
    }
    all.hull1[(Ntip(phylo1)+1):(length(int.names1))] <- all.hull1[(Ntip(phylo1)+2):(length(int.names1)+1)]
    all.hull1 <- all.hull1[1:length(int.names1)]
    names(all.hull1)[(Ntip(phylo1)+1) : (length(all.hull1))] <- int.names1[(Ntip(phylo1)+1):(length(int.names1))]
    
# make a loop that creates distributions of all tips FOR THE SECOND TREE
    for (p in 1:Ntip(phylo2)) {
      #cat("iteration", p, "of", length(unique(map$Name_in_Tree)), "\n") #keep track of what tree/loop# we're on
      current.taxon <- phylo2$tip.label[[p]]
      current.data <- subset(map, map[,1] == current.taxon)
      if (nrow(current.data) > 1000) {
        current.data <- current.data[sample(nrow(current.data), 1000, replace=F), ]
      } 
      points <- current.data[,c(2,3)]
      all.sp2[[p]] <- distribution.sp <- SpatialPoints(points)
      all.hull2[[p]] <- distribution.hull <- gBuffer(distribution.sp, width=point.width) # you can adjust the buffer of each SP with 'point.width'
    }
    names(all.sp2) <- phylo2$tip.label
    names(all.hull2) <- phylo2$tip.label
    
    # then a loop that creates distributions at nodes by combining daughter distributions
    for (kk in (length(int.names2)+1):(length(phylo2$tip.label)+2)){
      descendantz <- Descendants(phylo2, kk, type="children")
      combined.range <- gUnion(all.hull2[[descendantz[1]]], all.hull2[[descendantz[2]]])
      all.hull2[[kk]] <- combined.range
    }
    all.hull2[(Ntip(phylo2)+1):(length(int.names2))] <- all.hull2[(Ntip(phylo2)+2):(length(int.names2)+1)]
    all.hull2 <- all.hull2[1:length(int.names2)]
    names(all.hull2)[(Ntip(phylo2)+1) : (length(all.hull2))] <- int.names2[(Ntip(phylo2)+1):(length(int.names2))]
    
    all.hull <- append(all.hull1, all.hull2) # combine the two hull lists together
    
    unique.combo <-  combn(int.names, m = 2) # unique combination of all taxon names for pairwise comparison
    range_overlap <- rep(NA, ncol(unique.combo)) 
    arrange <- data.frame(species1 = unique.combo[1,], 
                          species2 = unique.combo[2,], 
                          range_overlap, stringsAsFactors=F)
    for (t in 1:length(arrange[,1])) {
      overlap <- gOverlaps(all.hull[[arrange[t,1]]], all.hull[[arrange[t,2]]])
      arrange[t,3] <- overlap
    }
    arrange[,3][arrange[,3] == "TRUE"] <- 1
    
    master.matrix <- matrix(nrow=length(int.names), ncol=length(int.names))
    master.matrix <- as.data.frame(master.matrix)
    
    ordered.names <- mat[,2] # name the rows/columns of the matrix by the properly ordered branching times
    rownames(master.matrix) <- ordered.names
    colnames(master.matrix) <- ordered.names
    
# loop through and populate the matrix according to our pairwise comparisons
    for (j in 1:length(ordered.names)) {
      for (k in 1:length(ordered.names)) {
        if (j == k) {
          master.matrix[j,k] <- 1
        } else {
          master.matrix[j,k] <- filter(arrange, (species1==ordered.names[j] & species2==ordered.names[k]) |
                                         (species2==ordered.names[j]) & species1==ordered.names[k])$range_overlap
        }
      }
    }
  } 
###############################################################################################################
# this is the alternative, if the trees are of different heights, it makes a dummy name/range for the root node 
  # of the shorter/younger tree, so that it gets processed properly by the 'CreateModelCoevolution' function, 
  # all else is the same so you can refer to the above for notes if necessary
  else {
    mat2 <- matrix(nrow = 0, ncol = 7)
    last.count <- counter_three_letters
    counter_three_letters <- counter_three_letters + 1
    for (i in 1:phylo2$Nnode) {
      other <- phylo2$edge[phylo2$edge[, 1] == i + totlen2, 2]
      for (b in other) {
        int <- matrix(ncol = 6)
        int[1] <- i + totlen2
        if (b > totlen2) {
          counter_three_letters <- counter_three_letters + 
            1
          int[2] <- paste(".", THREELETTERS[counter_three_letters], 
                          sep = "")
          int[3] <- b
        }
        else {
          int[2] <- phylo2$tip.label[[b]]
          int[3] <- 0
        }
        int[4] <- tree.height.diff + nodeheight(phylo2, (phylo2$edge[phylo2$edge[, 1] == i + totlen2, 1][1]))
        int[5] <- tree.height.diff + nodeheight(phylo2, b)
        int[6] <- b
        int[7] <- "tree2"
        mat2 <- rbind(mat2, int)
      }
    }
    mat2 <- mat2[match(phylo2$edge[,2],mat2[,6]),]
    first.name <- paste(".", THREELETTERS[last.count+1],sep="")
        root.branch <- as.data.frame(t(c(0, first.name, Ntip(phylo2)+1, 0, tree.height.diff, Ntip(phylo2)+1, "tree2")))
    mat2 <- rbind(root.branch, mat2)    
    
    
    # [,1] is the parent node number
    # [,2] are the daughter nodes or tips by name
    # [,3] are the daughter nodes (#>0) or tips (0)
    # [,4] are the starting ages of the branches (above the root)
    # [,5] are the maximum ages of the branches (above the root)
    # [,6] are the original daughter nodes or tips
    # [,7] is the tree the tips/nodes came from (silly I know)
    mat <- as.data.frame(rbind(mat1, mat2))
        mat[,4] <- round(as.numeric(as.character(mat[,4])), digits=6)
           mat[,5] <- round(as.numeric(as.character(mat[,5])), digits=6)

### MY PROBLEMS ARE COMING FROM HERE (SUPRISE)
    if (is.data.frame(map)) {
      nat <- list()
      for (i in 1:(length(allnodes)-1)) {
        if (i == 1) {
          nat[[i]] <- list(mat[mat[,4]==0,2])
        } else {
          starts <- filter(mat, mat[,4]>=allnodes[1] & mat[,4]<=allnodes[i])
          ends <- starts[starts[,5]>allnodes[i],2]
          nat[[i]] <- list(ends)
        }
      }
    }
    
    int.names1 <- rownames(vcvPhylo(phylo1)) # this does not include a crown taxon representing all tips
    counter_three_letters <- 1
    for (p in (Ntip(phylo1)+1):length(int.names1)) {
      int.names1[p] <- paste(".", THREELETTERS[counter_three_letters], 
                             sep = "")
      counter_three_letters <- counter_three_letters+1
    } 
    int.names <- int.names1
    # CREATE ANOTHER LOOP FOR THE SECOND TREE
    int.names2 <- c(rownames(vcvPhylo(phylo2)), "root")
    for (p in (Ntip(phylo2)+1):length(int.names2)) {
      int.names2[p] <- paste(".", THREELETTERS[counter_three_letters], 
                             sep = "")
      counter_three_letters <- counter_three_letters+1
    } 
    # THEN APPEND THE NAMES TO THE int.names OBJECT
    int.names <- append(int.names, int.names2)
    
    all.sp1 <- list(); all.sp2 <- list()
    all.hull1 <- list(); all.hull2 <- list()
    
    # make a loop that creates distributions of all tips FOR THE FIRST TREE
    for (p in 1:Ntip(phylo1)) {
      current.taxon <- phylo1$tip.label[[p]]
      current.data <- subset(map, map[,1] == current.taxon)
      if (nrow(current.data) > 1000) {
        current.data <- current.data[sample(nrow(current.data), 1000, replace=F), ]
      } 
      points <- current.data[,c(2,3)]
      all.sp1[[p]] <- distribution.sp <- SpatialPoints(points)
      all.hull1[[p]] <- distribution.hull <- gBuffer(distribution.sp, width=point.width)
    }
    names(all.sp1) <- phylo1$tip.label
    names(all.hull1) <- phylo1$tip.label
    
    # make a loop that creates distributions at nodes by combining daughter distributions
    # this is only a rough way of getting ancestral ranges, if you haven't run RASE
    if (is.null(rase.obj1)) {
      for (kk in (length(int.names1)+1):(length(phylo1$tip.label)+2)){
        descendantz <- Descendants(phylo1, kk, type="children")
        combined.range <- gUnion(all.hull1[[descendantz[1]]], all.hull1[[descendantz[2]]])
        all.hull1[[kk]] <- combined.range
      }
      all.hull1[(Ntip(phylo1)+1):(length(int.names1))] <- all.hull1[(Ntip(phylo1)+2):(length(int.names1)+1)] # paste over the root node distribution (we don't need it for the older tree!)
      all.hull1 <- all.hull1[1:length(int.names1)] # remove the duplicate last distribution
      names(all.hull1)[(Ntip(phylo1)+1) : (length(all.hull1))] <- int.names1[(Ntip(phylo1)+1):(length(int.names1))]
    }
    
    if (!is.null(rase.obj1)) {
      name.list <- NULL
      for (j in 1:length(rase1$ConvexHulls)) {
        nnames <- mat1[which(names(rase1$ConvexHulls)[j] == paste0("n",mat1[,3])),2]
        name.list[[j]] <- as.character(nnames[which(!nnames %in% phylo1$tip.label)])
      }
      names(rase1$SpatialPoints) <- name.list; rase1$SpatialPoints[[1]] <- NULL # drop the root node distribution
          all.sp1 <- append(all.sp1, rase1$ConvexHulls)
      names(rase1$ConvexHulls) <- name.list; rase1$ConvexHulls[[1]] <- NULL # drop the root node distribution
          all.hull1 <- append(all.hull1, rase1$ConvexHulls)
    }
    
    
    # make a loop that creates distributions of all tips FOR THE SECOND TREE
    for (p in 1:Ntip(phylo2)) {
      current.taxon <- phylo2$tip.label[[p]]
      current.data <- subset(map, map[,1] == current.taxon)
      if (nrow(current.data) > 1000) {
        current.data <- current.data[sample(nrow(current.data), 1000, replace=F), ]
      } 
      points <- current.data[,c(2,3)]
      all.sp2[[p]] <- distribution.sp <- SpatialPoints(points)
      all.hull2[[p]] <- distribution.hull <- gBuffer(distribution.sp, width=point.width)
    }
    names(all.sp2) <- phylo2$tip.label
    names(all.hull2) <- phylo2$tip.label

    # make a loop that creates distributions at nodes by combining daughter distributions
    # this works backwards from nodes nearest the tips, back to the root (but doesn't do the root)
    # we can take advantage of this naming method for the RASE output
    if (is.null(rase.obj2)) {
      for (kk in (length(int.names2)):(length(phylo2$tip.label)+1)){
        descendantz <- Descendants(phylo2, kk, type="children")
        combined.range <- gUnion(all.hull2[[descendantz[1]]], all.hull2[[descendantz[2]]])
        all.hull2[[kk]] <- combined.range
      }
      names(all.hull2)[(Ntip(phylo2)+1) : (length(all.hull2))] <- int.names2[(Ntip(phylo2)+1):(length(int.names2))]
    }
    if (!is.null(rase.obj2)) {
      name.list <- NULL
      for (j in 1:length(rase2$ConvexHulls)) {
        nnames <- mat2[which(names(rase2$ConvexHulls)[j] == paste0("n",mat2[,3])),2]
        name.list[[j]] <- as.character(nnames[which(!nnames %in% phylo2$tip.label)])
      }
      names(rase2$SpatialPoints) <- name.list
          all.sp2 <- append(all.sp2, rase2$ConvexHulls)
      names(rase2$ConvexHulls) <- name.list
          all.hull2 <- append(all.hull2, rase2$ConvexHulls)
    }
    
    all.hull <- append(all.hull1, all.hull2)
    
    unique.combo <-  combn(int.names, m = 2)
    range_overlap <- rep(NA, ncol(unique.combo))
    arrange <- data.frame(species1 = unique.combo[1,], 
                          species2 = unique.combo[2,], 
                          range_overlap, stringsAsFactors=F)
    for (t in 1:length(arrange[,1])) {
      overlap <- gOverlaps(all.hull[[arrange[t,1]]], all.hull[[arrange[t,2]]])
      arrange[t,3] <- overlap
    }
    arrange[,3][arrange[,3] == "TRUE"] <- 1
    
    master.matrix <- matrix(nrow=length(int.names), ncol=length(int.names))
    master.matrix <- as.data.frame(master.matrix)

    ordered.names <- mat[,2]
    rownames(master.matrix) <- ordered.names
    colnames(master.matrix) <- ordered.names
    
    for (j in 1:length(ordered.names)) {
      for (k in 1:length(ordered.names)) {
        if (j == k) {
          master.matrix[j,k] <- 1
        } else {
          master.matrix[j,k] <- filter(arrange, (species1==ordered.names[j] & species2==ordered.names[k]) |
                                         (species2==ordered.names[j]) & species1==ordered.names[k])$range_overlap
        }
      }
    }
  }
  
    geography.matrix <- list()
    for (i in 1:length(nat)) { 
      var.list <- unlist(nat[[i]])
      int.mat <- master.matrix[which(rownames(master.matrix)%in%var.list),]
      int.mat <- int.mat[,which(colnames(int.mat)%in%var.list)]
      int.mat <- as.matrix(int.mat) 
      geography.matrix[[i]] <- int.mat
    }
    
  end <- Sys.time()
  duration <- format(end-beginning)
  print(paste("Computation time :", duration))
  
  return(list(geography.object = geography.matrix, times = allnodes, spans = diff(allnodes),
              startingTimes = mat[,4], endTimes=mat[,5]))
}

#test <- CreateCoEvoGeoObject_SP(tree_1, tree_2, distribution)

# KEEPING THIS HERE TO TEST, BUT NEEDS TO BE MOVED TO A NEW FUNCTION
        #names(rase.out$ConvexHulls)
        #testo <- rase.out$ConvexHulls[2:length(rase.out$ConvexHulls)]
        #names(testo) <- int.names1[(Ntip(tree_1)+1) : length(int.names1)]
# if trees are uneven, we want to drop the first convex hull/spatial data from the
# older tree, if they're they same age, I'm not sure. 
# now I have to think about how the second one is going to work.
# the function would need to take in two RASE outputs, unless I can combine em.

# the nodes for each tree from rase are labelled "n9", "n10", ...
# the nodes in the above are just "9", "10", ...
# luckily we should be able to use "n9" == paste0("n", 9) to our advantage!
# mat1[,6]
# I want to now bring in two rase objects and be able to call the node distributions
# from the rase objects, so let's take it from there.
