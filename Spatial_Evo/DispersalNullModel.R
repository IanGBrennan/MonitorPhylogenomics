# Big thanks to Alex Skeels (https://alexskeels.com/) for help here!

DNM <- function (picante.cdm, tree, distances.among, abundance.matters = FALSE, 
          abundance.assigned = "directly") 
{
  tempCheck <- checkCDM(picante.cdm)
  if (tempCheck == "fail") {
    stop("CDM incompatible with dispersalNull model. See 'checkCDM' for details")
  }
  #if (length(setdiff(names(picante.cdm), tree$tip.label)) != 
  #    0) {
  #  stop("You have included species in your cdm that are not in your phylogeny")
  #}
  if (length(setdiff(row.names(picante.cdm), row.names(distances.among))) != 
      0 & length(setdiff(row.names(distances.among), row.names(picante.cdm))) != 
      0) {
    stop("Your cdm plot names and distance matrix names do no match")
  }
  if (any(row.names(picante.cdm) != row.names(distances.among))) {
    stop("Your cdm and distance matrix are not in the same plot order")
  }
  richness <- apply(picante.cdm, 1, lengthNonZeros)
  overallAbundance <- picante.cdm[picante.cdm != 0]
  replacementList <- list()
  for (i in 1:dim(picante.cdm)[1]) {
    phylocom <- matrix(ncol = 3, nrow = richness[i], 0)
    phylocom <- as.data.frame(phylocom)
    names(phylocom) <- c("plot", "abund", "id")
    j <- 0
    while (length(phylocom[phylocom$plot == row.names(picante.cdm)[i], 
                           ]$id) < richness[i]) {
      selectedPlot <- selectNear(distances.among[, i])
      if (abundance.matters) {
        temp <- sample(x = picante.cdm[selectedPlot, 
                                       ], size = 1, prob = picante.cdm[selectedPlot, 
                                                                       ])
      }
      else {
        possible <- picante.cdm[selectedPlot, ][picante.cdm[selectedPlot, 
                                                            ] != 0]
        names(possible) <- names(picante.cdm[selectedPlot, 
                                             ])[picante.cdm[selectedPlot, ] != 0]
        
        if(length(possible) == 1){ temp <- possible} else { temp <- sample(x = possible, size = 1)}
        
        
      }
      if (!(names(temp) %in% phylocom[phylocom$plot == 
                                      row.names(picante.cdm)[i], ]$id)) {
        j <- j + 1
        phylocom[j, 1] <- row.names(picante.cdm)[i]
        if (abundance.assigned == "directly") {
          phylocom[j, 2] <- temp
        }
        else if (abundance.assigned == "explore") {
          distribution <- round(rnorm(n = 100, mean = as.numeric(temp), 
                                      sd = 1))
          distribution[distribution < 0] <- 1
          chosen <- sample(distribution, 1)
          phylocom[j, 2] <- chosen
        }
        else if (abundance.assigned == "overall") {
          chosen <- sample(overallAbundance, 1)
          phylocom[j, 2] <- chosen
        }
        else {
          stop("abundance.assigned argument set to unrecognized value")
        }
        phylocom[j, 3] <- names(temp)
      }
    }
    replacementList[[i]] <- phylocom
  }
  newCDM <- Reduce(rbind, replacementList)
  newCDM <- sample2matrix(newCDM); print(c(ncol(newCDM), nrow(newCDM)))
  notFound <- setdiff(colnames(picante.cdm), names(newCDM))
  print(c(length(colnames(picante.cdm)), length(names(newCDM))))
  if (length(notFound > 0)) {
    toBind <- matrix(nrow = dim(newCDM)[[1]], ncol = length(notFound), 0)
    toBind[1,] <- 1 # this is because otherwise the taxa with no occurrences will ruin the dbFD function
    colnames(toBind) <- notFound
    newCDM <- cbind(newCDM, toBind)
  }
  print(c(ncol(newCDM), nrow(newCDM)))
  newCDM <- as.matrix(newCDM)
  newCDM <- newCDM[row.names(picante.cdm), ]
  newCDM <- newCDM[ , colnames(picante.cdm)]
  return(newCDM)
}


selectNear <- function(distances.between)
{
  #distances.between is a vector of distances between the focal cell and other cells
  #first exclude distances to the focal cell (and any other with distance = 0)
  distances.between <- distances.between[distances.between != 0]
  
  #now sample a cell name with a probability proportional to the inverse of the distance
  #from the focal cell
  newCell <- sample(x=names(distances.between), size=1, prob=1/distances.between)
  
  newCell
}
