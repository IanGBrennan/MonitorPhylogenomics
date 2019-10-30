####################################################
#    This file contains additional RPANDA functions
#    either from H.Morlon's GitHub, or I made up
####################################################

# Russell suggested that I could use the GenSA function instead of optim, and it might get a better model fit (though slower)
# I'd have to dig out where the actual optimizing function is first!

# FOR COEVOLUTION MODELS I PROBABLY NEED TO CORRECT THE BLOCKS IF TREE1 KEEPS BRANCHING BUT TREE2 HASN'T BRANCHED AT ALL
# THIS WOULD CAUSE ISSUES WITH THE BEGINNING MATRICES (SEE CoEvo_Split, et al.)

# I ALSO NEED TO CORRECT THE DIVISION OF S BY OVERLAPPING LINEAGES, NEED TO INCLUDE THE FOCAL TAXON (NOT MAKE DIAG(0) TO START)
# BECAUSE OF THIS AND THE LACK OF SYMMETRY, THE COEVO_SPLIT MODEL STILL ISN'T CORRECT!

createModel <- function(tree, keyword){
  
  if(keyword == "BM" || keyword == "BMbis"){
    
    comment <- "Brownian Motion model with linear drift.\nStarts with two lineages having the same value X_0 ~ Normal(m0,v0).\nOne trait in each lineage, all lineages evolving independently after branching.\ndX_t = d dt + sigma dW_t"
    paramsNames <- c("m0", "v0", "d", "sigma")
    params0 <- c(0,0,0,1)
    
    periodizing <- periodizeOneTree(tree)
    eventEndOfPeriods <- endOfPeriods(periodizing, tree)
    
    initialCondition <- function(params) return( list(mean=c(params[1]), var=matrix(c(params[2]))) )
    
    aAGamma <- function(i, params){
      vectorU <- getLivingLineages(i, eventEndOfPeriods)
      vectorA <- function(t) return(params[3]*vectorU)
      matrixGamma <- function(t) return(params[4]*diag(vectorU))
      matrixA <- diag(0, length(vectorU))
      
      return(list(a=vectorA, A=matrixA, Gamma=matrixGamma))
    }
    
    constraints <- function(params) return(params[2]>=0 && params[4]>=0)
    
    if( keyword == "BM" ){
      model <- new(Class="PhenotypicBM", name=keyword, period=periodizing$periods, aAGamma=aAGamma, numbersCopy=eventEndOfPeriods$copy, numbersPaste=eventEndOfPeriods$paste, initialCondition=initialCondition, paramsNames=paramsNames, constraints=constraints, params0=params0, tipLabels=tree$tip.label, tipLabelsSimu=eventEndOfPeriods$labeling, comment=comment, matrixCoalescenceTimes=findMRCA(tree, type="height"))
    }else{
      model <- new(Class="PhenotypicModel", name=keyword, period=periodizing$periods, aAGamma=aAGamma, numbersCopy=eventEndOfPeriods$copy, numbersPaste=eventEndOfPeriods$paste, initialCondition=initialCondition, paramsNames=paramsNames, constraints=constraints, params0=params0, tipLabels=eventEndOfPeriods$labeling, tipLabelsSimu=eventEndOfPeriods$labeling, comment=comment)
    }
    
  }
  else if(keyword == "BM_from0"){
    
    comment <- "Brownian Motion model with linear drift.\nStarts with two lineages having the same value X_0 = (0,0).\nOne trait in each lineage, all lineages evolving independently after branching.\ndX_t = d dt + sigma dW_t"
    paramsNames <- c("d", "sigma")
    params0 <- c(0,1)
    
    periodizing <- periodizeOneTree(tree)
    eventEndOfPeriods <- endOfPeriods(periodizing, tree)
    
    initialCondition <- function(params) return( list( mean=c(0), var=matrix(c(0)) ) )
    
    aAGamma <- function(i, params){
      vectorU <- getLivingLineages(i, eventEndOfPeriods)
      vectorA <- function(t) return(params[1]*vectorU)
      matrixGamma <- function(t) return(params[2]*diag(vectorU))
      matrixA <- diag(0, length(vectorU))
      
      return(list(a=vectorA, A=matrixA, Gamma=matrixGamma))
    }
    
    constraints <- function(params) return(params[2]>=0)
    
    model <- new(Class="PhenotypicBM", name=keyword, period=periodizing$periods, aAGamma=aAGamma, numbersCopy=eventEndOfPeriods$copy, numbersPaste=eventEndOfPeriods$paste, initialCondition=initialCondition, paramsNames=paramsNames, constraints=constraints, params0=params0, tipLabels=tree$tip.label, tipLabelsSimu=eventEndOfPeriods$labeling, comment=comment, matrixCoalescenceTimes=findMRCA(tree, type="height"))
    
  }
  else if(keyword == "BM_from0_driftless"){
    
    comment <- "Brownian Motion model without drift.\nStarts with two lineages having the same value X_0 = (0,0).\nOne trait in each lineage, all lineages evolving independently after branching.\ndX_t = sigma dW_t"
    paramsNames <- c("sigma")
    params0 <- c(1)
    
    periodizing <- periodizeOneTree(tree)
    eventEndOfPeriods <- endOfPeriods(periodizing, tree)
    
    initialCondition <- function(params) return( list( mean=c(0), var=matrix(c(0)) ) )
    
    aAGamma <- function(i, params){
      vectorU <- getLivingLineages(i, eventEndOfPeriods)
      vectorA <- function(t) return(0*vectorU)
      matrixGamma <- function(t) return(params[1]*diag(vectorU))
      matrixA <- diag(0, length(vectorU))
      
      return(list(a=vectorA, A=matrixA, Gamma=matrixGamma))
    }
    
    constraints <- function(params) return(params[1]>=0)
    
    model <- new(Class="PhenotypicBM", name=keyword, period=periodizing$periods, aAGamma=aAGamma, numbersCopy=eventEndOfPeriods$copy, numbersPaste=eventEndOfPeriods$paste, initialCondition=initialCondition, paramsNames=paramsNames, constraints=constraints, params0=params0, tipLabels=tree$tip.label, tipLabelsSimu=eventEndOfPeriods$labeling, comment=comment, matrixCoalescenceTimes=findMRCA(tree, type="height"))
    
  }
  else if(keyword == "OU" || keyword == "OUbis" || keyword == "OUter"){
    
    comment <- "Ornstein-Uhlenbeck model.\nStarts with two lineages having the same value X_0 ~ Normal(m0,v0).\nOne trait in each lineage, all lineages evolving independently after branching.\ndX_t = psi(theta- X_t) dt + sigma dW_t"
    paramsNames <- c("m0", "v0", "psi", "theta", "sigma")
    params0 <- c(0,0,1,0,1)
    
    periodizing <- periodizeOneTree(tree)
    eventEndOfPeriods <- endOfPeriods(periodizing, tree)
    
    initialCondition <- function(params) return( list(mean=c(params[1]), var=matrix(c(params[2]))) )
    
    aAGamma <- function(i, params){
      vectorU <- getLivingLineages(i, eventEndOfPeriods)
      vectorA <- function(t) return(params[3]*params[4]*vectorU)
      matrixGamma <- function(t) return(params[5]*diag(vectorU))
      matrixA <- params[3]*diag(vectorU)
      
      return(list(a=vectorA, A=matrixA, Gamma=matrixGamma))
    }
    
    constraints <- function(params) return(params[2]>=0 && params[5]>=0 && params[3]!=0)
    
    if( keyword == "OU" ){
      model <- new(Class="PhenotypicOU", name=keyword, period=periodizing$periods, aAGamma=aAGamma, numbersCopy=eventEndOfPeriods$copy, numbersPaste=eventEndOfPeriods$paste, initialCondition=initialCondition, paramsNames=paramsNames, constraints=constraints, params0=params0, tipLabels=tree$tip.label, tipLabelsSimu=eventEndOfPeriods$labeling, comment=comment, matrixCoalescenceTimes=findMRCA(tree, type="height"))
    }else if( keyword == "OUbis" ){
      model <- new(Class="PhenotypicADiag", name=keyword, period=periodizing$periods, aAGamma=aAGamma, numbersCopy=eventEndOfPeriods$copy, numbersPaste=eventEndOfPeriods$paste, initialCondition=initialCondition, paramsNames=paramsNames, constraints=constraints, params0=params0, tipLabels=eventEndOfPeriods$labeling, tipLabelsSimu=eventEndOfPeriods$labeling, comment=comment)
    }else{
      model <- new(Class="PhenotypicModel", name=keyword, period=periodizing$periods, aAGamma=aAGamma, numbersCopy=eventEndOfPeriods$copy, numbersPaste=eventEndOfPeriods$paste, initialCondition=initialCondition, paramsNames=paramsNames, constraints=constraints, params0=params0, tipLabels=eventEndOfPeriods$labeling, tipLabelsSimu=eventEndOfPeriods$labeling, comment=comment)
    }
    
  }
  else if(keyword == "OU_from0"){
    
    comment <- "Ornstein-Uhlenbeck model.\nStarts with two lineages having the same value X_0 = (0,0).\nOne trait in each lineage, all lineages evolving independently after branching.\ndX_t = psi(theta- X_t) dt + sigma dW_t"
    paramsNames <- c("psi", "theta", "sigma")
    params0 <- c(0.01,0,1)
    # This model requires the following list of parameters, in this order : psi, theta, sigma
    # One trait in each lineage
    # starts with two lineages with the same value X_0 ~ Normal(0,0)
    # dX_t = psi(theta- X_t) dt + sigma dW_t
    # Tree is an object of type "phylo" (cf. ape)
    
    periodizing <- periodizeOneTree(tree)
    eventEndOfPeriods <- endOfPeriods(periodizing, tree)
    
    initialCondition <- function(params) return( list(mean=c(0), var=matrix(c(0))) )
    
    aAGamma <- function(i, params){
      vectorU <- getLivingLineages(i, eventEndOfPeriods)
      vectorA <- function(t) return(params[1]*params[2]*vectorU)
      matrixGamma <- function(t) return(params[3]*diag(vectorU))
      matrixA <- params[1]*diag(vectorU)
      
      return(list(a=vectorA, A=matrixA, Gamma=matrixGamma))
    }
    
    constraints <- function(params) return(params[3]>=0 && params[1]!=0)
    
    model <- new(Class="PhenotypicOU", name=keyword, period=periodizing$periods, aAGamma=aAGamma, numbersCopy=eventEndOfPeriods$copy, numbersPaste=eventEndOfPeriods$paste, initialCondition=initialCondition, paramsNames=paramsNames, constraints=constraints, params0=params0, tipLabels=tree$tip.label, tipLabelsSimu=eventEndOfPeriods$labeling, comment=comment, matrixCoalescenceTimes=findMRCA(tree, type="height"))
    
  }
  else if(keyword == "ACDC" || keyword == "ACDCbis"){
    
    comment <- "ACcelerating or DeCelerating rate of evolution.\nStarts with two lineages having the same value X_0 ~ Normal(m0,v0).\nOne trait in each lineage, all lineages evolving independently after branching.\ndX_t = sigma0 exp(rt) dW_t"
    paramsNames <- c("m0", "v0", "sigma0", "r")
    params0 <- c(0,0,100,1)
    
    periodizing <- periodizeOneTree(tree)
    eventEndOfPeriods <- endOfPeriods(periodizing, tree)
    
    initialCondition <- function(params) return( list(mean=c(params[1]), var=matrix(c(params[2]))) )
    
    aAGamma <- function(i, params){
      vectorU <- getLivingLineages(i, eventEndOfPeriods)
      vectorA <- function(t) return(rep(0, length(vectorU)))
      matrixGamma <- function(t) return(params[3]*exp(params[4]*t)*diag(vectorU))
      matrixA <- diag(0, length(vectorU))
      
      return(list(a=vectorA, A=matrixA, Gamma=matrixGamma))
    }
    
    constraints <- function(params) return(params[2]>=0 && params[3]>0)
    
    if( keyword == "ACDC" ){
      model <- new(Class="PhenotypicACDC", name=keyword, period=periodizing$periods, aAGamma=aAGamma, numbersCopy=eventEndOfPeriods$copy, numbersPaste=eventEndOfPeriods$paste, initialCondition=initialCondition, paramsNames=paramsNames, constraints=constraints, params0=params0, tipLabels=tree$tip.label, tipLabelsSimu=eventEndOfPeriods$labeling, comment=comment, matrixCoalescenceTimes=findMRCA(tree, type="height"))
    }else{
      model <- new(Class="PhenotypicModel", name=keyword, period=periodizing$periods, aAGamma=aAGamma, numbersCopy=eventEndOfPeriods$copy, numbersPaste=eventEndOfPeriods$paste, initialCondition=initialCondition, paramsNames=paramsNames, constraints=constraints, params0=params0, tipLabels=eventEndOfPeriods$labeling, tipLabelsSimu=eventEndOfPeriods$labeling, comment=comment)
    }
    class(model)[1] <- "PhenotypicModel"
  }
  else if(keyword == "DD" || keyword == "DDbis"){
    
    comment <- "Diversity-dependent model.\nStarts with two lineages having the same value X_0 ~ Normal(m0,v0).\nOne trait in each lineage, all lineages evolving independently after branching.\ndX_t = sigma0 exp(r n_t) dW_t"
    paramsNames <- c("m0", "v0", "r", "sigma0")
    params0 <- c(0,0,100,1)
    
    periodizing <- periodizeOneTree(tree)
    eventEndOfPeriods <- endOfPeriods(periodizing, tree)
    
    initialCondition <- function(params) return( list(mean=c(params[1]), var=matrix(c(params[2]))) )
    
    aAGamma <- function(i, params){
      vectorU <- getLivingLineages(i, eventEndOfPeriods)
      vectorA <- function(t) return(rep(0, length(vectorU)))
      matrixGamma <- function(t) return(params[4]*exp(params[4]*sum(vectorU))*diag(vectorU))
      matrixA <- diag(0, length(vectorU))
      
      return(list(a=vectorA, A=matrixA, Gamma=matrixGamma))
    }
    
    constraints <- function(params) return(params[2]>=0 && params[4]>0)
    
    if( keyword == "DD" ){
      model <- new(Class="PhenotypicDD", name=keyword, period=periodizing$periods, aAGamma=aAGamma, numbersCopy=eventEndOfPeriods$copy, numbersPaste=eventEndOfPeriods$paste, initialCondition=initialCondition, paramsNames=paramsNames, constraints=constraints, params0=params0, tipLabels=tree$tip.label, tipLabelsSimu=eventEndOfPeriods$labeling, comment=comment, matrixCoalescenceJ=getMatrixCoalescenceJ(tree, periodizing$periods), nLivingLineages=eventEndOfPeriods$nLivingLineages)
    }else{
      model <- new(Class="PhenotypicModel", name=keyword, period=periodizing$periods, aAGamma=aAGamma, numbersCopy=eventEndOfPeriods$copy, numbersPaste=eventEndOfPeriods$paste, initialCondition=initialCondition, paramsNames=paramsNames, constraints=constraints, params0=params0, tipLabels=eventEndOfPeriods$labeling, tipLabelsSimu=eventEndOfPeriods$labeling, comment=comment)
    }
    class(model)[1] <- "PhenotypicModel"
  }
  else if(keyword == "PM" || keyword == "PMbis" || keyword == "PMter"){
    
    comment <- "Phenotype Matching model.\nStarts with two lineages having the same value X_0 ~ Normal(m0,v0).\nOne trait in each lineage, all lineages evolving then non-independtly according to the Phenotype Matching expression."
    print("PM: Phenotypic Matching model")
    paramsNames <- c("m0", "v0", "theta", "psi", "S", "sigma")
    params0 <- c(0,0,0,0.2,0.5,1)
    
    periodizing <- periodizeOneTree(tree)
    eventEndOfPeriods <- endOfPeriods(periodizing, tree)
    
    initialCondition <- function(params) return( list(mean=c(params[1]), var=matrix(c(params[2]))) ) 
    
    aAGamma <- function(i, params){
      vectorU <- getLivingLineages(i, eventEndOfPeriods)
      vectorA <- function(t) return(params[3]*params[4]*vectorU)
      matrixGamma <- function(t) return(params[6]*diag(vectorU))
      matrixA <- (params[4]+params[5])*diag(vectorU) - (params[5]/(sum(vectorU)-1)) * outer(vectorU,vectorU) 
      # above: matrix of('attraction towards optimum' + 'repulsion from others') - matrix of('repulsion' / # of lineages)
      
      return(list(a=vectorA, A=matrixA, Gamma=matrixGamma, u=vectorU, OU=TRUE))
    }
    
    constraints <- function(params) return(params[2]>=0 && params[6]>=0)
    
    if( keyword == "PM" ){
      model <- new(Class="PhenotypicPM", name=keyword, period=periodizing$periods, aAGamma=aAGamma, numbersCopy=eventEndOfPeriods$copy, numbersPaste=eventEndOfPeriods$paste, initialCondition=initialCondition, paramsNames=paramsNames, constraints=constraints, params0=params0, tipLabels=eventEndOfPeriods$labeling, tipLabelsSimu=eventEndOfPeriods$labeling, comment=comment)
    }else if( keyword == "PMbis" ){
      model <- new(Class="PhenotypicADiag", name=keyword, period=periodizing$periods, aAGamma=aAGamma, numbersCopy=eventEndOfPeriods$copy, numbersPaste=eventEndOfPeriods$paste, initialCondition=initialCondition, paramsNames=paramsNames, constraints=constraints, params0=params0, tipLabels=eventEndOfPeriods$labeling, tipLabelsSimu=eventEndOfPeriods$labeling, comment=comment)
    }else{
      model <- new(Class="PhenotypicModel", name=keyword, period=periodizing$periods, aAGamma=aAGamma, numbersCopy=eventEndOfPeriods$copy, numbersPaste=eventEndOfPeriods$paste, initialCondition=initialCondition, paramsNames=paramsNames, constraints=constraints, params0=params0, tipLabels=eventEndOfPeriods$labeling, tipLabelsSimu=eventEndOfPeriods$labeling, comment=comment)
    }
    
  } 
  else if(keyword == "MC") {
    
    comment <- "Matching Competition model\n Implemented as in Drury et al. Systematic Biology."
    paramsNames <- c("m0","logsigma","S")
    params0 <- c(0,1,0.1)
    
    periodizing <- periodizeOneTree(tree) 
    eventEndOfPeriods <- endOfPeriods(periodizing, tree)
    
    initialCondition <- function(params) return( list(mean=c(params[1]), var=matrix(c(0))) ) 
    
    ###is this where the A matrix incorporating geography needs to go? if so, what is the order in which lineage sympatry data need to be introduced
    
    aAGamma <- function(i, params){
      vectorU <- getLivingLineages(i, eventEndOfPeriods)
      vectorA <- function(t) return(0*vectorU)
      matrixGamma <- function(t) return(params[2]*diag(vectorU))
      matrixA <- params[3]*diag(vectorU) - (params[3]/sum(vectorU)) * outer(vectorU,vectorU) 
      return(list(a=vectorA, A=matrixA, Gamma=matrixGamma, u=vectorU, OU=FALSE))
    }
    constraints <- function(params) return(params[3]<=0)
    model <- new(Class="PhenotypicADiag", name=keyword, period=periodizing$periods, aAGamma=aAGamma, numbersCopy=eventEndOfPeriods$copy, numbersPaste=eventEndOfPeriods$paste, initialCondition=initialCondition, paramsNames=paramsNames, constraints=constraints, params0=params0, tipLabels=eventEndOfPeriods$labeling,  comment=comment)

  } 
  else if(keyword == "PM_OUless" || keyword == "PM_OUlessbis"){
    
    comment <- "Simplified Phenotype Matching model.\nStarts with two lineages having the same value X_0 ~ Normal(m0,v0).\nOne trait in each lineage, all lineages evolving then non-independtly according to the Phenotype Matching expression, without the OU term."
    print("PM_OUless: fit the simplified version of the Phenotypic Matching model (should be nearly identical to the MC model)")
    paramsNames <- c("m0", "v0", "S", "sigma")
    params0 <- c(0,0,0.5,1)
    
    periodizing <- periodizeOneTree(tree)
    eventEndOfPeriods <- endOfPeriods(periodizing, tree)
    
    initialCondition <- function(params) return( list(mean=c(params[1]), var=matrix(c(params[2]))) ) 
    
    aAGamma <- function(i, params){
      vectorU <- getLivingLineages(i, eventEndOfPeriods)
      vectorA <- function(t) return(0*vectorU)
      matrixGamma <- function(t) return(params[4]*diag(vectorU))
      matrixA <- params[3]*diag(vectorU) - (params[3]/sum(vectorU)) * outer(vectorU,vectorU) 
      
      return(list(a=vectorA, A=matrixA, Gamma=matrixGamma, u=vectorU, OU=FALSE))
    }
    
    constraints <- function(params) return(params[2]>=0 && params[4]>=0)
    
    if( keyword == "PM_OUless" ){
      model <- new(Class="PhenotypicPM", name=keyword, period=periodizing$periods, aAGamma=aAGamma, numbersCopy=eventEndOfPeriods$copy, numbersPaste=eventEndOfPeriods$paste, initialCondition=initialCondition, paramsNames=paramsNames, constraints=constraints, params0=params0, tipLabels=eventEndOfPeriods$labeling, tipLabelsSimu=eventEndOfPeriods$labeling, comment=comment)
    }else{
      model <- new(Class="PhenotypicADiag", name=keyword, period=periodizing$periods, aAGamma=aAGamma, numbersCopy=eventEndOfPeriods$copy, numbersPaste=eventEndOfPeriods$paste, initialCondition=initialCondition, paramsNames=paramsNames, constraints=constraints, params0=params0, tipLabels=eventEndOfPeriods$labeling, tipLabelsSimu=eventEndOfPeriods$labeling, comment=comment)
    }
    
  }
  else{
    stop("Keyword does not correspond to any model in the model bank")
  }
  
  
  return(model)
}

createGeoModel <- function(tree, geo.object, keyword){
  
  if(keyword == "MC" || keyword == "MC+geo"){
    
    comment <- "Matching competition model with biogeography\n Implemented as in Drury et al. Systematic Biology."
    paramsNames <- c("m0","logsigma","S")
    params0 <- c(0,log(1),0)
    
    resgeo.object <- resortGeoObject(tree, geo.object)
    periodizing <- periodizeOneTree_geo(tree,resgeo.object) 
    eventEndOfPeriods <- endOfPeriods(periodizing, tree)
    
    initialCondition <- function(params) return( list(mean=c(params[1]), var=matrix(c(0))) ) 
    
    ###is this where the A matrix incorporating geography needs to go? if so, what is the order in which lineage sympatry data need to be introduced
    
    aAGamma <- function(i, params){
      vectorU <- getLivingLineages(i, eventEndOfPeriods)
      vectorA <- function(t) return(0*vectorU)
      matrixGamma <- function(t) return(exp(params[2])*diag(vectorU))
      nij <- colSums(resgeo.object$geography.object[[i]])
      matrixA <- params[3]*diag(vectorU) -(resgeo.object$geography.object[[i]]*(params[3]/nij))
      return(list(a=vectorA, A=matrixA, Gamma=matrixGamma))
    }
    constraints <- function(params) return(params[3]<=0)
    model <- new(Class="PhenotypicModel", name="MC+geo", period=periodizing$periods, aAGamma=aAGamma, numbersCopy=eventEndOfPeriods$copy, numbersPaste=eventEndOfPeriods$paste, initialCondition=initialCondition, paramsNames=paramsNames, constraints=constraints, params0=params0, tipLabels=eventEndOfPeriods$labeling,  comment=comment)
    return(model)
  }
  
  else if(keyword == "PM+geo"){
    comment <- "Phenotype Matching model with biogeography.\nStarts with two lineages having the same value X_0 ~ Normal(m0,v0).\nOne trait in each lineage, all lineages evolving then non-independtly according to the Phenotype Matching expression."
    print("PM+geo: Phenotypic Matching model with biogeography.")
    paramsNames <- c("m0", "v0", "theta", "psi", "S", "sigma")
    params0 <- c(0,0,0,0.2,0.5,1)
    
    resgeo.object <- resortGeoObject(tree, geo.object)
    periodizing <- periodizeOneTree_geo(tree, resgeo.object)
    eventEndOfPeriods <- endOfPeriods(periodizing, tree)
    
    initialCondition <- function(params) return( list(mean=c(params[1]), var=matrix(c(params[2]))) ) 
    
    aAGamma <- function(i, params){
      vectorU <- getLivingLineages(i, eventEndOfPeriods)
      vectorA <- function(t) return(params[3]*params[4]*vectorU)
      matrixGamma <- function(t) return(params[6]*diag(vectorU))
      #nij <- colSums(resgeo.object$geography.object[[i]])
      # matrixA <- (params[4]+params[5])*diag(vectorU) - (params[5]/sum(vectorU)) * outer(vectorU,vectorU) # from the PM model
      # matrixA <- params[3]*diag(vectorU) -(geo.object$geography.object[[i]]*(params[3]/nij)) # from the MC_geo model
      # matrixA <- (params[4]+params[5])*diag(vectorU) - (params[5]/nij) * resgeo.object$geography.object[[i]] # second attempt
      matrixA <- (params[4]+params[5]) * diag(vectorU) - (params[5]/sum(vectorU)) * resgeo.object$geography.object[[i]] # third attempt
      
      return(list(a=vectorA, A=matrixA, Gamma=matrixGamma, u=vectorU, OU=TRUE))
    }
    
    constraints <- function(params) return(params[2]>=0 && params[6]>=0)
    
    if( keyword == "PM" ){
      model <- new(Class="PhenotypicModel", name="PM+geo", period=periodizing$periods, aAGamma=aAGamma, numbersCopy=eventEndOfPeriods$copy, numbersPaste=eventEndOfPeriods$paste, initialCondition=initialCondition, paramsNames=paramsNames, constraints=constraints, params0=params0, tipLabels=eventEndOfPeriods$labeling, tipLabelsSimu=eventEndOfPeriods$labeling, comment=comment)
    }else{
      model <- new(Class="PhenotypicModel", name=keyword, period=periodizing$periods, aAGamma=aAGamma, numbersCopy=eventEndOfPeriods$copy, numbersPaste=eventEndOfPeriods$paste, initialCondition=initialCondition, paramsNames=paramsNames, constraints=constraints, params0=params0, tipLabels=eventEndOfPeriods$labeling, tipLabelsSimu=eventEndOfPeriods$labeling, comment=comment)
    }
  } 
  
  else if(keyword == "PMOU+geo"){
    
    comment <- "Simplified Phenotype Matching model.\nStarts with two lineages having the same value X_0 ~ Normal(m0,v0).\nOne trait in each lineage, all lineages evolving then non-independtly according to the Phenotype Matching expression, without the OU term."
    print("PMOUless+geo: fit the simplified version of the Phenotypic Matching model accounting for biogeography. This correctly estimates S from only distributionally overlapping taxa")
    paramsNames <- c("m0", "v0", "S", "sigma")
    params0 <- c(0,0,0.5,1)
    
    resgeo.object <- resortGeoObject(tree, geo.object)
    periodizing <- periodizeOneTree_geo(tree, resgeo.object)
    eventEndOfPeriods <- endOfPeriods(periodizing, tree)
    
    initialCondition <- function(params) return( list(mean=c(params[1]), var=matrix(c(params[2]))) ) 
    
    aAGamma <- function(i, params){
      vectorU <- getLivingLineages(i, eventEndOfPeriods)
      vectorA <- function(t) return(0*vectorU)
      matrixGamma <- function(t) return(params[4]*diag(vectorU))
      # matrixA <- params[3]*diag(vectorU) - (params[3]/sum(vectorU)) * outer(vectorU,vectorU) 
      # matrixA <- (params[4]+params[5]) * diag(vectorU) - (params[5]/sum(vectorU)) * resgeo.object$geography.object[[i]] # third attempt
      nij <- colSums(resgeo.object$geography.object[[i]])
      matrixA <- params[3]*diag(vectorU) -(resgeo.object$geography.object[[i]]*(params[3]/nij))
      
      return(list(a=vectorA, A=matrixA, Gamma=matrixGamma, u=vectorU, OU=FALSE))
    }
    
    constraints <- function(params) return(params[2]>=0 && params[4]>=0)
    
    if( keyword == "PMOU+geo" ){
      model <- new(Class="PhenotypicPM", name=keyword, period=periodizing$periods, aAGamma=aAGamma, numbersCopy=eventEndOfPeriods$copy, numbersPaste=eventEndOfPeriods$paste, initialCondition=initialCondition, paramsNames=paramsNames, constraints=constraints, params0=params0, tipLabels=eventEndOfPeriods$labeling, tipLabelsSimu=eventEndOfPeriods$labeling, comment=comment)
    }else{
      model <- new(Class="PhenotypicADiag", name=keyword, period=periodizing$periods, aAGamma=aAGamma, numbersCopy=eventEndOfPeriods$copy, numbersPaste=eventEndOfPeriods$paste, initialCondition=initialCondition, paramsNames=paramsNames, constraints=constraints, params0=params0, tipLabels=eventEndOfPeriods$labeling, tipLabelsSimu=eventEndOfPeriods$labeling, comment=comment)
    }
    
  }
  
  else{
    stop("Keyword does not correspond to any model in the model bank")
  }
}

createModelCoevolution <- function(tree.1, tree.2, geo.object=NULL, keyword = "GMM"){
  
  # make sure we have the trees in the right order (oldest first), it matters for the matrixA construction (aAGamma$A)
  if (max(nodeHeights(tree.1))>=max(nodeHeights(tree.2))) { tree1 <- tree.1; tree2 <- tree.2 } 
  else if (max(nodeHeights(tree.2))>max(nodeHeights(tree.1))) { tree1 <- tree.2; tree2 <- tree.1 } 
  
  # this is the list of worthwhile models (working, mathematically tractable, or biologically sensible)
  if(keyword == "GMM" || keyword == "GMMbis"){
    if (!is.null(geo.object)) { print("ignoring geo.object, using GMM instead")}
    
    comment <- "The Generalist Matching Mutualism model. Assumes equal interaction (S) between all inter-clade lineages, but no interaction (0) among lineages within a tree (intra-clade)"
    print(comment)
    #comment <- "Generalist Matching Mutualism model.\nStarts with 3 or 4 lineages having the same value X_0 ~ Normal(m0,v0).\nOne trait in each lineage, all lineages evolving then non-independtly according to the GMM expression."
    paramsNames <- c("m0", "v0", "d1", "d2", "S", "sigma")
    params0 <- c(0,0,1,-1,0.5,1)
    
    eventEndOfPeriods <- endOfPeriodsGMM(tree1, tree2)
    n <- eventEndOfPeriods$nLineages1[1] + eventEndOfPeriods$nLineages2[1] - 1
    
    initialCondition <- function(params) return( list(mean=rep(params[1], times=n), var=matrix(rep(params[2], times=n*n), nrow=n ) ) ) 
    
    aAGamma <- function(i, params){
      vectorA <- function(t) return( c( rep(params[3]*params[5], times=eventEndOfPeriods$nLineages1[i]), rep(params[4]*params[5], times=eventEndOfPeriods$nLineages2[i]) ) )
      matrixGamma <- function(t) return(diag(params[6], eventEndOfPeriods$nLineages1[i] + eventEndOfPeriods$nLineages2[i]))
      
      bloc1 <- diag(params[5], eventEndOfPeriods$nLineages1[i])
      bloc2 <- matrix(rep(-params[5]/eventEndOfPeriods$nLineages2[i], times=eventEndOfPeriods$nLineages1[i]*eventEndOfPeriods$nLineages2[i]), nrow=eventEndOfPeriods$nLineages1[i])
      bloc3 <- matrix(rep(-params[5]/eventEndOfPeriods$nLineages1[i], times=eventEndOfPeriods$nLineages1[i]*eventEndOfPeriods$nLineages2[i]), nrow=eventEndOfPeriods$nLineages2[i])
      bloc4 <- diag(params[5], eventEndOfPeriods$nLineages2[i])
      matrixA <- rbind(cbind(bloc1, bloc2), cbind(bloc3, bloc4))
      
      return(list(a=vectorA, A=matrixA, Gamma=matrixGamma))
    } 
    
    constraints <- function(params) return(params[2]>=0 && params[6]>=0)
    
    if( keyword == "GMM" ){
      model <- new(Class="PhenotypicGMM", name=keyword, period=eventEndOfPeriods$periods, aAGamma=aAGamma, numbersCopy=eventEndOfPeriods$copy, numbersPaste=eventEndOfPeriods$paste, initialCondition=initialCondition, paramsNames=paramsNames, constraints=constraints, params0=params0, tipLabels=eventEndOfPeriods$labeling, tipLabelsSimu=eventEndOfPeriods$labeling, comment=comment, n1=eventEndOfPeriods$nLineages1, n2=eventEndOfPeriods$nLineages2)
    }else{
      model <- new(Class="PhenotypicModel", name=keyword, period=eventEndOfPeriods$periods, aAGamma=aAGamma, numbersCopy=eventEndOfPeriods$copy, numbersPaste=eventEndOfPeriods$paste, initialCondition=initialCondition, paramsNames=paramsNames, constraints=constraints, params0=params0, tipLabels=eventEndOfPeriods$labeling, tipLabelsSimu=eventEndOfPeriods$labeling, comment=comment)
    }
  }
  
  else if(keyword == "GMM_all"){
    if (!is.null(geo.object)) { print("ignoring geo.object, using GMM instead")}
    
    comment <- "The Generalist Matching Mutualism 'ALL' model. Assumes equal interaction (S) between all taxa in both trees (inter- and intra-clade)"
    print(comment)
    #comment <- "Generalist Matching Mutualism model.\nStarts with 3 or 4 lineages having the same value X_0 ~ Normal(m0,v0).\nOne trait in each lineage, all lineages evolving then non-independtly according to the GMM expression."
    paramsNames <- c("m0", "v0", "d1", "d2", "S", "sigma")
    params0 <- c(0,0,1,-1,0.5,1)
    
    eventEndOfPeriods <- endOfPeriodsGMM(tree1, tree2)
    n <- eventEndOfPeriods$nLineages1[1] + eventEndOfPeriods$nLineages2[1] - 1
    
    initialCondition <- function(params) return( list(mean=rep(params[1], times=n), var=matrix(rep(params[2], times=n*n), nrow=n ) ) ) 
    
    aAGamma <- function(i, params){
      vectorA <- function(t) return( c( rep(params[3]*params[5], times=eventEndOfPeriods$nLineages1[i]), rep(params[4]*params[5], times=eventEndOfPeriods$nLineages2[i]) ) )
      matrixGamma <- function(t) return(diag(params[6], eventEndOfPeriods$nLineages1[i] + eventEndOfPeriods$nLineages2[i]))

      nLineagesTotal <- eventEndOfPeriods$nLineages1[i] + eventEndOfPeriods$nLineages2[i]
      block <- matrix(1, nLineagesTotal, nLineagesTotal)
      matrixA <- -params[5]*(block/(rowSums(block)-1))
      diag(matrixA) <- params[5]
      
      return(list(a=vectorA, A=matrixA, Gamma=matrixGamma))
    } 
    
    constraints <- function(params) return(params[2]>=0 && params[6]>=0)
    
    if( keyword == "GMM_all" ){
      model <- new(Class="PhenotypicGMM", name=keyword, period=eventEndOfPeriods$periods, aAGamma=aAGamma, numbersCopy=eventEndOfPeriods$copy, numbersPaste=eventEndOfPeriods$paste, initialCondition=initialCondition, paramsNames=paramsNames, constraints=constraints, params0=params0, tipLabels=eventEndOfPeriods$labeling, tipLabelsSimu=eventEndOfPeriods$labeling, comment=comment, n1=eventEndOfPeriods$nLineages1, n2=eventEndOfPeriods$nLineages2)
    }else{
      model <- new(Class="PhenotypicModel", name=keyword, period=eventEndOfPeriods$periods, aAGamma=aAGamma, numbersCopy=eventEndOfPeriods$copy, numbersPaste=eventEndOfPeriods$paste, initialCondition=initialCondition, paramsNames=paramsNames, constraints=constraints, params0=params0, tipLabels=eventEndOfPeriods$labeling, tipLabelsSimu=eventEndOfPeriods$labeling, comment=comment)
    }
    class(model)[1] <- "PhenotypicModel"
  }
  
  else if(keyword == "CoEvo") {
    if (is.null(geo.object)) { stop("this model requires a geo.object")}
    
    comment <- "The CoEvo model. An extension of the GMM model, accounting for interactions only between geographic co-occurring lineages. As with the GMM model, it only estimates interaction (S) between taxa across trees (inter-clade, NOT intra-clade). This model also properly accounts for the number of co-occuring lineages by dividing S (Pk/l) using rowsums (see Manceau et al. pg.559, equation 7)."
    print(comment)
    # comment <- "Generalist Matching Mutualism model.
    # \nStarts with 3 or 4 lineages having the same value X_0 ~ Normal(m0,v0).
    # \nOne trait in each lineage, all lineages evolving then non-independtly according to the GMM expression.
    # This was my second attempt at including geography. It assumes no competition between lineages
    # within the same tree (clade), and only competition between trees.
    # It does properly divide S (Pk/l) by the number of overlapping lineages (rowsums).
    # (see Manceau et al. pg.559, equation 7)."
    paramsNames <- c("m0", "v0", "d1", "d2", "S", "sigma")
    params0 <- c(0,0,1,-1,0.5,1)
    
    eventEndOfPeriods <- endOfPeriodsGMMgeo(tree1, tree2, geo.object)
    n <- eventEndOfPeriods$nLineages1[1] + eventEndOfPeriods$nLineages2[1] - 1
    
    initialCondition <- function(params) return( list(mean=rep(params[1], times=n), var=matrix(rep(params[2], times=n*n), nrow=n ) ) ) 
    
    aAGamma <- function(i, params){
      vectorA <- function(t) return( c( rep(params[3]*params[5], times=eventEndOfPeriods$nLineages1[i]), rep(params[4]*params[5], times=eventEndOfPeriods$nLineages2[i]) ) )
      matrixGamma <- function(t) return(diag(params[6], eventEndOfPeriods$nLineages1[i] + eventEndOfPeriods$nLineages2[i]))
      
      block1 <- as.matrix(diag(0, eventEndOfPeriods$nLineages1[i])) # should this actually be diag(1,...) to account for the interaction of a lineage with its own value?
      block2.int <- geo.object$geography.object[[i]][(1:(eventEndOfPeriods$nLineages1[i])),
                                                     ((eventEndOfPeriods$nLineages1[i]+1):(eventEndOfPeriods$nLineages1[i]+eventEndOfPeriods$nLineages2[i]))]
      block2 <- t(t(block2.int))
      
      block3.int <- as.matrix(geo.object$geography.object[[i]][((eventEndOfPeriods$nLineages1[i]+1):(eventEndOfPeriods$nLineages1[i]+eventEndOfPeriods$nLineages2[i])), 
                                                               (1:(eventEndOfPeriods$nLineages1[i]))])
      if(ncol(block3.int)<=1){
        block3 <- t(block3.int)
      } else{block3 <- block3.int}
      
      block4 <- as.matrix(diag(0, eventEndOfPeriods$nLineages2[i])) # should this actually be diag(1,...) to account for the interaction of a lineage with its own value?
      
      matrixA.int <- rbind(cbind(block1, block2), cbind(block3, block4))
      matrixA.int <- -params[5]*sweep(matrixA.int, 1, rowSums(matrixA.int), "/")
      diag(matrixA.int) <- params[5]
      #diag(matrixA.int) <- 1 # was just using this to test the models
      matrixA.int[is.na(matrixA.int)] <- 0;
      matrixA <- matrixA.int
      
      return(list(a=vectorA, A=matrixA, Gamma=matrixGamma))
    }
    
    constraints <- function(params) return(params[2]>=0 && params[6]>=0)
    
    if( keyword == "CoEvo" ){
      model <- new(Class="PhenotypicGMM", name=keyword, period=eventEndOfPeriods$periods, aAGamma=aAGamma, numbersCopy=eventEndOfPeriods$copy, numbersPaste=eventEndOfPeriods$paste, initialCondition=initialCondition, paramsNames=paramsNames, constraints=constraints, params0=params0, tipLabels=eventEndOfPeriods$labeling, tipLabelsSimu=eventEndOfPeriods$labeling, comment=comment, n1=eventEndOfPeriods$nLineages1, n2=eventEndOfPeriods$nLineages2)
    }else{
      model <- new(Class="PhenotypicModel", name=keyword, period=eventEndOfPeriods$periods, aAGamma=aAGamma, numbersCopy=eventEndOfPeriods$copy, numbersPaste=eventEndOfPeriods$paste, initialCondition=initialCondition, paramsNames=paramsNames, constraints=constraints, params0=params0, tipLabels=eventEndOfPeriods$labeling, tipLabelsSimu=eventEndOfPeriods$labeling, comment=comment)
    }
    class(model)[1] <- "PhenotypicModel" # not necessary for this model
  }
  
  else if(keyword == "CoEvo_all"){
    
    comment <- "The CoEvo 'ALL' model. An extension of the GMM model, accounting for interactions only between geographic co-occurring lineages. It is also an extension of the CoEvo model, estimating interaction (S) between all co-occurring taxa (inter-clade AND intra-clade). This model also properly accounts for the number of co-occuring lineages by dividing S (Pk/l) using rowsums (see Manceau et al. pg.559, equation 7)."
    print(comment)
    #comment <- "Generalist Matching Mutualism model.\nStarts with 3 or 4 lineages having the same value X_0 ~ Normal(m0,v0).\nOne trait in each lineage, all lineages evolving then non-independtly according to the GMM expression."
    paramsNames <- c("m0", "v0", "d1", "d2", "S", "sigma")
    params0 <- c(0,0,1,-1,0.5,1)
    
    eventEndOfPeriods <- endOfPeriodsGMMgeo(tree1, tree2, geo.object)
    n <- eventEndOfPeriods$nLineages1[1] + eventEndOfPeriods$nLineages2[1] - 1
    
    initialCondition <- function(params) return( list(mean=rep(params[1], times=n), var=matrix(rep(params[2], times=n*n), nrow=n ) ) ) 
    
    aAGamma <- function(i, params){
      geo.mat <- geo.object$geography.object[[i]]
      vectorA <- function(t) return( c( rep(params[3]*params[5], times=eventEndOfPeriods$nLineages1[i]), rep(params[4]*params[5], times=eventEndOfPeriods$nLineages2[i]) ) )
      matrixGamma <- function(t) return(diag(params[6], eventEndOfPeriods$nLineages1[i] + eventEndOfPeriods$nLineages2[i]))
      
      matrixA <- -params[5]*(geo.mat/(rowSums(geo.mat)-1))
      diag(matrixA) <- params[5]
      matrixA[is.na(matrixA)] <- 0
      
      return(list(a=vectorA, A=matrixA, Gamma=matrixGamma, P=geo.mat))
    } 
    
    constraints <- function(params) return(params[2]>=0 && params[6]>=0)
    
    if( keyword == "CoEvo_all" ){
      model <- new(Class="PhenotypicGMM", name=keyword, period=eventEndOfPeriods$periods, aAGamma=aAGamma, numbersCopy=eventEndOfPeriods$copy, numbersPaste=eventEndOfPeriods$paste, initialCondition=initialCondition, paramsNames=paramsNames, constraints=constraints, params0=params0, tipLabels=eventEndOfPeriods$labeling, tipLabelsSimu=eventEndOfPeriods$labeling, comment=comment, n1=eventEndOfPeriods$nLineages1, n2=eventEndOfPeriods$nLineages2)
    }else{
      model <- new(Class="PhenotypicModel", name=keyword, period=eventEndOfPeriods$periods, aAGamma=aAGamma, numbersCopy=eventEndOfPeriods$copy, numbersPaste=eventEndOfPeriods$paste, initialCondition=initialCondition, paramsNames=paramsNames, constraints=constraints, params0=params0, tipLabels=eventEndOfPeriods$labeling, tipLabelsSimu=eventEndOfPeriods$labeling, comment=comment)
    }
    class(model)[1] <- "PhenotypicModel"
  }
  
  else if(keyword == "CoEvo_Split") {
    if (is.null(geo.object)) { stop("this model requires a geo.object")}
    
    comment <- "The CoEvo 'SPLIT' model. Again, an extension of the GMM model, accounting for interactions only between geographic co-occurring lineages. It accounts for interactions between all taxa like the CoEvo_all model, but estimates a different interaction parameter for intra-clade (S2) and inter-clade (S1) interactions. This model also properly accounts for the number of co-occuring lineages by dividing S (Pk/l) using rowsums (see Manceau et al. pg.559, equation 7)."
    print(comment)
    # comment <- "Generalist Matching Mutualism model.\nStarts with 3 or 4 lineages having the same value X_0 ~ Normal(m0,v0).
    # \nOne trait in each lineage, all lineages evolving then non-independtly according to the GMM expression.
    # This was my second attempt at including competition within a tree, as well as between trees. 
    # It estimates a single S (S2) value across the pair of trees (competition between trees), as 
    # as well as a second S value (S1) for competition within trees.
    # It also properly divides S (Pk,l) by the number of overlapping lineages (rowsums).
    # (see Manceau et al. pg.559, equation 7)."
    paramsNames <- c("m0", "v0", "d1", "d2", "S1", "sigma", "S2")
    params0 <- c(0,0,1,-1,0.5,1,0.5)
    
    eventEndOfPeriods <- endOfPeriodsGMMgeo(tree1, tree2, geo.object)
    n <- eventEndOfPeriods$nLineages1[1] + eventEndOfPeriods$nLineages2[1] - 1
    
    initialCondition <- function(params) return( list(mean=rep(params[1], times=n), var=matrix(rep(params[2], times=n*n), nrow=n ) ) ) 
    
    # THIS MODEL IS NOT ACTUALLY FIXED, IT SHOULD BE SYMMETRIC!
    aAGamma <- function(i, params){
      Pmat <- geo.object$geography.object[[i]]
      vectorA <- function(t) return( c( rep(params[3]*params[5], times=eventEndOfPeriods$nLineages1[i]), rep(params[4]*params[5], times=eventEndOfPeriods$nLineages2[i]) ) )
      matrixGamma <- function(t) return(diag(params[6], eventEndOfPeriods$nLineages1[i] + eventEndOfPeriods$nLineages2[i]))
      
      block1.int <- Pmat[(1:(eventEndOfPeriods$nLineages1[i])), 
                         1:(eventEndOfPeriods$nLineages1[i])]
      block1 <- (block1.int/(rowSums(block1.int)-1)) * -params[7]
      
      block2.int <- Pmat[(1:(eventEndOfPeriods$nLineages1[i])),
                         ((eventEndOfPeriods$nLineages1[i]+1):(eventEndOfPeriods$nLineages1[i]+eventEndOfPeriods$nLineages2[i]))]
      block2.int <- t(t(block2.int))
      block2 <- (block2.int/rowSums(block2.int)) * -params[5]
      
      block3.int <- as.matrix(Pmat[((eventEndOfPeriods$nLineages1[i]+1):(eventEndOfPeriods$nLineages1[i]+eventEndOfPeriods$nLineages2[i])), 
                                   (1:(eventEndOfPeriods$nLineages1[i]))])
      if(ncol(block3.int)<=1){
        block3.int <- t(block3.int); block3 <- (block3.int/rowSums(block3.int) * -params[5])
      } else{
        block3 <- (block3.int/rowSums(block3.int) * -params[5])}
      
      block4.int <- Pmat[((eventEndOfPeriods$nLineages1[i]+1):(eventEndOfPeriods$nLineages1[i]+eventEndOfPeriods$nLineages2[i])),
                         ((eventEndOfPeriods$nLineages1[i]+1):(eventEndOfPeriods$nLineages1[i]+eventEndOfPeriods$nLineages2[i]))]
      if(is.null(nrow(block4.int))){block4 <- as.matrix(block4.int)
      } else{block4 <- (block4.int/(rowSums(block4.int)-1)) * -params[7]}
      
      matrixA.int <- rbind(cbind(block1, block2), cbind(block3, block4))
      matrixA.int[is.na(matrixA.int)] <- 0;
      diag(matrixA.int) <- 0
      diag(matrixA.int) <- -rowSums(matrixA.int)
      matrixA <- matrixA.int
      
      return(list(a=vectorA, A=matrixA, Gamma=matrixGamma, P=Pmat))
    }
    
    constraints <- function(params) return(params[2]>=0 && params[6]>=0) # && params[5]<=0 && params[7]<=0)
    
    if( keyword == "CoEvo_Split" ){
      model <- new(Class="PhenotypicModel", name=keyword, period=eventEndOfPeriods$periods, aAGamma=aAGamma, numbersCopy=eventEndOfPeriods$copy, numbersPaste=eventEndOfPeriods$paste, initialCondition=initialCondition, paramsNames=paramsNames, constraints=constraints, params0=params0, tipLabels=eventEndOfPeriods$labeling, tipLabelsSimu=eventEndOfPeriods$labeling, comment=comment)
    }else{
      model <- new(Class="PhenotypicGMM", name=keyword, period=eventEndOfPeriods$periods, aAGamma=aAGamma, numbersCopy=eventEndOfPeriods$copy, numbersPaste=eventEndOfPeriods$paste, initialCondition=initialCondition, paramsNames=paramsNames, constraints=constraints, params0=params0, tipLabels=eventEndOfPeriods$labeling, tipLabelsSimu=eventEndOfPeriods$labeling, comment=comment, n1=eventEndOfPeriods$nLineages1, n2=eventEndOfPeriods$nLineages2)
    }
  }
  
  else if(keyword == "CoPM_geo") {
    if (is.null(geo.object)) { stop("this model requires a geo.object")}
    
    comment <- "The CoPM 'GEO' model. This is an extension of the CoPM model, which is a joint estimation of the PM model for two trees. It estimates single interaction (S) and rate (sigma) values for both trees, but S is estimated solely from intra-clade interactions (no interaction between trees). It correctly accounts for interaction only among geographic overlapping lineages, and corrects the interaction estimate for the number of overlapping lineages."
    print(comment)

    paramsNames <- c("m0", "v0", "d1", "d2", "S1", "sigma")
    params0 <- c(0,0,1,-1,0.5,1)
    
    eventEndOfPeriods <- endOfPeriodsGMMgeo(tree1, tree2, geo.object)
    n <- eventEndOfPeriods$nLineages1[1] + eventEndOfPeriods$nLineages2[1] - 1
    
    initialCondition <- function(params) return( list(mean=rep(params[1], times=n), var=matrix(rep(params[2], times=n*n), nrow=n ) ) ) 
    
    aAGamma <- function(i, params){
      Pmat <- geo.object$geography.object[[i]]
      vectorA <- function(t) return( c( rep(params[3]*params[5], times=eventEndOfPeriods$nLineages1[i]), rep(params[4]*params[5], times=eventEndOfPeriods$nLineages2[i]) ) )
      matrixGamma <- function(t) return(diag(params[6], eventEndOfPeriods$nLineages1[i] + eventEndOfPeriods$nLineages2[i]))
      
      block1.int <- Pmat[(1:(eventEndOfPeriods$nLineages1[i])), 
                         1:(eventEndOfPeriods$nLineages1[i])]
      block1 <- (block1.int/(rowSums(block1.int)-1)) * -params[5]
      
      block2 <- matrix(0, eventEndOfPeriods$nLineages1[i], eventEndOfPeriods$nLineages2[i])

      block3 <- matrix(0, eventEndOfPeriods$nLineages2[i], eventEndOfPeriods$nLineages1[i])
      
      block4.int <- Pmat[((eventEndOfPeriods$nLineages1[i]+1):(eventEndOfPeriods$nLineages1[i]+eventEndOfPeriods$nLineages2[i])),
                         ((eventEndOfPeriods$nLineages1[i]+1):(eventEndOfPeriods$nLineages1[i]+eventEndOfPeriods$nLineages2[i]))]
      if(is.null(nrow(block4.int))){block4 <- as.matrix(block4.int)
      } else{block4 <- (block4.int/(rowSums(block4.int)-1)) * -params[5]}
      
      matrixA.int <- rbind(cbind(block1, block2), cbind(block3, block4))
      matrixA.int[is.na(matrixA.int)] <- 0;
      diag(matrixA.int) <- 0
      diag(matrixA.int) <- -rowSums(matrixA.int)
      matrixA <- matrixA.int
      
      return(list(a=vectorA, A=matrixA, Gamma=matrixGamma, P=Pmat))
    }
    
    constraints <- function(params) return(params[2]>=0 && params[6]>=0) # && params[5]<=0 && params[7]<=0)
    
    if( keyword == "CoPM_geo" ){
      model <- new(Class="PhenotypicModel", name=keyword, period=eventEndOfPeriods$periods, aAGamma=aAGamma, numbersCopy=eventEndOfPeriods$copy, numbersPaste=eventEndOfPeriods$paste, initialCondition=initialCondition, paramsNames=paramsNames, constraints=constraints, params0=params0, tipLabels=eventEndOfPeriods$labeling, tipLabelsSimu=eventEndOfPeriods$labeling, comment=comment)
    }else{
      model <- new(Class="PhenotypicGMM", name=keyword, period=eventEndOfPeriods$periods, aAGamma=aAGamma, numbersCopy=eventEndOfPeriods$copy, numbersPaste=eventEndOfPeriods$paste, initialCondition=initialCondition, paramsNames=paramsNames, constraints=constraints, params0=params0, tipLabels=eventEndOfPeriods$labeling, tipLabelsSimu=eventEndOfPeriods$labeling, comment=comment, n1=eventEndOfPeriods$nLineages1, n2=eventEndOfPeriods$nLineages2)
    }
  }
  
  else if(keyword == "CoPM") {

    comment <- "The CoPM model. This is a joint estimation of the PM model for two trees. It estimates single interaction (S) and rate (sigma) values for both trees, but S is estimated solely from intra-clade interactions (no interaction between trees). All lineages in a tree are assumed to interact with ALL other lineages in that tree"
    print(comment)

    paramsNames <- c("m0", "v0", "d1", "d2", "S", "sigma")
    params0 <- c(0,0,1,-1,0.5,1)
    
    eventEndOfPeriods <- endOfPeriodsGMM(tree1, tree2)
    n <- eventEndOfPeriods$nLineages1[1] + eventEndOfPeriods$nLineages2[1] - 1
    
    initialCondition <- function(params) return( list(mean=rep(params[1], times=n), var=matrix(rep(params[2], times=n*n), nrow=n ) ) ) 
    
    aAGamma <- function(i, params){
      vectorA <- function(t) return( c( rep(params[3]*params[5], times=eventEndOfPeriods$nLineages1[i]), rep(params[4]*params[5], times=eventEndOfPeriods$nLineages2[i]) ) )
      matrixGamma <- function(t) return(diag(params[6], eventEndOfPeriods$nLineages1[i] + eventEndOfPeriods$nLineages2[i]))
      
      block1.int <- matrix(1, eventEndOfPeriods$nLineages1[i], eventEndOfPeriods$nLineages1[i])
      block1 <- (block1.int/(rowSums(block1.int)-1)) * -params[5]
      
      block2 <- matrix(0, eventEndOfPeriods$nLineages1[i], eventEndOfPeriods$nLineages2[i])
      
      block3 <- matrix(0, eventEndOfPeriods$nLineages2[i], eventEndOfPeriods$nLineages1[i])
      
      block4.int <- matrix(1, eventEndOfPeriods$nLineages2[i], eventEndOfPeriods$nLineages2[i])
        
      if((nrow(block4.int)==1)){block4 <- block4.int # not necessary to do " * -params[5]" because we'll zero the diagonals later anyway
      } else{block4 <- (block4.int/(rowSums(block4.int)-1)) * -params[5]}
      
      matrixA.int <- rbind(cbind(block1, block2), cbind(block3, block4))
      matrixA.int[is.na(matrixA.int)] <- 0; # don't think this is necessary, leaving it anyway
      diag(matrixA.int) <- 0
      diag(matrixA.int) <- -rowSums(matrixA.int)
      matrixA <- matrixA.int
      
      return(list(a=vectorA, A=matrixA, Gamma=matrixGamma))
    }
    
    constraints <- function(params) return(params[2]>=0 && params[6]>=0) # && params[5]<=0 && params[7]<=0)
    
    if( keyword == "CoPM" ){
      model <- new(Class="PhenotypicModel", name=keyword, period=eventEndOfPeriods$periods, aAGamma=aAGamma, numbersCopy=eventEndOfPeriods$copy, numbersPaste=eventEndOfPeriods$paste, initialCondition=initialCondition, paramsNames=paramsNames, constraints=constraints, params0=params0, tipLabels=eventEndOfPeriods$labeling, tipLabelsSimu=eventEndOfPeriods$labeling, comment=comment)
    }else{
      model <- new(Class="PhenotypicGMM", name=keyword, period=eventEndOfPeriods$periods, aAGamma=aAGamma, numbersCopy=eventEndOfPeriods$copy, numbersPaste=eventEndOfPeriods$paste, initialCondition=initialCondition, paramsNames=paramsNames, constraints=constraints, params0=params0, tipLabels=eventEndOfPeriods$labeling, tipLabelsSimu=eventEndOfPeriods$labeling, comment=comment, n1=eventEndOfPeriods$nLineages1, n2=eventEndOfPeriods$nLineages2)
    }
  }
  
  else if(keyword == "JointPM") {
    
    comment <- "The JointPM model. This is a joint estimation of the PM model for two trees. It differs from the CoPM model by estimating separate interaction values for each clade (tree1 = S1; tree2 = S2). All lineages in a tree are assumed to interact with ALL other lineages in that tree."
    print(comment)

    paramsNames <- c("m0", "v0", "d1", "d2", "S1", "sigma", "S2")
    params0 <- c(0,0,1,-1,0.5,1, 0.5)
    
    eventEndOfPeriods <- endOfPeriodsGMM(tree1, tree2)
    n <- eventEndOfPeriods$nLineages1[1] + eventEndOfPeriods$nLineages2[1] - 1
    
    initialCondition <- function(params) return( list(mean=rep(params[1], times=n), var=matrix(rep(params[2], times=n*n), nrow=n ) ) ) 
    
    aAGamma <- function(i, params){
      vectorA <- function(t) return( c( rep(params[3]*params[5], times=eventEndOfPeriods$nLineages1[i]), rep(params[4]*params[5], times=eventEndOfPeriods$nLineages2[i]) ) )
      matrixGamma <- function(t) return(diag(params[6], eventEndOfPeriods$nLineages1[i] + eventEndOfPeriods$nLineages2[i]))
      
      block1.int <- matrix(1, eventEndOfPeriods$nLineages1[i], eventEndOfPeriods$nLineages1[i])
      block1 <- (block1.int/(rowSums(block1.int)-1)) * -params[5]
      
      block2 <- matrix(0, eventEndOfPeriods$nLineages1[i], eventEndOfPeriods$nLineages2[i])
      
      block3 <- matrix(0, eventEndOfPeriods$nLineages2[i], eventEndOfPeriods$nLineages1[i])
      
      block4.int <- matrix(1, eventEndOfPeriods$nLineages2[i], eventEndOfPeriods$nLineages2[i])
      
      if((nrow(block4.int)==1)){block4 <- block4.int # not necessary to do " * -params[5]" because we'll zero the diagonals later anyway
      } else{block4 <- (block4.int/(rowSums(block4.int)-1)) * -params[7]}
      
      matrixA.int <- rbind(cbind(block1, block2), cbind(block3, block4))
      matrixA.int[is.na(matrixA.int)] <- 0; # don't think this is necessary, leaving it anyway
      diag(matrixA.int) <- 0
      diag(matrixA.int) <- -rowSums(matrixA.int)
      matrixA <- matrixA.int
      
      return(list(a=vectorA, A=matrixA, Gamma=matrixGamma))
    }
    
    constraints <- function(params) return(params[2]>=0 && params[6]>=0) # && params[5]<=0 && params[7]<=0)
    
    if( keyword == "JointPM" ){
      model <- new(Class="PhenotypicModel", name=keyword, period=eventEndOfPeriods$periods, aAGamma=aAGamma, numbersCopy=eventEndOfPeriods$copy, numbersPaste=eventEndOfPeriods$paste, initialCondition=initialCondition, paramsNames=paramsNames, constraints=constraints, params0=params0, tipLabels=eventEndOfPeriods$labeling, tipLabelsSimu=eventEndOfPeriods$labeling, comment=comment)
    }else{
      model <- new(Class="PhenotypicGMM", name=keyword, period=eventEndOfPeriods$periods, aAGamma=aAGamma, numbersCopy=eventEndOfPeriods$copy, numbersPaste=eventEndOfPeriods$paste, initialCondition=initialCondition, paramsNames=paramsNames, constraints=constraints, params0=params0, tipLabels=eventEndOfPeriods$labeling, tipLabelsSimu=eventEndOfPeriods$labeling, comment=comment, n1=eventEndOfPeriods$nLineages1, n2=eventEndOfPeriods$nLineages2)
    }
  }
  
  else if(keyword == "JointPM_geo") {
    if (is.null(geo.object)) { stop("this model requires a geo.object")}
    
    comment <- "The JointPM 'GEO' model. This is a joint estimation of the PM model for two trees. It differs from the CoPM_geo model by estimating separate interaction values for each clade (tree1 = S1; tree2 = S2). Like the CoPM_geo (unlike JointPM) it correctly estimates the interaction parameters (S1,S2) for only geographic overlapping taxa (it also corrects for the number of taxa overlapping)."
    print(comment)

    paramsNames <- c("m0", "v0", "d1", "d2", "S1", "sigma", "S2")
    params0 <- c(0,0,1,-1,0.5,1,0.5)
    
    eventEndOfPeriods <- endOfPeriodsGMMgeo(tree1, tree2, geo.object)
    n <- eventEndOfPeriods$nLineages1[1] + eventEndOfPeriods$nLineages2[1] - 1
    
    initialCondition <- function(params) return( list(mean=rep(params[1], times=n), var=matrix(rep(params[2], times=n*n), nrow=n ) ) ) 
    
    aAGamma <- function(i, params){
      Pmat <- geo.object$geography.object[[i]]
      vectorA <- function(t) return( c( rep(params[3]*params[5], times=eventEndOfPeriods$nLineages1[i]), rep(params[4]*params[5], times=eventEndOfPeriods$nLineages2[i]) ) )
      matrixGamma <- function(t) return(diag(params[6], eventEndOfPeriods$nLineages1[i] + eventEndOfPeriods$nLineages2[i]))
      
      block1.int <- Pmat[(1:(eventEndOfPeriods$nLineages1[i])), 
                         1:(eventEndOfPeriods$nLineages1[i])]
      block1 <- (block1.int/(rowSums(block1.int)-1)) * -params[5]
      
      block2 <- matrix(0, eventEndOfPeriods$nLineages1[i], eventEndOfPeriods$nLineages2[i])
      
      block3 <- matrix(0, eventEndOfPeriods$nLineages2[i], eventEndOfPeriods$nLineages1[i])
      
      block4.int <- Pmat[((eventEndOfPeriods$nLineages1[i]+1):(eventEndOfPeriods$nLineages1[i]+eventEndOfPeriods$nLineages2[i])),
                         ((eventEndOfPeriods$nLineages1[i]+1):(eventEndOfPeriods$nLineages1[i]+eventEndOfPeriods$nLineages2[i]))]
      if(is.null(nrow(block4.int))){block4 <- as.matrix(block4.int)
      } else{block4 <- (block4.int/(rowSums(block4.int)-1)) * -params[7]}
      
      matrixA.int <- rbind(cbind(block1, block2), cbind(block3, block4))
      matrixA.int[is.na(matrixA.int)] <- 0;
      diag(matrixA.int) <- 0
      diag(matrixA.int) <- -rowSums(matrixA.int)
      matrixA <- matrixA.int
      
      return(list(a=vectorA, A=matrixA, Gamma=matrixGamma, P=Pmat))
    }
    
    constraints <- function(params) return(params[2]>=0 && params[6]>=0) # && params[5]<=0 && params[7]<=0)
    
    if( keyword == "JointPM_geo" ){
      model <- new(Class="PhenotypicModel", name=keyword, period=eventEndOfPeriods$periods, aAGamma=aAGamma, numbersCopy=eventEndOfPeriods$copy, numbersPaste=eventEndOfPeriods$paste, initialCondition=initialCondition, paramsNames=paramsNames, constraints=constraints, params0=params0, tipLabels=eventEndOfPeriods$labeling, tipLabelsSimu=eventEndOfPeriods$labeling, comment=comment)
    }else{
      model <- new(Class="PhenotypicGMM", name=keyword, period=eventEndOfPeriods$periods, aAGamma=aAGamma, numbersCopy=eventEndOfPeriods$copy, numbersPaste=eventEndOfPeriods$paste, initialCondition=initialCondition, paramsNames=paramsNames, constraints=constraints, params0=params0, tipLabels=eventEndOfPeriods$labeling, tipLabelsSimu=eventEndOfPeriods$labeling, comment=comment, n1=eventEndOfPeriods$nLineages1, n2=eventEndOfPeriods$nLineages2)
    }
  }
  
  else if(keyword == "CoBM"){
    
    comment <- "Brownian Motion model with linear drift.\nStarts with two lineages having the same value X_0 ~ Normal(m0,v0).\nOne trait in each lineage, all lineages evolving independently after branching.\ndX_t = d dt + sigma dW_t"
    paramsNames <- c("m0", "v0", "d", "sigma")
    params0 <- c(0,0,0,1)
    
    #periodizing <- periodizeOneTree(tree)
    #eventEndOfPeriods <- endOfPeriods(periodizing, tree)
    
    eventEndOfPeriods <- endOfPeriodsGMM(tree1, tree2)
    n <- eventEndOfPeriods$nLineages1[1] + eventEndOfPeriods$nLineages2[1] - 1
    
    #initialCondition <- function(params) return( list(mean=c(params[1]), var=matrix(c(params[2]))) )
    initialCondition <- function(params) return( list(mean=rep(params[1], times=n), var=matrix(rep(params[2], times=n*n), nrow=n ) ) ) 
    
    
    aAGamma <- function(i, params){
      #vectorU <- getLivingLineages(i, eventEndOfPeriods)
      vectorU <- rep(1, (eventEndOfPeriods$nLineages1[i] + eventEndOfPeriods$nLineages2[i]))
      vectorA <- function(t) return(params[3]*vectorU)
      matrixGamma <- function(t) return(params[4]*diag(vectorU))
      matrixA <- diag(0, length(vectorU))
      
      return(list(a=vectorA, A=matrixA, Gamma=matrixGamma))
    }
    
    constraints <- function(params) return(params[2]>=0 && params[4]>=0)
    
    if( keyword == "CoBM" ){
      model <- new(Class="PhenotypicBM", name=keyword, period=eventEndOfPeriods$periods, aAGamma=aAGamma, numbersCopy=eventEndOfPeriods$copy, numbersPaste=eventEndOfPeriods$paste, initialCondition=initialCondition, paramsNames=paramsNames, constraints=constraints, params0=params0, tipLabels=eventEndOfPeriods$labeling, tipLabelsSimu=eventEndOfPeriods$labeling, comment=comment)
    }else{
      model <- new(Class="PhenotypicModel", name=keyword, period=eventEndOfPeriods$periods, aAGamma=aAGamma, numbersCopy=eventEndOfPeriods$copy, numbersPaste=eventEndOfPeriods$paste, initialCondition=initialCondition, paramsNames=paramsNames, constraints=constraints, params0=params0, tipLabels=eventEndOfPeriods$labeling, tipLabelsSimu=eventEndOfPeriods$labeling, comment=comment)
    }
    class(model)[1] <- "PhenotypicModel"
  }
  
  else if(keyword == "CoOU"){
    
    comment <- "Ornstein-Uhlenbeck model.\nStarts with two lineages having the same value X_0 ~ Normal(m0,v0).\nOne trait in each lineage, all lineages evolving independently after branching.\ndX_t = psi(theta- X_t) dt + sigma dW_t"
    paramsNames <- c("m0", "v0", "psi", "theta", "sigma")
    params0 <- c(0,0,1,0,1)
    
    #periodizing <- periodizeOneTree(tree)
    #eventEndOfPeriods <- endOfPeriods(periodizing, tree)
    
    eventEndOfPeriods <- endOfPeriodsGMM(tree1, tree2)
    n <- eventEndOfPeriods$nLineages1[1] + eventEndOfPeriods$nLineages2[1] - 1
    
    #initialCondition <- function(params) return( list(mean=c(params[1]), var=matrix(c(params[2]))) )
    initialCondition <- function(params) return( list(mean=rep(params[1], times=n), var=matrix(rep(params[2], times=n*n), nrow=n ) ) ) 
    
    aAGamma <- function(i, params){
      #vectorU <- getLivingLineages(i, eventEndOfPeriods)
      vectorU <- rep(1, (eventEndOfPeriods$nLineages1[i] + eventEndOfPeriods$nLineages2[i]))
      vectorA <- function(t) return(params[3]*params[4]*vectorU)
      matrixGamma <- function(t) return(params[5]*diag(vectorU))
      matrixA <- params[3]*diag(vectorU)
      
      return(list(a=vectorA, A=matrixA, Gamma=matrixGamma))
    }
    
    #constraints <- function(params) return(params[2]>=0 && params[5]>=0 && params[3]!=0)
    constraints <- function(params) return(params[2]>=0 && params[5]>=0 && params[3]>0)
    
    
    if( keyword == "CoOU" ){
      model <- new(Class="PhenotypicOU", name=keyword, period=eventEndOfPeriods$periods, aAGamma=aAGamma, numbersCopy=eventEndOfPeriods$copy, numbersPaste=eventEndOfPeriods$paste, initialCondition=initialCondition, paramsNames=paramsNames, constraints=constraints, params0=params0, tipLabels=eventEndOfPeriods$labeling, tipLabelsSimu=eventEndOfPeriods$labeling, comment=comment)
    }else{
      model <- new(Class="PhenotypicModel", name=keyword, period=eventEndOfPeriods$periods, aAGamma=aAGamma, numbersCopy=eventEndOfPeriods$copy, numbersPaste=eventEndOfPeriods$paste, initialCondition=initialCondition, paramsNames=paramsNames, constraints=constraints, params0=params0, tipLabels=eventEndOfPeriods$labeling, tipLabelsSimu=eventEndOfPeriods$labeling, comment=comment)
    }
    class(model)[1] <- "PhenotypicModel"
  }
  
  
  # this is the list of models I've developed but probably shouldn't be fit for one reason or another (I've added "_" to avoid using them)
  else if(keyword == "_GMM_driftless"){
    if (!is.null(geo.object)) { print("ignoring geo.object, using GMM instead")}
    
    print("The GMM model, assumes equal interaction (S) between all lineages of both trees, but no interaction among lineages within a tree. This version does not estimate the drift parameters (d1,d2) which are hard to fit on small trees, or those with few interactions")
    comment <- "Generalist Matching Mutualism model.\nStarts with 3 or 4 lineages having the same value X_0 ~ Normal(m0,v0).\nOne trait in each lineage, all lineages evolving then non-independtly according to the GMM expression."
    #paramsNames <- c("m0", "v0", "d1", "d2", "S", "sigma")
    #params0 <- c(0,0,1,-1,0.5,1)
    
    paramsNames <- c("m0", "v0", "S", "sigma")
    params0 <- c(0,0,0.5,1)
    
    eventEndOfPeriods <- endOfPeriodsGMM(tree1, tree2)
    n <- eventEndOfPeriods$nLineages1[1] + eventEndOfPeriods$nLineages2[1] - 1
    
    initialCondition <- function(params) return( list(mean=rep(params[1], times=n), var=matrix(rep(params[2], times=n*n), nrow=n ) ) ) 
    
    aAGamma <- function(i, params){
      vectorA <- function(t) return( c( rep(params[3], times=eventEndOfPeriods$nLineages1[i]), rep(params[3], times=eventEndOfPeriods$nLineages2[i]) ) )
      matrixGamma <- function(t) return(diag(params[4], eventEndOfPeriods$nLineages1[i] + eventEndOfPeriods$nLineages2[i]))
      
      bloc1 <- diag(params[3], eventEndOfPeriods$nLineages1[i])
      bloc2 <- matrix(rep(-params[3]/eventEndOfPeriods$nLineages2[i], times=eventEndOfPeriods$nLineages1[i]*eventEndOfPeriods$nLineages2[i]), nrow=eventEndOfPeriods$nLineages1[i])
      bloc3 <- matrix(rep(-params[3]/eventEndOfPeriods$nLineages1[i], times=eventEndOfPeriods$nLineages1[i]*eventEndOfPeriods$nLineages2[i]), nrow=eventEndOfPeriods$nLineages2[i])
      bloc4 <- diag(params[3], eventEndOfPeriods$nLineages2[i])
      matrixA <- rbind(cbind(bloc1, bloc2), cbind(bloc3, bloc4))
      
      return(list(a=vectorA, A=matrixA, Gamma=matrixGamma))
    } 
    
    constraints <- function(params) return(params[2]>=0 && params[4]>=0)
    
    if( keyword == "GMM_driftless" ){
      model <- new(Class="PhenotypicGMM_driftless", name=keyword, period=eventEndOfPeriods$periods, aAGamma=aAGamma, numbersCopy=eventEndOfPeriods$copy, numbersPaste=eventEndOfPeriods$paste, initialCondition=initialCondition, paramsNames=paramsNames, constraints=constraints, params0=params0, tipLabels=eventEndOfPeriods$labeling, tipLabelsSimu=eventEndOfPeriods$labeling, comment=comment, n1=eventEndOfPeriods$nLineages1, n2=eventEndOfPeriods$nLineages2)
    }else{
      model <- new(Class="PhenotypicModel", name=keyword, period=eventEndOfPeriods$periods, aAGamma=aAGamma, numbersCopy=eventEndOfPeriods$copy, numbersPaste=eventEndOfPeriods$paste, initialCondition=initialCondition, paramsNames=paramsNames, constraints=constraints, params0=params0, tipLabels=eventEndOfPeriods$labeling, tipLabelsSimu=eventEndOfPeriods$labeling, comment=comment)
    }
  }
  
  else if(keyword == "_GMM+geo") {
    if (is.null(geo.object)) { stop("this model requires a geo.object")}
    
    print("GMM+geo: A GMM model incorporating geography, but assumes equal interaction with all overlapping lineages, and estimates S between trees only. Probably incorrect implementation of geography, see notes in code 'RPANDA_extras' for specifics")
    comment <- "Generalist Matching Mutualism model.
    \nStarts with 3 or 4 lineages having the same value X_0 ~ Normal(m0,v0).
    \nOne trait in each lineage, all lineages evolving then non-independtly according to the GMM expression.
    This was my original attempt at including geography, but all it does is multiply the A matrix by the geo.object.
    It works in a basic way, but it assumes that if a species overlaps with 1 lineage or 20 lineages, 
    that they all get the same S value, instead of dividing the S value by the number of overlapping
    lineages (see Manceau et al. pg.559, equation 7)."
    paramsNames <- c("m0", "v0", "d1", "d2", "S", "sigma")
    params0 <- c(0,0,1,-1,0.5,1)

    eventEndOfPeriods <- endOfPeriodsGMMgeo(tree1, tree2, geo.object)
    n <- eventEndOfPeriods$nLineages1[1] + eventEndOfPeriods$nLineages2[1] - 1
    
    initialCondition <- function(params) return( list(mean=rep(params[1], times=n), var=matrix(rep(params[2], times=n*n), nrow=n ) ) ) 
    
    aAGamma <- function(i, params){
      vectorA <- function(t) return( c( rep(params[3]*params[5], times=eventEndOfPeriods$nLineages1[i]), rep(params[4]*params[5], times=eventEndOfPeriods$nLineages2[i]) ) )
      matrixGamma <- function(t) return(diag(params[6], eventEndOfPeriods$nLineages1[i] + eventEndOfPeriods$nLineages2[i]))
      
      bloc1 <- diag(params[5], eventEndOfPeriods$nLineages1[i])
      bloc2 <- matrix(rep(-params[5]/eventEndOfPeriods$nLineages2[i], times=eventEndOfPeriods$nLineages1[i]*eventEndOfPeriods$nLineages2[i]), nrow=eventEndOfPeriods$nLineages1[i])
      bloc3 <- matrix(rep(-params[5]/eventEndOfPeriods$nLineages1[i], times=eventEndOfPeriods$nLineages1[i]*eventEndOfPeriods$nLineages2[i]), nrow=eventEndOfPeriods$nLineages2[i])
      bloc4 <- diag(params[5], eventEndOfPeriods$nLineages2[i])
      matrixA.int <- rbind(cbind(bloc1, bloc2), cbind(bloc3, bloc4))
      # I think all I'd need to do is:
      matrixA <- matrixA.int * geo.object$geography.object[[i]]
      #diag(matrixA) <- 0
      
      
      return(list(a=vectorA, A=matrixA, Gamma=matrixGamma))
    }
    
    constraints <- function(params) return(params[2]>=0 && params[6]>=0)
    
    if( keyword == "GMM+geo" ){
      model <- new(Class="PhenotypicGMM", name=keyword, period=eventEndOfPeriods$periods, aAGamma=aAGamma, numbersCopy=eventEndOfPeriods$copy, numbersPaste=eventEndOfPeriods$paste, initialCondition=initialCondition, paramsNames=paramsNames, constraints=constraints, params0=params0, tipLabels=eventEndOfPeriods$labeling, tipLabelsSimu=eventEndOfPeriods$labeling, comment=comment, n1=eventEndOfPeriods$nLineages1, n2=eventEndOfPeriods$nLineages2)
    }else{
      model <- new(Class="PhenotypicModel", name=keyword, period=eventEndOfPeriods$periods, aAGamma=aAGamma, numbersCopy=eventEndOfPeriods$copy, numbersPaste=eventEndOfPeriods$paste, initialCondition=initialCondition, paramsNames=paramsNames, constraints=constraints, params0=params0, tipLabels=eventEndOfPeriods$labeling, tipLabelsSimu=eventEndOfPeriods$labeling, comment=comment)
    }
    #class(model)[1] <- "PhenotypicModel"
  }
  
  else if (keyword == "_JointMC") {
    
    print("JointMC: this doesn't work yet, but might be interesting to test if a single S value fits both trees")
    comment = "this is experimental"
    paramsNames <- c("m0","logsigma","S")
    params0 <- c(0,log(1),0)
    
    eventEndOfPeriods <- endOfPeriodsGMM(tree1, tree2)
    n <- eventEndOfPeriods$nLineages1[1] + eventEndOfPeriods$nLineages2[1] - 1
    
    initialCondition <- function(params) return( list(mean=c(params[1]), var=matrix(c(0))) ) 
    
    ###is this where the A matrix incorporating geography needs to go? if so, what is the order in which lineage sympatry data need to be introduced

    aAGamma <- function(i, params){
      vectorU <- rep(1, sum(eventEndOfPeriods$nLineages1[i] + eventEndOfPeriods$nLineages2[i]))
      vectorA <- function(t) return(0*vectorU)
      matrixGamma <- function(t) return(exp(params[2])*diag(vectorU))
      # nij <- colSums(geo.object$geography.object[[i]])
      # matrixA <- params[3]*diag(vectorU) - ((params[3]/nij))
      matrixA <- params[3]*diag(vectorU) - (params[3]/sum(vectorU)) * outer(vectorU,vectorU)
      return(list(a=vectorA, A=matrixA, Gamma=matrixGamma))
    }
    constraints <- function(params) return(params[3]<=0)
    #model <- new(Class="PhenotypicGMM", name=keyword, period=eventEndOfPeriods$periods, aAGamma=aAGamma, numbersCopy=eventEndOfPeriods$copy, numbersPaste=eventEndOfPeriods$paste, initialCondition=initialCondition, paramsNames=paramsNames, constraints=constraints, params0=params0, tipLabels=eventEndOfPeriods$labeling,                                           comment=comment, n1=eventEndOfPeriods$nLineages1, n2=eventEndOfPeriods$nLineages2)
    model <- new(Class="PhenotypicGMM", name=keyword, period=eventEndOfPeriods$periods, aAGamma=aAGamma, numbersCopy=eventEndOfPeriods$copy, numbersPaste=eventEndOfPeriods$paste, initialCondition=initialCondition, paramsNames=paramsNames, constraints=constraints, params0=params0, tipLabels=eventEndOfPeriods$labeling, tipLabelsSimu=eventEndOfPeriods$labeling, comment=comment, n1=eventEndOfPeriods$nLineages1, n2=eventEndOfPeriods$nLineages2)
    return(model)
  }
  
  else if(keyword == "_CoEvo") {
    if (is.null(geo.object)) { stop("this model requires a geo.object")}
    print("CoEvo: A GMM model incorporating geography (correcting for the number of interacting lineages), but estimating S between trees only (no within-tree interaction).")
    comment <- "Generalist Matching Mutualism model.
    \nStarts with 3 or 4 lineages having the same value X_0 ~ Normal(m0,v0).
    \nOne trait in each lineage, all lineages evolving then non-independtly according to the GMM expression.
    This was my second attempt at including geography. It assumes no competition between lineages
    within the same tree (clade), and only competition between trees.
    It does properly divide S (Pk/l) by the number of overlapping lineages (rowsums).
    (see Manceau et al. pg.559, equation 7)."
    paramsNames <- c("m0", "v0", "d1", "d2", "S", "sigma")
    params0 <- c(0,0,1,-1,0.5,1)

    eventEndOfPeriods <- endOfPeriodsGMMgeo(tree1, tree2, geo.object)
    n <- eventEndOfPeriods$nLineages1[1] + eventEndOfPeriods$nLineages2[1] - 1
    
    initialCondition <- function(params) return( list(mean=rep(params[1], times=n), var=matrix(rep(params[2], times=n*n), nrow=n ) ) ) 
    
    aAGamma <- function(i, params){
      vectorA <- function(t) return( c( rep(params[3]*params[5], times=eventEndOfPeriods$nLineages1[i]), rep(params[4]*params[5], times=eventEndOfPeriods$nLineages2[i]) ) )
      matrixGamma <- function(t) return(diag(params[6], eventEndOfPeriods$nLineages1[i] + eventEndOfPeriods$nLineages2[i]))
      
      bloc1 <- diag(params[5], eventEndOfPeriods$nLineages1[i])
      # bloc1 <- params[5]*(geo.object$geography.object[[i]][(1:(eventEndOfPeriods$nLineages1[i])),
      #                                                      (1:(eventEndOfPeriods$nLineages1[i]))])
      
      #bloc2 <- matrix(rep(-params[5]/eventEndOfPeriods$nLineages2[i], times=eventEndOfPeriods$nLineages1[i]*eventEndOfPeriods$nLineages2[i]), nrow=eventEndOfPeriods$nLineages1[i])
      bloc2.int <- geo.object$geography.object[[i]][(1:(eventEndOfPeriods$nLineages1[i])),
                                                    ((eventEndOfPeriods$nLineages1[i]+1):(eventEndOfPeriods$nLineages1[i]+eventEndOfPeriods$nLineages2[i]))]
      if(is.null(nrow(bloc2.int))){
        bloc2 <- as.matrix(params[5]*(bloc2.int/bloc2.int)); bloc2 <- as.matrix(bloc2)
      } else {
        bloc2 <- params[5]*(sweep(bloc2.int, 2, (colSums(bloc2.int)-1), "/"))
      }
      bloc2[is.na(bloc2)] <- 0; bloc2[is.infinite(bloc2)] <- 0
      
      #bloc3 <- matrix(rep(-params[5]/eventEndOfPeriods$nLineages1[i], times=eventEndOfPeriods$nLineages1[i]*eventEndOfPeriods$nLineages2[i]), nrow=eventEndOfPeriods$nLineages2[i])
      bloc3.int <- geo.object$geography.object[[i]][((eventEndOfPeriods$nLineages1[i]+1):(eventEndOfPeriods$nLineages1[i]+eventEndOfPeriods$nLineages2[i])), 
                                                    (1:(eventEndOfPeriods$nLineages1[i]))]
      if(is.null(nrow(bloc3.int))){
        bloc3 <- t(params[5]*(bloc3.int/bloc3.int)); bloc3 <- as.matrix(bloc3)
      } else {
        bloc3 <- params[5]*(sweep(bloc3.int, 1, (rowSums(bloc3.int)-1), "/"))
      }
      bloc3[is.na(bloc3)] <- 0; bloc3[is.infinite(bloc3)] <- 0
      
      bloc4 <- diag(params[5], eventEndOfPeriods$nLineages2[i])
      matrixA.int <- rbind(cbind(bloc1, bloc2), cbind(bloc3, bloc4))
      diag(matrixA.int) <- 0; diag(matrixA.int) <- abs(rowSums(matrixA.int))
      matrixA <- matrixA.int

      return(list(a=vectorA, A=matrixA, Gamma=matrixGamma))
    }
    
    constraints <- function(params) return(params[2]>=0 && params[6]>=0)
    
    if( keyword == "CoEvo" ){
      model <- new(Class="PhenotypicGMM", name=keyword, period=eventEndOfPeriods$periods, aAGamma=aAGamma, numbersCopy=eventEndOfPeriods$copy, numbersPaste=eventEndOfPeriods$paste, initialCondition=initialCondition, paramsNames=paramsNames, constraints=constraints, params0=params0, tipLabels=eventEndOfPeriods$labeling, tipLabelsSimu=eventEndOfPeriods$labeling, comment=comment, n1=eventEndOfPeriods$nLineages1, n2=eventEndOfPeriods$nLineages2)
    }else{
      model <- new(Class="PhenotypicModel", name=keyword, period=eventEndOfPeriods$periods, aAGamma=aAGamma, numbersCopy=eventEndOfPeriods$copy, numbersPaste=eventEndOfPeriods$paste, initialCondition=initialCondition, paramsNames=paramsNames, constraints=constraints, params0=params0, tipLabels=eventEndOfPeriods$labeling, tipLabelsSimu=eventEndOfPeriods$labeling, comment=comment)
    }
    #class(model)[1] <- "PhenotypicModel" # not necessary for this model
  }
  
  else if(keyword == "_oldCoEvo") {
    if (is.null(geo.object)) { stop("this model requires a geo.object")}
    print("CoEvo: A GMM model incorporating geography (correcting for the number of interacting lineages), but estimating S between trees only (no within-tree interaction)")
    comment <- "Generalist Matching Mutualism model.
    \nStarts with 3 or 4 lineages having the same value X_0 ~ Normal(m0,v0).
    \nOne trait in each lineage, all lineages evolving then non-independtly according to the GMM expression.
    This was my second attempt at including geography. It assumes no competition between lineages
    within the same tree (clade), and only competition between trees.
    It does properly divide S (Pk/l) by the number of overlapping lineages (rowsums).
    (see Manceau et al. pg.559, equation 7)."
    paramsNames <- c("m0", "v0", "d1", "d2", "S", "sigma")
    params0 <- c(0,0,1,-1,0.5,1)
    
    eventEndOfPeriods <- endOfPeriodsGMMgeo(tree1, tree2, geo.object)
    n <- eventEndOfPeriods$nLineages1[1] + eventEndOfPeriods$nLineages2[1] - 1
    
    initialCondition <- function(params) return( list(mean=rep(params[1], times=n), var=matrix(rep(params[2], times=n*n), nrow=n ) ) ) 
    
    aAGamma <- function(i, params){
      vectorA <- function(t) return( c( rep(params[3]*params[5], times=eventEndOfPeriods$nLineages1[i]), rep(params[4]*params[5], times=eventEndOfPeriods$nLineages2[i]) ) )
      matrixGamma <- function(t) return(diag(params[6], eventEndOfPeriods$nLineages1[i] + eventEndOfPeriods$nLineages2[i]))
      
      bloc1 <- diag(params[5], eventEndOfPeriods$nLineages1[i])
      # bloc1 <- params[5]*(geo.object$geography.object[[i]][(1:(eventEndOfPeriods$nLineages1[i])),
      #                                                      (1:(eventEndOfPeriods$nLineages1[i]))])
      
      #bloc2 <- matrix(rep(-params[5]/eventEndOfPeriods$nLineages2[i], times=eventEndOfPeriods$nLineages1[i]*eventEndOfPeriods$nLineages2[i]), nrow=eventEndOfPeriods$nLineages1[i])
      bloc2.int <- geo.object$geography.object[[i]][(1:(eventEndOfPeriods$nLineages1[i])),
                                                    ((eventEndOfPeriods$nLineages1[i]+1):(eventEndOfPeriods$nLineages1[i]+eventEndOfPeriods$nLineages2[i]))]
      if(is.null(nrow(bloc2.int))){
        bloc2 <- as.matrix(t(-params[5]*(bloc2.int/bloc2.int)))
      } else {
        bloc2 <- -params[5]*(sweep(bloc2.int, 2, (colSums(bloc2.int)-1), "/"))
      }
      bloc2[is.na(bloc2)] <- 0; bloc2[is.infinite(bloc2)] <- 0
      
      #bloc3 <- matrix(rep(-params[5]/eventEndOfPeriods$nLineages1[i], times=eventEndOfPeriods$nLineages1[i]*eventEndOfPeriods$nLineages2[i]), nrow=eventEndOfPeriods$nLineages2[i])
      bloc3.int <- geo.object$geography.object[[i]][((eventEndOfPeriods$nLineages1[i]+1):(eventEndOfPeriods$nLineages1[i]+eventEndOfPeriods$nLineages2[i])), 
                                                    (1:(eventEndOfPeriods$nLineages1[i]))]
      if(is.null(nrow(bloc3.int))){
        bloc3 <- t(t(-params[5]*(bloc3.int/bloc3.int)))
      } else {
        bloc3 <- -params[5]*(sweep(bloc3.int, 1, (rowSums(bloc3.int)-1), "/"))
      }
      bloc3[is.na(bloc3)] <- 0; bloc3[is.infinite(bloc3)] <- 0
      
      bloc4 <- diag(params[5], eventEndOfPeriods$nLineages2[i])
      matrixA <- rbind(cbind(bloc1, bloc2), cbind(bloc3, bloc4))
      
      
      return(list(a=vectorA, A=matrixA, Gamma=matrixGamma))
    }
    
    constraints <- function(params) return(params[2]>=0 && params[6]>=0)
    
    if( keyword == "oldCoEvo" ){
      model <- new(Class="PhenotypicGMM", name=keyword, period=eventEndOfPeriods$periods, aAGamma=aAGamma, numbersCopy=eventEndOfPeriods$copy, numbersPaste=eventEndOfPeriods$paste, initialCondition=initialCondition, paramsNames=paramsNames, constraints=constraints, params0=params0, tipLabels=eventEndOfPeriods$labeling, tipLabelsSimu=eventEndOfPeriods$labeling, comment=comment, n1=eventEndOfPeriods$nLineages1, n2=eventEndOfPeriods$nLineages2)
    }else{
      model <- new(Class="PhenotypicModel", name=keyword, period=eventEndOfPeriods$periods, aAGamma=aAGamma, numbersCopy=eventEndOfPeriods$copy, numbersPaste=eventEndOfPeriods$paste, initialCondition=initialCondition, paramsNames=paramsNames, constraints=constraints, params0=params0, tipLabels=eventEndOfPeriods$labeling, tipLabelsSimu=eventEndOfPeriods$labeling, comment=comment)
    }
    class(model)[1] <- "PhenotypicModel"
  }
  
  else if(keyword == "_CoEvoII") {
    if (is.null(geo.object)) { stop("this model requires a geo.object")}
    print("CoEvo: A GMM model incorporating geography (correcting for the number of interacting lineages), but estimating S between trees only (no within-tree interaction).")
    comment <- "Generalist Matching Mutualism model.
    \nStarts with 3 or 4 lineages having the same value X_0 ~ Normal(m0,v0).
    \nOne trait in each lineage, all lineages evolving then non-independtly according to the GMM expression.
    This was my second attempt at including geography. It assumes no competition between lineages
    within the same tree (clade), and only competition between trees.
    It does properly divide S (Pk/l) by the number of overlapping lineages (rowsums).
    (see Manceau et al. pg.559, equation 7)."
    paramsNames <- c("m0", "v0", "d1", "d2", "S", "sigma")
    params0 <- c(0,0,1,-1,0.5,1)
    
    eventEndOfPeriods <- endOfPeriodsGMMgeo(tree1, tree2, geo.object)
    n <- eventEndOfPeriods$nLineages1[1] + eventEndOfPeriods$nLineages2[1] - 1
    
    initialCondition <- function(params) return( list(mean=rep(params[1], times=n), var=matrix(rep(params[2], times=n*n), nrow=n ) ) ) 
    
    aAGamma <- function(i, params){
      vectorA <- function(t) return( c( rep(params[3]*params[5], times=eventEndOfPeriods$nLineages1[i]), rep(params[4]*params[5], times=eventEndOfPeriods$nLineages2[i]) ) )
      matrixGamma <- function(t) return(diag(params[6], eventEndOfPeriods$nLineages1[i] + eventEndOfPeriods$nLineages2[i]))
      
      block1 <- diag(params[5], eventEndOfPeriods$nLineages1[i])
      
      block2.int <- geo.object$geography.object[[i]][(1:(eventEndOfPeriods$nLineages1[i])),
                                                     ((eventEndOfPeriods$nLineages1[i]+1):(eventEndOfPeriods$nLineages1[i]+eventEndOfPeriods$nLineages2[i]))]
      block2 <- t(t(block2.int))
      
      block3.int <- as.matrix(geo.object$geography.object[[i]][((eventEndOfPeriods$nLineages1[i]+1):(eventEndOfPeriods$nLineages1[i]+eventEndOfPeriods$nLineages2[i])), 
                                                               (1:(eventEndOfPeriods$nLineages1[i]))])
      if(ncol(block3.int)<=1){
        block3 <- t(block3.int)
      } else{block3 <- block3.int}
      
      block4 <- diag(params[5], eventEndOfPeriods$nLineages2[i])
      
      matrixA.int <- rbind(cbind(block1, block2), cbind(block3, block4))
      matrixA.int <- -params[5]*sweep(matrixA.int, 1, rowSums(matrixA.int)-params[5], "/")
      diag(matrixA.int) <- params[5]
      matrixA <- matrixA.int
      
      return(list(a=vectorA, A=matrixA, Gamma=matrixGamma))
    }
    
    constraints <- function(params) return(params[2]>=0 && params[6]>=0)
    
    if( keyword == "CoEvoII" ){
      model <- new(Class="PhenotypicGMM", name=keyword, period=eventEndOfPeriods$periods, aAGamma=aAGamma, numbersCopy=eventEndOfPeriods$copy, numbersPaste=eventEndOfPeriods$paste, initialCondition=initialCondition, paramsNames=paramsNames, constraints=constraints, params0=params0, tipLabels=eventEndOfPeriods$labeling, tipLabelsSimu=eventEndOfPeriods$labeling, comment=comment, n1=eventEndOfPeriods$nLineages1, n2=eventEndOfPeriods$nLineages2)
    }else{
      model <- new(Class="PhenotypicModel", name=keyword, period=eventEndOfPeriods$periods, aAGamma=aAGamma, numbersCopy=eventEndOfPeriods$copy, numbersPaste=eventEndOfPeriods$paste, initialCondition=initialCondition, paramsNames=paramsNames, constraints=constraints, params0=params0, tipLabels=eventEndOfPeriods$labeling, tipLabelsSimu=eventEndOfPeriods$labeling, comment=comment)
    }
    class(model)[1] <- "PhenotypicModel" # not necessary for this model
  }
  
  else if(keyword == "_CoEvo_driftless") {
    if (is.null(geo.object)) { stop("this model requires a geo.object")}
    print("CoEvo: A GMM model incorporating geography (correcting for the number of interacting lineages), but estimating S between trees only (no within-tree interaction). It is a simplified version of the CoEvo model, which doesn't attempt to fit drift (d1/d2)")
    comment <- "Generalist Matching Mutualism model.
    \nStarts with 3 or 4 lineages having the same value X_0 ~ Normal(m0,v0).
    \nOne trait in each lineage, all lineages evolving then non-independtly according to the GMM expression.
    This was my second attempt at including geography. It assumes no competition between lineages
    within the same tree (clade), and only competition between trees.
    It does properly divide S (Pk/l) by the number of overlapping lineages (rowsums).
    (see Manceau et al. pg.559, equation 7)."
    # paramsNames <- c("m0", "v0", "d1", "d2", "S", "sigma")
    # params0 <- c(0,0,1,-1,0.5,1)
    paramsNames <- c("m0", "v0", "S", "sigma")
    params0 <- c(0,0,0.5,1)
    
    eventEndOfPeriods <- endOfPeriodsGMMgeo(tree1, tree2, geo.object)
    n <- eventEndOfPeriods$nLineages1[1] + eventEndOfPeriods$nLineages2[1] - 1
    
    initialCondition <- function(params) return( list(mean=rep(params[1], times=n), var=matrix(rep(params[2], times=n*n), nrow=n ) ) ) 
    
    aAGamma <- function(i, params){
      vectorA <- function(t) return( c( rep(params[3], times=eventEndOfPeriods$nLineages1[i]), rep(params[3], times=eventEndOfPeriods$nLineages2[i]) ) )
      matrixGamma <- function(t) return(diag(params[4], eventEndOfPeriods$nLineages1[i] + eventEndOfPeriods$nLineages2[i]))
      
      bloc1 <- diag(params[3], eventEndOfPeriods$nLineages1[i])
      # bloc1 <- params[5]*(geo.object$geography.object[[i]][(1:(eventEndOfPeriods$nLineages1[i])),
      #                                                      (1:(eventEndOfPeriods$nLineages1[i]))])
      
      #bloc2 <- matrix(rep(-params[5]/eventEndOfPeriods$nLineages2[i], times=eventEndOfPeriods$nLineages1[i]*eventEndOfPeriods$nLineages2[i]), nrow=eventEndOfPeriods$nLineages1[i])
      bloc2.int <- geo.object$geography.object[[i]][(1:(eventEndOfPeriods$nLineages1[i])),
                                                    ((eventEndOfPeriods$nLineages1[i]+1):(eventEndOfPeriods$nLineages1[i]+eventEndOfPeriods$nLineages2[i]))]
      if(is.null(nrow(bloc2.int))){
        bloc2 <- as.matrix(-params[3]*(bloc2.int/bloc2.int)); bloc2 <- as.matrix(bloc2)
      } else {
        bloc2 <- -params[3]*(sweep(bloc2.int, 2, (colSums(bloc2.int)-1), "/"))
      }
      bloc2[is.na(bloc2)] <- 0; bloc2[is.infinite(bloc2)] <- 0
      
      #bloc3 <- matrix(rep(-params[5]/eventEndOfPeriods$nLineages1[i], times=eventEndOfPeriods$nLineages1[i]*eventEndOfPeriods$nLineages2[i]), nrow=eventEndOfPeriods$nLineages2[i])
      bloc3.int <- geo.object$geography.object[[i]][((eventEndOfPeriods$nLineages1[i]+1):(eventEndOfPeriods$nLineages1[i]+eventEndOfPeriods$nLineages2[i])), 
                                                    (1:(eventEndOfPeriods$nLineages1[i]))]
      if(is.null(nrow(bloc3.int))){
        bloc3 <- t(-params[3]*(bloc3.int/bloc3.int)); bloc3 <- as.matrix(bloc3)
      } else {
        bloc3 <- -params[3]*(sweep(bloc3.int, 1, (rowSums(bloc3.int)-1), "/"))
      }
      bloc3[is.na(bloc3)] <- 0; bloc3[is.infinite(bloc3)] <- 0
      
      bloc4 <- diag(params[3], eventEndOfPeriods$nLineages2[i])
      matrixA <- rbind(cbind(bloc1, bloc2), cbind(bloc3, bloc4))
      
      
      return(list(a=vectorA, A=matrixA, Gamma=matrixGamma))
    }
    
    constraints <- function(params) return(params[2]>=0 && params[4]>=0)
    
    if( keyword == "CoEvo_driftless" ){
      model <- new(Class="PhenotypicGMM_driftless", name=keyword, period=eventEndOfPeriods$periods, aAGamma=aAGamma, numbersCopy=eventEndOfPeriods$copy, numbersPaste=eventEndOfPeriods$paste, initialCondition=initialCondition, paramsNames=paramsNames, constraints=constraints, params0=params0, tipLabels=eventEndOfPeriods$labeling, tipLabelsSimu=eventEndOfPeriods$labeling, comment=comment, n1=eventEndOfPeriods$nLineages1, n2=eventEndOfPeriods$nLineages2)
    }else{
      model <- new(Class="PhenotypicModel",         name=keyword, period=eventEndOfPeriods$periods, aAGamma=aAGamma, numbersCopy=eventEndOfPeriods$copy, numbersPaste=eventEndOfPeriods$paste, initialCondition=initialCondition, paramsNames=paramsNames, constraints=constraints, params0=params0, tipLabels=eventEndOfPeriods$labeling, tipLabelsSimu=eventEndOfPeriods$labeling, comment=comment)
    }
    class(model)[1] <- "PhenotypicModel"
  }
  
  else if(keyword == "_CoEvo_driftless_fix") {
    if (is.null(geo.object)) { stop("this model requires a geo.object")}
    print("CoEvo: A GMM model incorporating geography (correcting for the number of interacting lineages), but estimating S between trees only (no within-tree interaction). It is a simplified version of the CoEvo model, which doesn't attempt to fit drift (d1/d2)")
    comment <- "Generalist Matching Mutualism model.
    \nStarts with 3 or 4 lineages having the same value X_0 ~ Normal(m0,v0).
    \nOne trait in each lineage, all lineages evolving then non-independtly according to the GMM expression.
    This was my second attempt at including geography. It assumes no competition between lineages
    within the same tree (clade), and only competition between trees.
    It does properly divide S (Pk/l) by the number of overlapping lineages (rowsums).
    (see Manceau et al. pg.559, equation 7)."
    # paramsNames <- c("m0", "v0", "d1", "d2", "S", "sigma")
    # params0 <- c(0,0,1,-1,0.5,1)
    paramsNames <- c("m0", "v0", "S", "sigma")
    params0 <- c(0,0,0.5,1)
    
    eventEndOfPeriods <- endOfPeriodsGMMgeo(tree1, tree2, geo.object)
    n <- eventEndOfPeriods$nLineages1[1] + eventEndOfPeriods$nLineages2[1] - 1
    
    initialCondition <- function(params) return( list(mean=rep(params[1], times=n), var=matrix(rep(params[2], times=n*n), nrow=n ) ) ) 
    
    aAGamma <- function(i, params){
      vectorA <- function(t) return( c( rep(params[3], times=eventEndOfPeriods$nLineages1[i]), rep(params[3], times=eventEndOfPeriods$nLineages2[i]) ) )
      matrixGamma <- function(t) return(diag(params[4], eventEndOfPeriods$nLineages1[i] + eventEndOfPeriods$nLineages2[i]))
      
      block1 <- as.matrix(diag(0, eventEndOfPeriods$nLineages1[i]))
      block2.int <- geo.object$geography.object[[i]][(1:(eventEndOfPeriods$nLineages1[i])),
                                                     ((eventEndOfPeriods$nLineages1[i]+1):(eventEndOfPeriods$nLineages1[i]+eventEndOfPeriods$nLineages2[i]))]
      block2 <- t(t(block2.int))
      
      block3.int <- as.matrix(geo.object$geography.object[[i]][((eventEndOfPeriods$nLineages1[i]+1):(eventEndOfPeriods$nLineages1[i]+eventEndOfPeriods$nLineages2[i])), 
                                                               (1:(eventEndOfPeriods$nLineages1[i]))])
      if(ncol(block3.int)<=1){
        block3 <- t(block3.int)
      } else{block3 <- block3.int}
      
      block4 <- as.matrix(diag(0, eventEndOfPeriods$nLineages2[i]))
      
      matrixA.int <- rbind(cbind(block1, block2), cbind(block3, block4))
      matrixA.int <- -params[5]*sweep(matrixA.int, 1, rowSums(matrixA.int), "/")
      diag(matrixA.int) <- params[5]
      matrixA <- matrixA.int
      
      return(list(a=vectorA, A=matrixA, Gamma=matrixGamma))
    }
    
    constraints <- function(params) return(params[2]>=0 && params[4]>=0)
    
    if( keyword == "CoEvo_driftless_fix" ){
      model <- new(Class="PhenotypicGMM_driftless", name=keyword, period=eventEndOfPeriods$periods, aAGamma=aAGamma, numbersCopy=eventEndOfPeriods$copy, numbersPaste=eventEndOfPeriods$paste, initialCondition=initialCondition, paramsNames=paramsNames, constraints=constraints, params0=params0, tipLabels=eventEndOfPeriods$labeling, tipLabelsSimu=eventEndOfPeriods$labeling, comment=comment, n1=eventEndOfPeriods$nLineages1, n2=eventEndOfPeriods$nLineages2)
    }else{
      model <- new(Class="PhenotypicModel",         name=keyword, period=eventEndOfPeriods$periods, aAGamma=aAGamma, numbersCopy=eventEndOfPeriods$copy, numbersPaste=eventEndOfPeriods$paste, initialCondition=initialCondition, paramsNames=paramsNames, constraints=constraints, params0=params0, tipLabels=eventEndOfPeriods$labeling, tipLabelsSimu=eventEndOfPeriods$labeling, comment=comment)
    }
    #class(model)[1] <- "PhenotypicModel"
  }
  
  else if(keyword == "_CoEvo_driftless_split") {
    if (is.null(geo.object)) { stop("this model requires a geo.object")}
    print("CoEvo: A GMM model incorporating geography (correcting for the number of interacting lineages), but estimating S between trees only (no within-tree interaction). It is a simplified version of the CoEvo model, which doesn't attempt to fit drift (d1/d2)")
    comment <- "Generalist Matching Mutualism model.
    \nStarts with 3 or 4 lineages having the same value X_0 ~ Normal(m0,v0).
    \nOne trait in each lineage, all lineages evolving then non-independtly according to the GMM expression.
    This was my second attempt at including geography. It assumes no competition between lineages
    within the same tree (clade), and only competition between trees.
    It does properly divide S (Pk/l) by the number of overlapping lineages (rowsums).
    (see Manceau et al. pg.559, equation 7)."
    # paramsNames <- c("m0", "v0", "d1", "d2", "S", "sigma")
    # params0 <- c(0,0,1,-1,0.5,1)
    paramsNames <- c("m0", "v0", "S1", "sigma", "S2")
    params0 <- c(0,0,0.5,1,0.4)
    
    eventEndOfPeriods <- endOfPeriodsGMMgeo(tree1, tree2, geo.object)
    n <- eventEndOfPeriods$nLineages1[1] + eventEndOfPeriods$nLineages2[1] - 1
    
    initialCondition <- function(params) return( list(mean=rep(params[1], times=n), var=matrix(rep(params[2], times=n*n), nrow=n ) ) ) 
    
    aAGamma <- function(i, params){
      vectorA <- function(t) return( c( rep(params[3], times=eventEndOfPeriods$nLineages1[i]), rep(params[3], times=eventEndOfPeriods$nLineages2[i]) ) )
      matrixGamma <- function(t) return(diag(params[4], eventEndOfPeriods$nLineages1[i] + eventEndOfPeriods$nLineages2[i]))
      
      # bloc1 <- diag(params[3], eventEndOfPeriods$nLineages1[i])
      # bloc1 <- params[5]*(geo.object$geography.object[[i]][(1:(eventEndOfPeriods$nLineages1[i])),
      #                                                      (1:(eventEndOfPeriods$nLineages1[i]))])
      bloc1.int <- (geo.object$geography.object[[i]][(1:(eventEndOfPeriods$nLineages1[i])),
                                                     (1:(eventEndOfPeriods$nLineages1[i]))])
      if(is.null(nrow(bloc1.int))){
        bloc1 <- bloc1.int/bloc1.int;
      } else {
        bloc1 <- params[3]*(sweep(bloc1.int, 1, (rowSums(bloc1.int)-1), "/"))
        bloc1[is.na(bloc1)]<-0; bloc1[is.infinite(bloc1)]<-0
      }
      
      #bloc2 <- matrix(rep(-params[5]/eventEndOfPeriods$nLineages2[i], times=eventEndOfPeriods$nLineages1[i]*eventEndOfPeriods$nLineages2[i]), nrow=eventEndOfPeriods$nLineages1[i])
      bloc2.int <- geo.object$geography.object[[i]][(1:(eventEndOfPeriods$nLineages1[i])),
                                                    ((eventEndOfPeriods$nLineages1[i]+1):(eventEndOfPeriods$nLineages1[i]+eventEndOfPeriods$nLineages2[i]))]
      if(is.null(nrow(bloc2.int))){
        bloc2 <- as.matrix(params[5]*(bloc2.int/bloc2.int));
      } else {
        bloc2 <- params[5]*(sweep(bloc2.int, 1, (rowSums(bloc2.int)), "/"))
      }
      bloc2[is.na(bloc2)] <- 0; bloc2[is.infinite(bloc2)] <- 0
      
      #bloc3 <- matrix(rep(-params[5]/eventEndOfPeriods$nLineages1[i], times=eventEndOfPeriods$nLineages1[i]*eventEndOfPeriods$nLineages2[i]), nrow=eventEndOfPeriods$nLineages2[i])
      bloc3.int <- geo.object$geography.object[[i]][((eventEndOfPeriods$nLineages1[i]+1):(eventEndOfPeriods$nLineages1[i]+eventEndOfPeriods$nLineages2[i])), 
                                                    (1:(eventEndOfPeriods$nLineages1[i]))]
      if(is.null(nrow(bloc3.int))){
        bloc3 <- (t(params[5]*(bloc3.int/bloc3.int))); bloc3 <- as.matrix(bloc3) 
      } else {
        bloc3 <- params[5]*(sweep(bloc3.int, 1, (rowSums(bloc3.int)), "/"))
      }
      bloc3[is.na(bloc3)] <- 0; bloc3[is.infinite(bloc3)] <- 0
      
      #bloc4 <- diag(params[3], eventEndOfPeriods$nLineages2[i])
      bloc4.int <- geo.object$geography.object[[i]][((eventEndOfPeriods$nLineages1[i]+1):(eventEndOfPeriods$nLineages1[i]+eventEndOfPeriods$nLineages2[i])),
                                                    ((eventEndOfPeriods$nLineages1[i]+1):(eventEndOfPeriods$nLineages1[i]+eventEndOfPeriods$nLineages2[i]))]
      if(is.null(nrow(bloc4.int))){
        bloc4 <- bloc4.int/bloc4.int; bloc4 <- as.matrix(bloc4)
      } else {
        bloc4 <- params[3]*(sweep(bloc4.int, 1, (rowSums(bloc4.int)-1), "/"))
        bloc4[is.na(bloc4)]<-0; bloc4[is.infinite(bloc4)]<-0
      } 
      
      #matrixA <- rbind(cbind(bloc1, bloc2), cbind(bloc3, bloc4))
      matrixA <- rbind(cbind(bloc1, bloc2), cbind(bloc3, bloc4))
      diag(matrixA) <- 0; diag(matrixA) <- abs(rowSums(matrixA))
      #matrixA <- -matrixA
      
      return(list(a=vectorA, A=matrixA, Gamma=matrixGamma))
    }
    
    constraints <- function(params) return(params[2]>=0 && params[4]>=0)
    
    if( keyword == "CoEvo_driftless_split" ){
      model <- new(Class="PhenotypicGMM_driftless_split", name=keyword, period=eventEndOfPeriods$periods, aAGamma=aAGamma, numbersCopy=eventEndOfPeriods$copy, numbersPaste=eventEndOfPeriods$paste, initialCondition=initialCondition, paramsNames=paramsNames, constraints=constraints, params0=params0, tipLabels=eventEndOfPeriods$labeling, tipLabelsSimu=eventEndOfPeriods$labeling, comment=comment, n1=eventEndOfPeriods$nLineages1, n2=eventEndOfPeriods$nLineages2)
    }else{
      model <- new(Class="PhenotypicModel", name=keyword, period=eventEndOfPeriods$periods, aAGamma=aAGamma, numbersCopy=eventEndOfPeriods$copy, numbersPaste=eventEndOfPeriods$paste, initialCondition=initialCondition, paramsNames=paramsNames, constraints=constraints, params0=params0, tipLabels=eventEndOfPeriods$labeling, tipLabelsSimu=eventEndOfPeriods$labeling, comment=comment)
    }
    class(model)[1] <- "PhenotypicModel"
  }
  
  else if(keyword == "_CoEvo_3S") {
    if (is.null(geo.object)) { stop("this model requires a geo.object")}
    print("CoEvo: A GMM model incorporating geography (correcting for the number of interacting lineages), but estimating S between trees only (no within-tree interaction). It is a simplified version of the CoEvo model, which doesn't attempt to fit drift (d1/d2)")
    comment <- "Generalist Matching Mutualism model.
    \nStarts with 3 or 4 lineages having the same value X_0 ~ Normal(m0,v0).
    \nOne trait in each lineage, all lineages evolving then non-independtly according to the GMM expression.
    This was my second attempt at including geography. It assumes no competition between lineages
    within the same tree (clade), and only competition between trees.
    It does properly divide S (Pk/l) by the number of overlapping lineages (rowsums).
    (see Manceau et al. pg.559, equation 7)."
    # paramsNames <- c("m0", "v0", "d1", "d2", "S", "sigma")
    # params0 <- c(0,0,1,-1,0.5,1)
    paramsNames <- c("m0", "v0", "S1", "sigma", "S2", "S3")
    params0 <- c(0,0,0.5,1,0.4,0.2)
    
    eventEndOfPeriods <- endOfPeriodsGMMgeo(tree1, tree2, geo.object)
    n <- eventEndOfPeriods$nLineages1[1] + eventEndOfPeriods$nLineages2[1] - 1
    
    initialCondition <- function(params) return( list(mean=rep(params[1], times=n), var=matrix(rep(params[2], times=n*n), nrow=n ) ) ) 
    
    aAGamma <- function(i, params){
      vectorA <- function(t) return( c( rep(params[3], times=eventEndOfPeriods$nLineages1[i]), rep(params[3], times=eventEndOfPeriods$nLineages2[i]) ) )
      matrixGamma <- function(t) return(diag(params[4], eventEndOfPeriods$nLineages1[i] + eventEndOfPeriods$nLineages2[i]))
      
      # bloc1 <- diag(params[3], eventEndOfPeriods$nLineages1[i])
      # bloc1 <- params[5]*(geo.object$geography.object[[i]][(1:(eventEndOfPeriods$nLineages1[i])),
      #                                                      (1:(eventEndOfPeriods$nLineages1[i]))])
      bloc1.int <- (geo.object$geography.object[[i]][(1:(eventEndOfPeriods$nLineages1[i])),
                                                     (1:(eventEndOfPeriods$nLineages1[i]))])
      if(is.null(nrow(bloc1.int))){
        bloc1 <- bloc1.int/bloc1.int
      } else {
        bloc1 <- params[3]*(sweep(bloc1.int, 1, (rowSums(bloc1.int)-1), "/"))
        bloc1[is.na(bloc1)]<-0; bloc1[is.infinite(bloc1)]<-0
      }
      
      #bloc2 <- matrix(rep(-params[5]/eventEndOfPeriods$nLineages2[i], times=eventEndOfPeriods$nLineages1[i]*eventEndOfPeriods$nLineages2[i]), nrow=eventEndOfPeriods$nLineages1[i])
      bloc2.int <- geo.object$geography.object[[i]][(1:(eventEndOfPeriods$nLineages1[i])),
                                                    ((eventEndOfPeriods$nLineages1[i]+1):(eventEndOfPeriods$nLineages1[i]+eventEndOfPeriods$nLineages2[i]))]
      if(is.null(nrow(bloc2.int))){
        bloc2 <- as.matrix(t(params[5]*(bloc2.int/bloc2.int)))
      } else {
        bloc2 <- params[5]*(sweep(bloc2.int, 1, (rowSums(bloc2.int)), "/"))
      }
      bloc2[is.na(bloc2)] <- 0; bloc2[is.infinite(bloc2)] <- 0
      
      #bloc3 <- matrix(rep(-params[5]/eventEndOfPeriods$nLineages1[i], times=eventEndOfPeriods$nLineages1[i]*eventEndOfPeriods$nLineages2[i]), nrow=eventEndOfPeriods$nLineages2[i])
      bloc3.int <- geo.object$geography.object[[i]][((eventEndOfPeriods$nLineages1[i]+1):(eventEndOfPeriods$nLineages1[i]+eventEndOfPeriods$nLineages2[i])), 
                                                    (1:(eventEndOfPeriods$nLineages1[i]))]
      if(is.null(nrow(bloc3.int))){
        bloc3 <- t(t(params[5]*(bloc3.int/bloc3.int)))
      } else {
        bloc3 <- params[5]*(sweep(bloc3.int, 1, (rowSums(bloc3.int)), "/"))
      }
      bloc3[is.na(bloc3)] <- 0; bloc3[is.infinite(bloc3)] <- 0
      
      #bloc4 <- diag(params[3], eventEndOfPeriods$nLineages2[i])
      bloc4.int <- geo.object$geography.object[[i]][((eventEndOfPeriods$nLineages1[i]+1):(eventEndOfPeriods$nLineages1[i]+eventEndOfPeriods$nLineages2[i])),
                                                    ((eventEndOfPeriods$nLineages1[i]+1):(eventEndOfPeriods$nLineages1[i]+eventEndOfPeriods$nLineages2[i]))]
      if(is.null(nrow(bloc4.int))){
        bloc4 <- bloc4.int/bloc4.int
      } else {
        bloc4 <- params[6]*(sweep(bloc4.int, 1, (rowSums(bloc4.int)-1), "/"))
        bloc4[is.na(bloc4)]<-0; bloc4[is.infinite(bloc4)]<-0
      } 
      
      #matrixA <- rbind(cbind(bloc1, bloc2), cbind(bloc3, bloc4))
      matrixA <- rbind(cbind(bloc1, bloc2), cbind(bloc3, bloc4))
      diag(matrixA) <- 0; diag(matrixA) <- abs(rowSums(matrixA))
      #matrixA <- -matrixA
      
      return(list(a=vectorA, A=matrixA, Gamma=matrixGamma))
    }
    
    constraints <- function(params) return(params[2]>=0 && params[4]>=0)
    
    if( keyword == "CoEvo_3S" ){
      model <- new(Class="PhenotypicGMM_driftless", name=keyword, period=eventEndOfPeriods$periods, aAGamma=aAGamma, numbersCopy=eventEndOfPeriods$copy, numbersPaste=eventEndOfPeriods$paste, initialCondition=initialCondition, paramsNames=paramsNames, constraints=constraints, params0=params0, tipLabels=eventEndOfPeriods$labeling, tipLabelsSimu=eventEndOfPeriods$labeling, comment=comment, n1=eventEndOfPeriods$nLineages1, n2=eventEndOfPeriods$nLineages2)
    }else{
      model <- new(Class="PhenotypicModel", name=keyword, period=eventEndOfPeriods$periods, aAGamma=aAGamma, numbersCopy=eventEndOfPeriods$copy, numbersPaste=eventEndOfPeriods$paste, initialCondition=initialCondition, paramsNames=paramsNames, constraints=constraints, params0=params0, tipLabels=eventEndOfPeriods$labeling, tipLabelsSimu=eventEndOfPeriods$labeling, comment=comment)
    }
    #class(model) <- "PhenotypicModel"
  }
  
  else if(keyword == "_CoEvo_Split") {
    if (is.null(geo.object)) { stop("this model requires a geo.object")}
    
    print("CoEvo_Split: A GMM model incorporating geography (correcting for the number of interacting lineages), and estimating two S values. S1 is the strength of interaction between lineages within trees, S2 is the strength of interaction between lineages among trees. This model isn't implemented correctly though, because S is different for the direction of interaction (goanna vs mammal = S1; mammal vs. goanna = S2). Needs to be fixed!")
    comment <- "Generalist Matching Mutualism model.\nStarts with 3 or 4 lineages having the same value X_0 ~ Normal(m0,v0).
    \nOne trait in each lineage, all lineages evolving then non-independtly according to the GMM expression.
    This was my second attempt at including competition within a tree, as well as between trees. 
    It estimates a single S (S2) value across the pair of trees (competition between trees), as 
    as well as a second S value (S1) for competition within trees.
    It also properly divides S (Pk,l) by the number of overlapping lineages (rowsums).
    (see Manceau et al. pg.559, equation 7)."
    paramsNames <- c("m0", "v0", "d1", "d2", "S1", "sigma", "S2")
    params0 <- c(0,0,1,-1,0.5,1,0.5)
    
    eventEndOfPeriods <- endOfPeriodsGMMgeo(tree1, tree2, geo.object)
    n <- eventEndOfPeriods$nLineages1[1] + eventEndOfPeriods$nLineages2[1] - 1
    
    initialCondition <- function(params) return( list(mean=rep(params[1], times=n), var=matrix(rep(params[2], times=n*n), nrow=n ) ) ) 
    
    aAGamma <- function(i, params){
      vectorA <- function(t) return( c( rep(params[3]*params[5], times=eventEndOfPeriods$nLineages1[i]), rep(params[4]*params[5], times=eventEndOfPeriods$nLineages2[i]) ) )
      matrixGamma <- function(t) return(diag(params[6], eventEndOfPeriods$nLineages1[i] + eventEndOfPeriods$nLineages2[i]))
      
      block1 <- (geo.object$geography.object[[i]][(1:(eventEndOfPeriods$nLineages1[i])),
                                                     (1:(eventEndOfPeriods$nLineages1[i]))])
      block2.int <- geo.object$geography.object[[i]][(1:(eventEndOfPeriods$nLineages1[i])),
                                                     ((eventEndOfPeriods$nLineages1[i]+1):(eventEndOfPeriods$nLineages1[i]+eventEndOfPeriods$nLineages2[i]))]
      block2 <- t(t(block2.int))
      block3.int <- as.matrix(geo.object$geography.object[[i]][((eventEndOfPeriods$nLineages1[i]+1):(eventEndOfPeriods$nLineages1[i]+eventEndOfPeriods$nLineages2[i])), 
                                                               (1:(eventEndOfPeriods$nLineages1[i]))])
      if(ncol(block3.int)<=1){
        block3 <- t(block3.int)
      } else{block3 <- block3.int}
      
      block4.int <- geo.object$geography.object[[i]][((eventEndOfPeriods$nLineages1[i]+1):(eventEndOfPeriods$nLineages1[i]+eventEndOfPeriods$nLineages2[i])),
                                                    ((eventEndOfPeriods$nLineages1[i]+1):(eventEndOfPeriods$nLineages1[i]+eventEndOfPeriods$nLineages2[i]))]
      if(is.null(nrow(block4.int))){
        block4 <- as.matrix(block4.int)
      } else{block4 <- block4.int}
      
      matrixA.int <- rbind(cbind(block1, block2), cbind(block3, block4))
      matrixA.int <- -params[5]*sweep(matrixA.int, 1, rowSums(matrixA.int)-1, "/")
      diag(matrixA.int) <- params[5]
      matrixA <- matrixA.int
      
      return(list(a=vectorA, A=matrixA, Gamma=matrixGamma))
    }
    
    constraints <- function(params) return(params[2]>=0 && params[6]>=0) # && params[5]<=0 && params[7]<=0)
    
    if( keyword == "CoEvo_Split" ){
      model <- new(Class="PhenotypicGMM", name=keyword, period=eventEndOfPeriods$periods, aAGamma=aAGamma, numbersCopy=eventEndOfPeriods$copy, numbersPaste=eventEndOfPeriods$paste, initialCondition=initialCondition, paramsNames=paramsNames, constraints=constraints, params0=params0, tipLabels=eventEndOfPeriods$labeling, tipLabelsSimu=eventEndOfPeriods$labeling, comment=comment, n1=eventEndOfPeriods$nLineages1, n2=eventEndOfPeriods$nLineages2)
    }else{
      model <- new(Class="PhenotypicModel", name=keyword, period=eventEndOfPeriods$periods, aAGamma=aAGamma, numbersCopy=eventEndOfPeriods$copy, numbersPaste=eventEndOfPeriods$paste, initialCondition=initialCondition, paramsNames=paramsNames, constraints=constraints, params0=params0, tipLabels=eventEndOfPeriods$labeling, tipLabelsSimu=eventEndOfPeriods$labeling, comment=comment)
    }
    #class(model)[1] <- "PhenotypicModel"
  }
  
  else if(keyword == "_JointPM") {
    #if (is.null(geo.object)) { stop("this model requires a geo.object")}
    
    comment <- "This Model incorporates the PM model, estimated jointly across two trees.
    Thus, it estimates only a single S (attraction/repulsion) parameter, and assumes each 
    lineage interacts with ALL other lineages (all A-matrix off-diagonal values are the same)."
    paramsNames <- c("m0", "v0", "d1", "d2", "S", "sigma")
    params0 <- c(0,0,1,-1,0.5,1)
    
    #eventEndOfPeriods <- endOfPeriodsGMMgeo(tree1, tree2, geo.object)
    eventEndOfPeriods <- endOfPeriodsGMM(tree1, tree2)
    n <- eventEndOfPeriods$nLineages1[1] + eventEndOfPeriods$nLineages2[1] - 1
    
    initialCondition <- function(params) return( list(mean=rep(params[1], times=n), var=matrix(rep(params[2], times=n*n), nrow=n ) ) ) 
    
    aAGamma <- function(i, params){
      vectorA <- function(t) return( c( rep(params[3]*params[5], times=eventEndOfPeriods$nLineages1[i]), rep(params[4]*params[5], times=eventEndOfPeriods$nLineages2[i]) ) )
      matrixGamma <- function(t) return(diag(params[6], eventEndOfPeriods$nLineages1[i] + eventEndOfPeriods$nLineages2[i]))
      
      matrixA <- -params[5]*(matrix(nrow=(eventEndOfPeriods$nLineages1[i]+eventEndOfPeriods$nLineages2[i]),
                                    ncol=(eventEndOfPeriods$nLineages1[i]+eventEndOfPeriods$nLineages2[i]),1))
      diag(matrixA) <- params[5]
      
      return(list(a=vectorA, A=matrixA, Gamma=matrixGamma))
    }
    
    constraints <- function(params) return(params[2]>=0 && params[6]>=0)
    
    if( keyword == "JointPM" ){
      model <- new(Class="PhenotypicGMM", name=keyword, period=eventEndOfPeriods$periods, aAGamma=aAGamma, numbersCopy=eventEndOfPeriods$copy, numbersPaste=eventEndOfPeriods$paste, initialCondition=initialCondition, paramsNames=paramsNames, constraints=constraints, params0=params0, tipLabels=eventEndOfPeriods$labeling, tipLabelsSimu=eventEndOfPeriods$labeling, comment=comment, n1=eventEndOfPeriods$nLineages1, n2=eventEndOfPeriods$nLineages2)
    }else{
      model <- new(Class="PhenotypicModel", name=keyword, period=eventEndOfPeriods$periods, aAGamma=aAGamma, numbersCopy=eventEndOfPeriods$copy, numbersPaste=eventEndOfPeriods$paste, initialCondition=initialCondition, paramsNames=paramsNames, constraints=constraints, params0=params0, tipLabels=eventEndOfPeriods$labeling, tipLabelsSimu=eventEndOfPeriods$labeling, comment=comment)
    }
    class(model)[1] <- "PhenotypicModel"
  }
  
  else if(keyword == "_JointPMOUless") {
    #if (is.null(geo.object)) { stop("this model requires a geo.object")}
    
    comment <- "This Model incorporates the PM model, estimated jointly across two trees.
    Thus, it estimates only a single S (attraction/repulsion) parameter, and assumes each 
    lineage interacts with ALL other lineages (all A-matrix off-diagonal values are the same)."
    paramsNames <- c("m0", "v0", "S", "sigma")
    params0 <- c(0,0,0.5,1)
    
    #eventEndOfPeriods <- endOfPeriodsGMMgeo(tree1, tree2, geo.object)
    eventEndOfPeriods <- endOfPeriodsGMM(tree1, tree2)
    n <- eventEndOfPeriods$nLineages1[1] + eventEndOfPeriods$nLineages2[1] - 1
    
    initialCondition <- function(params) return( list(mean=rep(params[1], times=n), var=matrix(rep(params[2], times=n*n), nrow=n ) ) ) 
    
    aAGamma <- function(i, params){
      vectorA <- function(t) return( rep(params[3], sum(eventEndOfPeriods$nLineages1[i] + eventEndOfPeriods$nLineages2[i])) )
      matrixGamma <- function(t) return(diag(params[4], eventEndOfPeriods$nLineages1[i] + eventEndOfPeriods$nLineages2[i]))
      
      matrixA <- -params[3]*(matrix(nrow=(eventEndOfPeriods$nLineages1[i]+eventEndOfPeriods$nLineages2[i]),
                                    ncol=(eventEndOfPeriods$nLineages1[i]+eventEndOfPeriods$nLineages2[i]),1))
      diag(matrixA) <- params[3]
      
      return(list(a=vectorA, A=matrixA, Gamma=matrixGamma))
    }
    
    constraints <- function(params) return(params[2]>=0 && params[4]>=0)
    
    if( keyword == "JointPMOUless" ){
      model <- new(Class="PhenotypicGMM", name=keyword, period=eventEndOfPeriods$periods, aAGamma=aAGamma, numbersCopy=eventEndOfPeriods$copy, numbersPaste=eventEndOfPeriods$paste, initialCondition=initialCondition, paramsNames=paramsNames, constraints=constraints, params0=params0, tipLabels=eventEndOfPeriods$labeling, tipLabelsSimu=eventEndOfPeriods$labeling, comment=comment, n1=eventEndOfPeriods$nLineages1, n2=eventEndOfPeriods$nLineages2)
    }else{
      model <- new(Class="PhenotypicModel", name=keyword, period=eventEndOfPeriods$periods, aAGamma=aAGamma, numbersCopy=eventEndOfPeriods$copy, numbersPaste=eventEndOfPeriods$paste, initialCondition=initialCondition, paramsNames=paramsNames, constraints=constraints, params0=params0, tipLabels=eventEndOfPeriods$labeling, tipLabelsSimu=eventEndOfPeriods$labeling, comment=comment)
    }
    class(model)[1] <- "PhenotypicModel"
  }
  
  else if(keyword == "_JointPM_geo") {
    if (is.null(geo.object)) { stop("this model requires a geo.object")}
    
    comment <- "This Model incorporates the PM model, estimated jointly across two trees, but
    also includes an aspect of geography. Thus, it estimates only a single S (attraction/repulsion) parameter, 
    but assumes that this interaction occurs only with co-occuring lineages (from geo.object).
    However, all interacting lineages are treated equally (all A-matrix off-diagonal values are the same).
    Really, this isn't correct, because we need to correct S for the number of co-occuring lineages (Pk,l)"
    paramsNames <- c("m0", "v0", "d1", "d2", "S", "sigma")
    params0 <- c(0,0,1,-1,0.5,1)

    eventEndOfPeriods <- endOfPeriodsGMMgeo(tree1, tree2, geo.object)
    n <- eventEndOfPeriods$nLineages1[1] + eventEndOfPeriods$nLineages2[1] - 1
    
    initialCondition <- function(params) return( list(mean=rep(params[1], times=n), var=matrix(rep(params[2], times=n*n), nrow=n ) ) ) 
    
    aAGamma <- function(i, params){
      vectorA <- function(t) return( c( rep(params[3]*params[5], times=eventEndOfPeriods$nLineages1[i]), rep(params[4]*params[5], times=eventEndOfPeriods$nLineages2[i]) ) )
      matrixGamma <- function(t) return(diag(params[6], eventEndOfPeriods$nLineages1[i] + eventEndOfPeriods$nLineages2[i]))
      
      matrixA <- -params[5]*geo.object$geography.object[[i]]
      diag(matrixA) <- params[5]

      return(list(a=vectorA, A=matrixA, Gamma=matrixGamma))
    }
    
    constraints <- function(params) return(params[2]>=0 && params[6]>=0)
    
    if( keyword == "JointPM_geo" ){
      model <- new(Class="PhenotypicGMM", name=keyword, period=eventEndOfPeriods$periods, aAGamma=aAGamma, numbersCopy=eventEndOfPeriods$copy, numbersPaste=eventEndOfPeriods$paste, initialCondition=initialCondition, paramsNames=paramsNames, constraints=constraints, params0=params0, tipLabels=eventEndOfPeriods$labeling, tipLabelsSimu=eventEndOfPeriods$labeling, comment=comment, n1=eventEndOfPeriods$nLineages1, n2=eventEndOfPeriods$nLineages2)
    }else{
      model <- new(Class="PhenotypicModel", name=keyword, period=eventEndOfPeriods$periods, aAGamma=aAGamma, numbersCopy=eventEndOfPeriods$copy, numbersPaste=eventEndOfPeriods$paste, initialCondition=initialCondition, paramsNames=paramsNames, constraints=constraints, params0=params0, tipLabels=eventEndOfPeriods$labeling, tipLabelsSimu=eventEndOfPeriods$labeling, comment=comment)
    }
    class(model)[1] <- "PhenotypicModel"
  }
  
  else if(keyword == "_JointS") {
    if (is.null(geo.object)) { stop("this model requires a geo.object")}
    
    print("This model is a GMM/CoEvo model that estimates ")
    comment <- "I haven't actually done this yet, it should be simple though!"
    paramsNames <- c("m0", "v0", "d1", "d2", "S", "sigma")
    params0 <- c(0,0,1,-1,0.5,1)
    
    eventEndOfPeriods <- endOfPeriodsGMMgeo(tree1, tree2, geo.object)
    n <- eventEndOfPeriods$nLineages1[1] + eventEndOfPeriods$nLineages2[1] - 1
    
    initialCondition <- function(params) return( list(mean=rep(params[1], times=n), var=matrix(rep(params[2], times=n*n), nrow=n ) ) ) 
    
    aAGamma <- function(i, params){
      vectorA <- function(t) return( c( rep(params[3]*params[5], times=eventEndOfPeriods$nLineages1[i]), rep(params[4]*params[5], times=eventEndOfPeriods$nLineages2[i]) ) )
      matrixGamma <- function(t) return(diag(params[6], eventEndOfPeriods$nLineages1[i] + eventEndOfPeriods$nLineages2[i]))
      
      block <- -params[5]*geo.object$geography.object[[i]]
      matrixA <- sweep(block, 1, (rowSums(geo.object$geography.object[[i]])-1), "/")      # I think all I'd need to do is:
      #diag(matrixA) <- params[5]
      diag(matrixA) <- 0; diag(matrixA) <- -(rowSums(matrixA))
      matrixA <- -matrixA

      return(list(a=vectorA, A=matrixA, Gamma=matrixGamma))
    }
    
    constraints <- function(params) return(params[2]>=0 && params[6]>=0)
    
    if( keyword == "JointS" ){
      model <- new(Class="PhenotypicGMM", name=keyword, period=eventEndOfPeriods$periods, aAGamma=aAGamma, numbersCopy=eventEndOfPeriods$copy, numbersPaste=eventEndOfPeriods$paste, initialCondition=initialCondition, paramsNames=paramsNames, constraints=constraints, params0=params0, tipLabels=eventEndOfPeriods$labeling, tipLabelsSimu=eventEndOfPeriods$labeling, comment=comment, n1=eventEndOfPeriods$nLineages1, n2=eventEndOfPeriods$nLineages2)
    }else{
      model <- new(Class="PhenotypicModel", name=keyword, period=eventEndOfPeriods$periods, aAGamma=aAGamma, numbersCopy=eventEndOfPeriods$copy, numbersPaste=eventEndOfPeriods$paste, initialCondition=initialCondition, paramsNames=paramsNames, constraints=constraints, params0=params0, tipLabels=eventEndOfPeriods$labeling, tipLabelsSimu=eventEndOfPeriods$labeling, comment=comment)
    }
    class(model)[1] <- "PhenotypicModel"
  }
  
  else if(keyword == "_coBM"){
    if (!is.null(geo.object)) { print("ignoring geo.object, using GMM instead")}
    
    comment <- "Coevolution model estimating a single Brownian Motion rate for two trees"
    paramsNames <- c("m0", "v0", "d1", "d2", "S", "sigma")
    params0 <- c(0,0,1,-1,0.5,1)
    
    eventEndOfPeriods <- endOfPeriodsGMM(tree1, tree2)
    n <- eventEndOfPeriods$nLineages1[1] + eventEndOfPeriods$nLineages2[1] - 1
    
    initialCondition <- function(params) return( list(mean=rep(params[1], times=n), var=matrix(rep(params[2], times=n*n), nrow=n ) ) ) 
    
    aAGamma <- function(i, params){
      vectorA <- function(t) return( c( rep(params[3]*params[5], times=eventEndOfPeriods$nLineages1[i]), rep(params[4]*params[5], times=eventEndOfPeriods$nLineages2[i]) ) )
      matrixGamma <- function(t) return(diag(params[6], eventEndOfPeriods$nLineages1[i] + eventEndOfPeriods$nLineages2[i]))
      
      bloc1 <- diag(params[5], eventEndOfPeriods$nLineages1[i])
      bloc2 <- matrix(rep(-params[5]/eventEndOfPeriods$nLineages2[i], times=eventEndOfPeriods$nLineages1[i]*eventEndOfPeriods$nLineages2[i]), nrow=eventEndOfPeriods$nLineages1[i])
      bloc3 <- matrix(rep(-params[5]/eventEndOfPeriods$nLineages1[i], times=eventEndOfPeriods$nLineages1[i]*eventEndOfPeriods$nLineages2[i]), nrow=eventEndOfPeriods$nLineages2[i])
      bloc4 <- diag(params[5], eventEndOfPeriods$nLineages2[i])
      matrixA <- rbind(cbind(bloc1, bloc2), cbind(bloc3, bloc4))
      # I think all I'd need to do is:
      #matrixA <- matrixA * geo.object$geography.object[[i]]
      
      return(list(a=vectorA, A=matrixA, Gamma=matrixGamma))
    } 
    
    constraints <- function(params) return(params[2]>=0 && params[6]>=0)
    
    if( keyword == "GMM" ){
      model <- new(Class="PhenotypicGMM", name=keyword, c, aAGamma=aAGamma, numbersCopy=eventEndOfPeriods$copy, numbersPaste=eventEndOfPeriods$paste, initialCondition=initialCondition, paramsNames=paramsNames, constraints=constraints, params0=params0, tipLabels=eventEndOfPeriods$labeling, tipLabelsSimu=eventEndOfPeriods$labeling, comment=comment, n1=eventEndOfPeriods$nLineages1, n2=eventEndOfPeriods$nLineages2)
    }else{
      model <- new(Class="PhenotypicModel", name=keyword, period=eventEndOfPeriods$periods, aAGamma=aAGamma, numbersCopy=eventEndOfPeriods$copy, numbersPaste=eventEndOfPeriods$paste, initialCondition=initialCondition, paramsNames=paramsNames, constraints=constraints, params0=params0, tipLabels=eventEndOfPeriods$labeling, tipLabelsSimu=eventEndOfPeriods$labeling, comment=comment)
    }
    class(model)[1] <- "PhenotypicModel"
  }
  
  return(model)
}

createGeoModelCoevolution <- function(tree1, tree2, geo.object, keyword = "GMM"){
  
  if(keyword == "GMM" || keyword == "GMMbis"){
    
    comment <- "Generalist Matching Mutualism model.\nStarts with 3 or 4 lineages having the same value X_0 ~ Normal(m0,v0).\nOne trait in each lineage, all lineages evolving then non-independtly according to the GMM expression."
    paramsNames <- c("m0", "v0", "d1", "d2", "S", "sigma")
    params0 <- c(0,0,1,-1,0.5,1)

# NEED TO CREATE A NEW 'resortGeoObject' FUNCTION FOR HANDLING A geo.object COMPOSED OF TWO TREES    
# EASIER SOLUTION WOULD BE TO FIX THE CONFUSION BETWEEN THE '$nLineages2[1]' OBJECT (1), WHEN IT SHOULD BE ZERO
# WHEN i=1, AND WE TRY TO MULTIPLY 'matrixA.int' BY 'geo.object$geography.object[[1]]' THEY'RE NONCONFORMBALE
# 'matrixA.int' IS 3X3 AND 'geo.object...' IS 2X2!
# MAYBE MAKE A STEM LINEAGE FOR THE YOUNGER TREE, UNTIL IT SPLITS
    #resgeo.object <- resortGeoObject(tree, geo.object)
    eventEndOfPeriods <- endOfPeriodsGMMgeo(tree1, tree2, geo.object)
    n <- eventEndOfPeriods$nLineages1[1] + eventEndOfPeriods$nLineages2[1] - 1
    
    initialCondition <- function(params) return( list(mean=rep(params[1], times=n), var=matrix(rep(params[2], times=n*n), nrow=n ) ) ) 
    
    aAGamma <- function(i, params){
      vectorA <- function(t) return( c( rep(params[3]*params[5], times=eventEndOfPeriods$nLineages1[i]), rep(params[4]*params[5], times=eventEndOfPeriods$nLineages2[i]) ) )
      matrixGamma <- function(t) return(diag(params[6], eventEndOfPeriods$nLineages1[i] + eventEndOfPeriods$nLineages2[i]))
      
      bloc1 <- diag(params[5], eventEndOfPeriods$nLineages1[i])
      bloc2 <- matrix(rep(-params[5]/eventEndOfPeriods$nLineages2[i], times=eventEndOfPeriods$nLineages1[i]*eventEndOfPeriods$nLineages2[i]), nrow=eventEndOfPeriods$nLineages1[i])
      bloc3 <- matrix(rep(-params[5]/eventEndOfPeriods$nLineages1[i], times=eventEndOfPeriods$nLineages1[i]*eventEndOfPeriods$nLineages2[i]), nrow=eventEndOfPeriods$nLineages2[i])
      bloc4 <- diag(params[5], eventEndOfPeriods$nLineages2[i])
      matrixA.int <- rbind(cbind(bloc1, bloc2), cbind(bloc3, bloc4))
      # I think all I'd need to do is:
      matrixA <- matrixA.int * geo.object$geography.object[[i]]
      
      return(list(a=vectorA, A=matrixA, Gamma=matrixGamma))
    }
    
    constraints <- function(params) return(params[2]>=0 && params[6]>=0)
    
    if( keyword == "GMM" ){
      model <- new(Class="PhenotypicGMM", name=keyword, period=eventEndOfPeriods$periods, aAGamma=aAGamma, numbersCopy=eventEndOfPeriods$copy, numbersPaste=eventEndOfPeriods$paste, initialCondition=initialCondition, paramsNames=paramsNames, constraints=constraints, params0=params0, tipLabels=eventEndOfPeriods$labeling, tipLabelsSimu=eventEndOfPeriods$labeling, comment=comment, n1=eventEndOfPeriods$nLineages1, n2=eventEndOfPeriods$nLineages2)
    }else{
      model <- new(Class="PhenotypicModel", name=keyword, period=eventEndOfPeriods$periods, aAGamma=aAGamma, numbersCopy=eventEndOfPeriods$copy, numbersPaste=eventEndOfPeriods$paste, initialCondition=initialCondition, paramsNames=paramsNames, constraints=constraints, params0=params0, tipLabels=eventEndOfPeriods$labeling, tipLabelsSimu=eventEndOfPeriods$labeling, comment=comment)
    }
    
  }
  
  return(model)
}



##################################################
#    Describe the periods on a 'phylo' tree
##################################################

getMatrixCoalescenceJ <- function(tree, periods){
  # The entry (k,l) of the matrix is the index j such that tau_j = t_{k,l}
  matrixCoalescenceTimes <- findMRCA(tree, type="height")
  n <- length(matrixCoalescenceTimes[,1])
  matrixCoalescenceJ <- diag(0, n)
  for(k in 1:n){
    for(l in 1:n){
      matrixCoalescenceJ[k,l] <- which(periods == matrixCoalescenceTimes[k,l])
    }
  }
  
  return(matrixCoalescenceJ)
}

isATip <- function(tree, branch_number){
  return(!(tree$edge[branch_number,2] %in% tree$edge[,1]))
}

periodizeOneTree <- function(tree){
  # Returns 3 vectors giving 
  # 1) the periods of the tree, 
  # 2) the starting times of all branches in the tree 
  # 3) the death time of all branches in the tree
  
  nodeheight <- nodeHeights(tree)
  startingTimes <- round(nodeheight[,1],6)
  endTimes <- round(nodeheight[,2], 6)
  all_time_events <- sort(c(startingTimes, endTimes))
  # the following removes identical entries in the vector
  periods <- unique(all_time_events)
  
  return(list(periods=periods, startingTimes=startingTimes, endTimes=endTimes))
}

periodizeOneTree_geo <- function(tree,geo.object){
  # Returns 3 vectors giving 
  # 1) the periods of the tree, 
  # 2) the starting times of all branches in the tree 
  # 3) the death time of all branches in the tree
  hold<-nodeHeights(tree)
  startingTimes <- hold[,1]
  endTimes <- hold[,2]
  all_time_events <- sort(c(startingTimes, endTimes))
  
  nodetimes=max(branching.times(tree))-sort(branching.times(tree),decreasing=TRUE)
  extv<-vapply(geo.object$geography.object,function(x)dim(x)[1],1)
  outv<-c(1)
  for(i in 2:length(extv)){
    if(extv[i]!=extv[i-1]){
      outv<-c(outv,i)
    }}
  
  chg.times=which(!1:length(geo.object$times)%in%c(outv,length(geo.object$times)))
  periods=sort(c(geo.object$times[chg.times],unique(startingTimes),max(endTimes)))
  return(list(periods=periods, startingTimes=startingTimes, endTimes=endTimes))
}

endOfPeriods <- function(periodizing, tree){
  # Returns the list of branching or dying lineages at the beginning of each period : copy
  # Together with the list of places where the new lineage is inserted (or zero if a lineage dies) : paste
  # And the number of lineages on the focal period : nLineages
  # The rule is : at each branching point, the first of the two new branches is assigned its mother label, and the new one takes the last label (n, where n is the number of lineages at that time)
  
  nBranch <- length(periodizing$startingTimes)
  nPeriods <- length(periodizing$periods)
  
  numbersCopy <- rep(0, times=nPeriods)
  numbersPaste <- rep(0, times=nPeriods)
  numbersLineages <- rep(0, times=nPeriods)
  numbersLivingLineages <- rep(0, times=nPeriods)
  
  # We initialize the labeling of branches in the tree
  labelingLineages <- rep(0, times=nBranch)
  initialBranches <- periodizing$startingTimes[periodizing$startingTimes==0]
  if(length(initialBranches) == 1){
    labelingLineages[1] <- 1
    n <- 1
  }else{
    labelingLineages[periodizing$startingTimes==0] <- c(1,2)
    n <- 2
  }
  numbersLineages[1] <- n
  numbersLivingLineages[1] <- n
  numbersCopy[1] <- 1
  numbersPaste[1] <- 2
  
  for(i in 2:nPeriods){
    tau_i <- periodizing$periods[i]
    newBranches <- which(tau_i == periodizing$startingTimes)
    # If tau_i is a birth time on the tree
    if(length(newBranches) == 2){
      n <- n+1
      labelingLineages[newBranches[1]] <- labelingLineages[newBranches[1]-1]
      labelingLineages[newBranches[2]] <- n
      numbersCopy[i] <- labelingLineages[newBranches[1]-1]
      numbersPaste[i] <- n
      numbersLivingLineages[i] <- numbersLivingLineages[i-1]+1
      # Else, tau_i is only a death time of one or many terminal branches.
    }else{
      deadBranches <- which(tau_i == periodizing$endTimes)
      numbersCopy[i] <- labelingLineages[ deadBranches[1] ]
      numbersPaste[i] <- 0
      numbersLivingLineages[i] <- numbersLivingLineages[i-1]-1
    }
    numbersLineages[i] <- n
  }
  
  permutationLabels <- labelingLineages[!(periodizing$endTimes %in% periodizing$startingTimes)]
  labeling <- tree$tip.label[order(permutationLabels)]
  
  return(list(copy=numbersCopy, paste=numbersPaste, nLineages=numbersLineages, labeling=labeling, nLivingLineages=numbersLivingLineages))
}

getLivingLineages <- function(i, eventEndOfPeriods){
  
  livingLineages <- rep(1, times=eventEndOfPeriods$nLineages[i])
  deads <- eventEndOfPeriods$copy[1:i][eventEndOfPeriods$paste[1:i] == 0]
  livingLineages[deads] <- 0
  
  return(livingLineages)
}

endOfPeriodsGMM <- function(tree1, tree2){
  # Warning !! This has to be used on ultrametric trees only ! It could be extended to take into account the death of lineages, though.
  # Returns the list of branching lineages at the beginning of each period : copy
  # Together with the list of places where the new lineage is inserted : paste
  # And the number of lineages in the first and second tree on the focal period : nLineages1, nLineages2
  # The rule is : If a lineage in the first tree gives birth, the first of the two new branches is assigned its mother label, and the new one takes the label nLineages1. All other lineages are then pushed back. If a lineage in the second tree gives birth, the first of the two new branches is assigned its mother label, and the last one takes the last label (nLineages1 + nLineages2)
  
  periodizing1 <- periodizeOneTree(tree1)
  periodizing2 <- periodizeOneTree(tree2)
  nBranch1 <- length(periodizing1$startingTimes)
  nBranch2 <- length(periodizing2$startingTimes)
  nPeriods <- length(periodizing1$periods) + length(periodizing2$periods) - 1
  
  numbersCopy <- rep(0, times=nPeriods)
  numbersPaste <- rep(0, times=nPeriods)
  numbersLineages1 <- rep(0, times=nPeriods)
  numbersLineages2 <- rep(0, times=nPeriods)
  periods <- rep(0, times=nPeriods)
  
  # We initialize the labeling of branches in the tree
  labelingLineages1 <- rep(0, times=nBranch1)
  labelingLineages2 <- rep(0, times=nBranch2)
  
  # The highest tree starts with two lineages (crown) the other one starts with one (root)
  Tmax1 <- periodizing1$periods[length(periodizing1$periods)]
  Tmax2 <- periodizing2$periods[length(periodizing2$periods)]
  if( Tmax1 < Tmax2 ){
    n1 <- 1
    n2 <- 2
    labelingLineages1[1] <- c(1)
    labelingLineages2[periodizing2$startingTimes==0] <- c(1,2)
    periodizing1$periods <- periodizing1$periods + (Tmax2-Tmax1)
    periodizing1$startingTimes <- periodizing1$startingTimes + (Tmax2-Tmax1)
    periodizing1$endTimes <- periodizing1$endTimes + (Tmax2-Tmax1)
    numbersCopy[1] <- n2
    numbersPaste[1] <- n2+n1
  }else if( Tmax2 < Tmax1 ){
    n1 <- 2
    n2 <- 1
    labelingLineages1[periodizing1$startingTimes==0] <- c(1,2)
    labelingLineages2[1] <- c(1)
    periodizing2$periods <- periodizing2$periods + (Tmax1-Tmax2)
    periodizing2$startingTimes <- periodizing2$startingTimes + (Tmax1-Tmax2)
    periodizing2$endTimes <- periodizing2$endTimes + (Tmax1-Tmax2)
    numbersCopy[1] <- 1
    numbersPaste[1] <- 2
  }else{
    n1 <- 2
    n2 <- 2
    labelingLineages1[periodizing1$startingTimes==0] <- c(1,2)
    labelingLineages2[periodizing2$startingTimes==0] <- c(1,2)
    numbersCopy[1] <- 1
    numbersPaste[1] <- 2
  }
  numbersLineages1[1] <- n1
  numbersLineages2[1] <- n2
  
  for(i in 2:(nPeriods-1)){
    tau_i1 <- periodizing1$periods[n1]
    tau_i2 <- periodizing2$periods[n2]
    if( tau_i1 < tau_i2 ){
      n1 <- n1 +1
      newBranches <- which(tau_i1 == periodizing1$startingTimes)
      if(n1 > 2){
        labelingLineages1[newBranches[1]] <- labelingLineages1[newBranches[1]-1]
        numbersCopy[i] <- labelingLineages1[newBranches[1]-1]
      }else{
        numbersCopy[i] <- 1
      }
      labelingLineages1[newBranches[2]] <- n1
      numbersPaste[i] <- n1
      periods[i] <- tau_i1
    }else{
      n2 <- n2 +1
      newBranches <- which(tau_i2 == periodizing2$startingTimes)
      if(n2 > 2){
        labelingLineages2[newBranches[1]] <- labelingLineages2[newBranches[1]-1]
        numbersCopy[i] <- n1 + labelingLineages2[newBranches[1]-1]
      }else{
        numbersCopy[i] <- n1 + 1
      }
      labelingLineages2[newBranches[2]] <- n2
      numbersPaste[i] <- n1+n2
      periods[i] <- tau_i2
    }
    numbersLineages1[i] <- n1
    numbersLineages2[i] <- n2
  }
  
  permutationLabels1 <- labelingLineages1[!(periodizing1$endTimes %in% periodizing1$startingTimes)]
  labeling1 <- tree1$tip.label[order(permutationLabels1)]
  permutationLabels2 <- labelingLineages2[!(periodizing2$endTimes %in% periodizing2$startingTimes)]
  labeling2 <- tree2$tip.label[order(permutationLabels2)]
  labeling <- c(labeling1, labeling2)
  
  periods[nPeriods] <- max(Tmax1, Tmax2)
  
  return(list(periods=periods, copy=numbersCopy, paste=numbersPaste, nLineages1=numbersLineages1, nLineages2=numbersLineages2, labeling=labeling))
}

endOfPeriodsGMMgeo <- function(tree1, tree2, geo.object){
  # Warning !! This has to be used on ultrametric trees only ! It could be extended to take into account the death of lineages, though.
  # Returns the list of branching lineages at the beginning of each period : copy
  # Together with the list of places where the new lineage is inserted : paste
  # And the number of lineages in the first and second tree on the focal period : nLineages1, nLineages2
  # The rule is : If a lineage in the first tree gives birth, the first of the two new branches is assigned its mother label, and the new one takes the label nLineages1. All other lineages are then pushed back. If a lineage in the second tree gives birth, the first of the two new branches is assigned its mother label, and the last one takes the last label (nLineages1 + nLineages2)
  
  periodizing1 <- periodizeOneTree_geo(tree1, geo.object)
  periodizing2 <- periodizeOneTree_geo(tree2, geo.object)
  nBranch1 <- length(periodizing1$startingTimes)
  nBranch2 <- length(periodizing2$startingTimes)
  nPeriods <- length(periodizing1$periods) + length(periodizing2$periods) - 1
  
  numbersCopy <- rep(0, times=nPeriods)
  numbersPaste <- rep(0, times=nPeriods)
  numbersLineages1 <- rep(0, times=nPeriods)
  numbersLineages2 <- rep(0, times=nPeriods)
  periods <- rep(0, times=nPeriods)
  
  # We initialize the labeling of branches in the tree
  labelingLineages1 <- rep(0, times=nBranch1)
  labelingLineages2 <- rep(0, times=nBranch2)
  
  # The highest tree starts with two lineages (crown) the other one starts with one (root)
  Tmax1 <- periodizing1$periods[length(periodizing1$periods)]
  Tmax2 <- periodizing2$periods[length(periodizing2$periods)]
  if( Tmax1 < Tmax2 ){
    n1 <- 1
    n2 <- 2
    labelingLineages1[1] <- c(1)
    labelingLineages2[periodizing2$startingTimes==0] <- c(1,2)
    periodizing1$periods <- periodizing1$periods + (Tmax2-Tmax1)
    periodizing1$startingTimes <- periodizing1$startingTimes + (Tmax2-Tmax1)
    periodizing1$endTimes <- periodizing1$endTimes + (Tmax2-Tmax1)
    numbersCopy[1] <- n2
    numbersPaste[1] <- n2+n1
  }else if( Tmax2 < Tmax1 ){
    n1 <- 2
    n2 <- 1
    labelingLineages1[periodizing1$startingTimes==0] <- c(1,2)
    labelingLineages2[1] <- c(1)
    periodizing2$periods <- periodizing2$periods + (Tmax1-Tmax2)
    periodizing2$startingTimes <- periodizing2$startingTimes + (Tmax1-Tmax2)
    periodizing2$endTimes <- periodizing2$endTimes + (Tmax1-Tmax2)
    numbersCopy[1] <- 1
    numbersPaste[1] <- 2
  }else{
    n1 <- 2
    n2 <- 2
    labelingLineages1[periodizing1$startingTimes==0] <- c(1,2)
    labelingLineages2[periodizing2$startingTimes==0] <- c(1,2)
    numbersCopy[1] <- 1
    numbersPaste[1] <- 2
  }
  numbersLineages1[1] <- n1
  numbersLineages2[1] <- n2
  
  for(i in 2:(nPeriods-1)){
    tau_i1 <- periodizing1$periods[n1]
    tau_i2 <- periodizing2$periods[n2]
    if( tau_i1 < tau_i2 ){
      n1 <- n1 +1
      newBranches <- which(tau_i1 == periodizing1$startingTimes)
      if(n1 > 2){
        labelingLineages1[newBranches[1]] <- labelingLineages1[newBranches[1]-1]
        numbersCopy[i] <- labelingLineages1[newBranches[1]-1]
      }else{
        numbersCopy[i] <- 1
      }
      labelingLineages1[newBranches[2]] <- n1
      numbersPaste[i] <- n1
      periods[i] <- tau_i1
    }else{
      n2 <- n2 +1
      newBranches <- which(tau_i2 == periodizing2$startingTimes)
      if(n2 > 2){
        labelingLineages2[newBranches[1]] <- labelingLineages2[newBranches[1]-1]
        numbersCopy[i] <- n1 + labelingLineages2[newBranches[1]-1]
      }else{
        numbersCopy[i] <- n1 + 1
      }
      labelingLineages2[newBranches[2]] <- n2
      numbersPaste[i] <- n1+n2
      periods[i] <- tau_i2
    }
    numbersLineages1[i] <- n1
    numbersLineages2[i] <- n2
  }
  
  permutationLabels1 <- labelingLineages1[!(periodizing1$endTimes %in% periodizing1$startingTimes)]
  labeling1 <- tree1$tip.label[order(permutationLabels1)]
  permutationLabels2 <- labelingLineages2[!(periodizing2$endTimes %in% periodizing2$startingTimes)]
  labeling2 <- tree2$tip.label[order(permutationLabels2)]
  labeling <- c(labeling1, labeling2)
  
  periods[nPeriods] <- max(Tmax1, Tmax2)
  
  return(list(periods=periods, copy=numbersCopy, paste=numbersPaste, nLineages1=numbersLineages1, nLineages2=numbersLineages2, labeling=labeling))
}

createModel_MC <- function(tree){
  comment <- "Matching competition model\n Implemented as in Drury et al. Systematic Biology."
  print("MC: Matching competition model as implemented in Drury et al. (2015)")
  paramsNames <- c("m0","logsigma","S")
  params0 <- c(0,log(1),0)
  
  periodizing <- periodizeOneTree(tree) 
  eventEndOfPeriods <- endOfPeriods(periodizing, tree)
  
  initialCondition <- function(params) return( list(mean=c(params[1]), var=matrix(c(0))) ) 
  
  ###is this where the A matrix incorporating geography needs to go? if so, what is the order in which lineage sympatry data need to be introduced
  
  aAGamma <- function(i, params){
    vectorU <- getLivingLineages(i, eventEndOfPeriods)
    vectorA <- function(t) return(0*vectorU)
    matrixGamma <- function(t) return(exp(params[2])*diag(vectorU))
    matrixA <- params[3]*diag(vectorU) - (params[3]/sum(vectorU)) * outer(vectorU,vectorU) 
    return(list(a=vectorA, A=matrixA, Gamma=matrixGamma))
  }
  constraints <- function(params) return(params[3]<=0)
  model <- new(Class="PhenotypicADiag", name="MC", period=periodizing$periods, aAGamma=aAGamma, numbersCopy=eventEndOfPeriods$copy, numbersPaste=eventEndOfPeriods$paste, initialCondition=initialCondition, paramsNames=paramsNames, constraints=constraints, params0=params0, tipLabels=eventEndOfPeriods$labeling,  comment=comment)
  return(model)
}

createModel_MC_geo <- function(tree,geo.object, keyword = "MC+geo"){
  
  if (keyword == "MC+geo" || keyword == "MC_geo") {
    comment <- "Matching competition model with biogeography\n Implemented as in Drury et al. Systematic Biology."
    print("MC+geo: Matching competition model with biogeography")
    paramsNames <- c("m0","logsigma","S")
    params0 <- c(0,log(1),0)
    
    periodizing <- periodizeOneTree_geo(tree,geo.object) 
    eventEndOfPeriods <- endOfPeriods(periodizing, tree)
    
    initialCondition <- function(params) return( list(mean=c(params[1]), var=matrix(c(0))) ) 
    
    ###is this where the A matrix incorporating geography needs to go? if so, what is the order in which lineage sympatry data need to be introduced
    # I think I can incorporate the geo matrix and "fix" this function at the NIJ step, then the matrix A step
    
    aAGamma <- function(i, params){
      vectorU <- getLivingLineages(i, eventEndOfPeriods)
      vectorA <- function(t) return(0*vectorU)
      matrixGamma <- function(t) return(exp(params[2])*diag(vectorU))
      nij <- colSums(geo.object$geography.object[[i]])
      matrixA <- params[3]*diag(vectorU) - (geo.object$geography.object[[i]]*(params[3]/nij))
      return(list(a=vectorA, A=matrixA, Gamma=matrixGamma))
    }
    constraints <- function(params) return(params[3]<=0)
    model <- new(Class="PhenotypicModel", name="MC+geo", period=periodizing$periods, aAGamma=aAGamma, numbersCopy=eventEndOfPeriods$copy, numbersPaste=eventEndOfPeriods$paste, initialCondition=initialCondition, paramsNames=paramsNames, constraints=constraints, params0=params0, tipLabels=eventEndOfPeriods$labeling,  comment=comment)
    return(model)
  }
  else if (keyword == "MC_geo2") {
    comment <- "Matching competition model with biogeography\n Implemented as in Drury et al. Systematic Biology."
    paramsNames <- c("m0","logsigma","S")
    params0 <- c(0,log(1),0)
    
    periodizing <- periodizeOneTree_geo(tree,geo.object) 
    eventEndOfPeriods <- endOfPeriods(periodizing, tree)
    
    initialCondition <- function(params) return( list(mean=c(params[1]), var=matrix(c(0))) ) 
    
    ###is this where the A matrix incorporating geography needs to go? if so, what is the order in which lineage sympatry data need to be introduced
    # I think I can incorporate the geo matrix and "fix" this function at the NIJ step, then the matrix A step
    
    aAGamma <- function(i, params){
      vectorU <- getLivingLineages(i, eventEndOfPeriods)
      vectorA <- function(t) return(0*vectorU)
      matrixGamma <- function(t) return(exp(params[2])*diag(vectorU))
      nij <- colSums(geo.object$geography.object[[i]])
      matrixA <- params[3]*diag(vectorU) - (geo.object$geography.object[[i]]*(params[3]/nij))
      return(list(a=vectorA, A=matrixA, Gamma=matrixGamma))
    }
    constraints <- function(params) return(params[3]<=0)
    model <- new(Class="PhenotypicModel", name="MC+geo2", period=periodizing$periods, aAGamma=aAGamma, numbersCopy=eventEndOfPeriods$copy, numbersPaste=eventEndOfPeriods$paste, initialCondition=initialCondition, paramsNames=paramsNames, constraints=constraints, params0=params0, tipLabels=eventEndOfPeriods$labeling,  comment=comment)
    return(model)
  }
  
}

createModel_PM_geo <- function(tree,geo.object,keyword){
  comment <- "Phenotype Matching model with biogeography.\nStarts with two lineages having the same value X_0 ~ Normal(m0,v0).\nOne trait in each lineage, all lineages evolving then non-independtly according to the Phenotype Matching expression."
  paramsNames <- c("m0", "v0", "theta", "psi", "S", "sigma")
  params0 <- c(0,0,0,0.2,0.5,1)
  
  periodizing <- periodizeOneTree_geo(tree, geo.object)
  eventEndOfPeriods <- endOfPeriods(periodizing, tree)
  
  initialCondition <- function(params) return( list(mean=c(params[1]), var=matrix(c(params[2]))) ) 
  
  aAGamma <- function(i, params){
    vectorU <- getLivingLineages(i, eventEndOfPeriods)
    vectorA <- function(t) return(params[3]*params[4]*vectorU)
    matrixGamma <- function(t) return(params[6]*diag(vectorU))
    nij <- colSums(geo.object$geography.object[[i]])
    # matrixA <- (params[4]+params[5])*diag(vectorU) - (params[5]/sum(vectorU)) * outer(vectorU,vectorU) # from the PM model
    # matrixA <- params[3]*diag(vectorU) -(geo.object$geography.object[[i]]*(params[3]/nij)) # from the MC_geo model
    # matrixA <- (params[4]+params[5])*diag(vectorU) - (geo.object$geography.object[[i]]*((params[4]+params[5])/nij)) # first attempt
    matrixA <- (params[4]+params[5])*diag(vectorU) - (params[5]/nij) * geo.object$geography.object[[i]] # second attempt
    matrixA <- 

    return(list(a=vectorA, A=matrixA, Gamma=matrixGamma, u=vectorU, OU=TRUE))
  }
  
  constraints <- function(params) return(params[2]>=0 && params[6]>=0)
  
  if( keyword == "PM" ){
    model <- new(Class="PhenotypicPM", name="PM+geo", period=periodizing$periods, aAGamma=aAGamma, numbersCopy=eventEndOfPeriods$copy, numbersPaste=eventEndOfPeriods$paste, initialCondition=initialCondition, paramsNames=paramsNames, constraints=constraints, params0=params0, tipLabels=eventEndOfPeriods$labeling, tipLabelsSimu=eventEndOfPeriods$labeling, comment=comment)
  }else if( keyword == "PMbis" ){
    model <- new(Class="PhenotypicADiag", name=keyword, period=periodizing$periods, aAGamma=aAGamma, numbersCopy=eventEndOfPeriods$copy, numbersPaste=eventEndOfPeriods$paste, initialCondition=initialCondition, paramsNames=paramsNames, constraints=constraints, params0=params0, tipLabels=eventEndOfPeriods$labeling, tipLabelsSimu=eventEndOfPeriods$labeling, comment=comment)
  }else{
    model <- new(Class="PhenotypicModel", name=keyword, period=periodizing$periods, aAGamma=aAGamma, numbersCopy=eventEndOfPeriods$copy, numbersPaste=eventEndOfPeriods$paste, initialCondition=initialCondition, paramsNames=paramsNames, constraints=constraints, params0=params0, tipLabels=eventEndOfPeriods$labeling, tipLabelsSimu=eventEndOfPeriods$labeling, comment=comment)
  }
  
  return(model)
#####################################################################  ##
  #comment <- "Phylogenetic Matching model with biogeography"
  #paramsNames <- c("m0","logsigma","S")
  #params0 <- c(0,log(1),0)
  #
  #periodizing <- periodizeOneTree_geo(tree,geo.object) 
  #eventEndOfPeriods <- endOfPeriods(periodizing, tree)
  #
  #initialCondition <- function(params) return( list(mean=c(params[1]), var=matrix(c(0))) ) 
  #
  ####is this where the A matrix incorporating geography needs to go? if so, what is the order in which lineage sympatry data need to be introduced
  #
  #aAGamma <- function(i, params){
  #  vectorU <- getLivingLineages(i, eventEndOfPeriods)
  #  vectorA <- function(t) return(0*vectorU)
  #  matrixGamma <- function(t) return(exp(params[2])*diag(vectorU))
  #  nij <- colSums(geo.object$geography.object[[i]])
  #  matrixA <- params[3]*diag(vectorU) -(geo.object$geography.object[[i]]*(params[3]/nij))
  #  return(list(a=vectorA, A=matrixA, Gamma=matrixGamma))
  #}
  #constraints <- function(params) return(params[3]<=0)
  #model <- new(Class="PhenotypicModel", name="MC_geo", period=periodizing$periods, aAGamma=aAGamma, numbersCopy=eventEndOfPeriods$copy, numbersPaste=eventEndOfPeriods$paste, initialCondition=initialCondition, paramsNames=paramsNames, constraints=constraints, params0=params0, tipLabels=eventEndOfPeriods$labeling,  comment=comment)
  #return(model)
}

resortGeoObject<-function(phylo1, geo.object, phylo2=NULL){
  gmat<-geo.object$geography.object
  if (is.null(phylo2)) {
    phylo <- phylo1
    if(any(grepl("___",phylo$tip.label))|any(grepl("-",phylo$tip.label))|any(grepl("/",phylo$tip.label))){stop("script will not work with '___', '-', '+', '*','/', or '^' in any tip labels; remove these characters")}
    paste(rep(LETTERS,each=26),LETTERS,sep="")->TWOLETTERS
    paste(rep(TWOLETTERS,each=26),LETTERS,sep="")->THREELETTERS
    nodeDist<-vector(mode = "numeric", length = phylo$Nnode)
    totlen<-length(phylo$tip.label)
    root <-totlen  + 1
    heights<-nodeHeights(phylo)
    for (i in 1:dim(phylo$edge)[1]){
      nodeDist[[phylo$edge[i, 1] - totlen]] <- heights[i]
    }
    nodeDist<-c(nodeDist,max(heights))
    nodeDiff<-diff(nodeDist)
    flag=0
    if(sum(nodeDiff<0)>0){  ##this loop renumbers the nodes if trees nodes are not placed in sequential order
      node.order<-match(rank(heights[,1],ties.method="min"),seq(1, by = 2, len = phylo$Nnode))
      node.order<-node.order+totlen
      old.edge<-phylo$edge
      old.phylo<-phylo
      phylo$edge[,1]<-node.order
      for(j in 1:length(phylo$edge[,2])){
        if(phylo$edge[j,2]>totlen){
          #match number order in old edge
          #lookup value in new edge
          #replace with value
          phylo$edge[j,2]<-phylo$edge[,1][match(phylo$edge[j,2],old.edge[,1])]
        }
      }
      nodeDist<-vector()
      for (i in 1:dim(phylo$edge)[1]){
        nodeDist[[phylo$edge[i, 1] - totlen]] <- heights[i]
      }
      nodeDist<-c(nodeDist,max(heights))
      nodeDiff<-diff(nodeDist)
      flag=1
    }
    
    
    tips<-1:length(phylo$tip.label) 
    order.vec<-rep(0,length=dim(phylo$edge)[1])
    
    counter=1
    for(i in (totlen+1):(totlen+phylo$Nnode)){
      for(j in 1:2){
        if(order.vec[which(phylo$edge[,1]==i)][j]==0){
          order.vec[which(phylo$edge[,1]==i)][j]<-counter
          m<-phylo$edge[which(phylo$edge[,1]==i)[j],2]
          while(!m%in%tips){
            order.vec[which(phylo$edge[,1]==m)][1]<-counter
            m<-phylo$edge[which(phylo$edge[,1]==m)[1],2]
          }
          counter=counter+1 
        } 
      }
    }
    
    
    newedge<-cbind(phylo$edge,order.vec)
    
    mat<-matrix(nrow=0, ncol=4)
    counter_three_letters <- 0
    for(i in 1:phylo$Nnode){
      other<-newedge[newedge[,1]==i+totlen, 2]
      for(b in other){
        int<-matrix(ncol=4)
        int[1]<-i+totlen
        if(b>totlen){
          counter_three_letters <- counter_three_letters + 1
          int[2]<-paste(".",THREELETTERS[counter_three_letters],sep="")
          int[3]<-b
          int[4]<-newedge[which(newedge[,2]==b),3]
        } else {
          int[2]<-phylo$tip.label[b]
          int[3]<-b
          int[4]<-newedge[which(newedge[,2]==b),3] 
        }
        mat<-rbind(mat,int)
      }
    }
    
    sorted.gmat<-list()
    for(i in 1:length(gmat)){
      sorted.gmat[[i]]<-gmat[[i]][order(as.numeric(mat[match(rownames(gmat[[i]]),mat[,2]),4])),order(as.numeric(mat[match(rownames(gmat[[i]]),mat[,2]),4]))]
    }
    
    return(list(geography.object=sorted.gmat,times=geo.object$times,spans=geo.object$spans))
    
  } else {
    if(any(grepl("___",phylo$tip.label))|any(grepl("-",phylo$tip.label))|any(grepl("/",phylo$tip.label))){stop("script will not work with '___', '-', '+', '*','/', or '^' in any tip labels; remove these characters")}
    paste(rep(LETTERS,each=26),LETTERS,sep="")->TWOLETTERS
    paste(rep(TWOLETTERS,each=26),LETTERS,sep="")->THREELETTERS
    nodeDist<-vector(mode = "numeric", length = phylo$Nnode)
    totlen<-length(phylo$tip.label)
    root <-totlen  + 1
    heights<-nodeHeights(phylo)
    for (i in 1:dim(phylo$edge)[1]){
      nodeDist[[phylo$edge[i, 1] - totlen]] <- heights[i]
    }
    nodeDist<-c(nodeDist,max(heights))
    nodeDiff<-diff(nodeDist)
    flag=0
    if(sum(nodeDiff<0)>0){  ##this loop renumbers the nodes if trees nodes are not placed in sequential order
      node.order<-match(rank(heights[,1],ties.method="min"),seq(1, by = 2, len = phylo$Nnode))
      node.order<-node.order+totlen
      old.edge<-phylo$edge
      old.phylo<-phylo
      phylo$edge[,1]<-node.order
      for(j in 1:length(phylo$edge[,2])){
        if(phylo$edge[j,2]>totlen){
          #match number order in old edge
          #lookup value in new edge
          #replace with value
          phylo$edge[j,2]<-phylo$edge[,1][match(phylo$edge[j,2],old.edge[,1])]
        }
      }
      nodeDist<-vector()
      for (i in 1:dim(phylo$edge)[1]){
        nodeDist[[phylo$edge[i, 1] - totlen]] <- heights[i]
      }
      nodeDist<-c(nodeDist,max(heights))
      nodeDiff<-diff(nodeDist)
      flag=1
    }
    
    
    tips<-1:length(phylo$tip.label) 
    order.vec<-rep(0,length=dim(phylo$edge)[1])
    
    counter=1
    for(i in (totlen+1):(totlen+phylo$Nnode)){
      for(j in 1:2){
        if(order.vec[which(phylo$edge[,1]==i)][j]==0){
          order.vec[which(phylo$edge[,1]==i)][j]<-counter
          m<-phylo$edge[which(phylo$edge[,1]==i)[j],2]
          while(!m%in%tips){
            order.vec[which(phylo$edge[,1]==m)][1]<-counter
            m<-phylo$edge[which(phylo$edge[,1]==m)[1],2]
          }
          counter=counter+1 
        } 
      }
    }
    
    
    newedge<-cbind(phylo$edge,order.vec)
    
    mat<-matrix(nrow=0, ncol=4)
    counter_three_letters <- 0
    for(i in 1:phylo$Nnode){
      other<-newedge[newedge[,1]==i+totlen, 2]
      for(b in other){
        int<-matrix(ncol=4)
        int[1]<-i+totlen
        if(b>totlen){
          counter_three_letters <- counter_three_letters + 1
          int[2]<-paste(".",THREELETTERS[counter_three_letters],sep="")
          int[3]<-b
          int[4]<-newedge[which(newedge[,2]==b),3]
        } else {
          int[2]<-phylo$tip.label[b]
          int[3]<-b
          int[4]<-newedge[which(newedge[,2]==b),3] 
        }
        mat<-rbind(mat,int)
      }
    }
    
    sorted.gmat<-list()
    for(i in 1:length(gmat)){
      sorted.gmat[[i]]<-gmat[[i]][order(as.numeric(mat[match(rownames(gmat[[i]]),mat[,2]),4])),order(as.numeric(mat[match(rownames(gmat[[i]]),mat[,2]),4]))]
    }
    
    
    return(list(geography.object=sorted.gmat,times=geo.object$times,spans=geo.object$spans))
  }
  
}






###########################################################