################
# Fitting the model only a single time appears a bit unpredictable
# The final implementation should fit the models multiple times, I'll keep them up here:

# Let's make a function that uses mclapply to fit the model a number of times with different starting params
# start by creating sets of plausible starting parameters from across the surface
search.surface <- function(model, n.iter = 10, traits, n.proc = 8, no.S=1, results=c("best", "all"), start.params=NULL) {
  beginning <- Sys.time()
  
  if(is.null(start.params)){
    init.params <- lapply(1:n.iter, function(x) {
      c(m0 = mean(traits), 
        v0 = runif(1, min = 1e-10, max = 0.01),
        d1 = runif(1, min=1e-10, max=0.1),
        d2 = -runif(1, min=1e-10, max=0.1),
        S = rnorm(1, 0, 0.25),
        sigma = runif(1, min=1e-10, max=0.1))
    })
  } else { init.params <- rep(list(start.params), n.iter)}
  
  
  if(no.S==2){init.params <- lapply(1:n.iter, function(x) init.params[[x]] = append(init.params[[x]], rnorm(1,0,0.25)))}
  
  if(model@name=="PM_OUless") {init.params <- lapply(1:n.iter, function(x) init.params[[x]] <- init.params[[x]][c("m0", "v0", "S", "sigma")])}
  if(model@name=="ACDC") {init.params <- lapply(1:n.iter, function(x) init.params[[x]] <- init.params[[x]][c("m0", "v0", "sigma", "S")])}
  
  res.list <- mclapply(1:n.iter, function(x) {
    fitTipData(model, traits, GLSstyle=T, params0 = init.params[[x]])}, mc.cores = n.proc)
  
  res.list <- res.list[order(sapply(res.list, '[[', 1))] # sort the model fits if you'd like
  
  res.conv <- unlist(lapply(res.list, function(x) x$convergence)) # make a vector of the convergence diagnostics, so we can get the index number to remove bad ones
  res.list[which(!res.conv==0)] <- NULL # remove any model fits that didn't converge
  
  if(length(res.list)==0) {res.list[[1]] <- fitTipData(model, traits, GLSstyle = T, params0 = model@params0)}
  
  res.values <- unlist(lapply(res.list, function(x) x$value)) # make a vector of the values, so we can get the index number of the best
  #res.values[which(abs(res.values) > abs(5 * median(res.values)))] <- NA # remove any nonsensical model fits from the options
  #res.values[which(abs(res.values) < median(res.values)/5)] <- NA
  best.res <- res.list[[which.min(res.values)]]
  
  end <- Sys.time()
  duration <- format(end-beginning)
  print(paste("Computation time to fit the", model@name, "model from", n.iter, "starting points:", duration))
  
  if(results=="all"){return(list(all.results = res.list, best.result = best.res))}
  else if(results=="best"){return(best.res)}
  
}
# The function lets us control a few things via the commands:
# model: just the model of interest, you have to have built it already
# n.iter: the number of model fittings you'd like completed, defaults to 10
# traits: the input traits for model fitting
# n.proc: number of processors. this function will fit the model in parallel.
# no.S: defaults to 1. if your model requires estimating/fitting more than 1 S parameter, say so
# results: would you like the function to report just the best fitting run, or all the results from each fit attempt