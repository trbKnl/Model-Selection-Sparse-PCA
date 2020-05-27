###########

rm(list=ls())
library(testPackage)
library(doParallel)
source("./modelSelection.R")
source("./tucker.R")
source("./makeData.R")


######## Script to run the simulation study.

#given the selected model determine the quality of the estimation
determineQuality <- function(datObject, comdis, comdisspar, p, ncomp, Wstart, groups, tuningSettings) {
    
    estimation <- newAlgoCpp(datObject$X, 
               ridge = tuningSettings[[1]],
               lasso = tuningSettings[[2]],
               constraints = matrix(1, p, ncomp),
               grouplasso = tuningSettings[[3]],
               elitistlasso = tuningSettings[[4]],
               groups = groups,
               Q = ncomp,
               itr = 100000,
               nStarts = 1,
               printLoss = FALSE,
               coorDec = FALSE,
               Wstart = Wstart)

    estimationRes <- estimation
    tuckerCongruenceResW <- tuckerCongruence(estimation$W,  datObject$P[, 1:ncomp])
    tuckerCongruenceResP <- tuckerCongruence(estimation$P,  datObject$P[, 1:ncomp])
    corClassRes <- corClass(W = datObject$P, estimation$W, comdisspar = comdisspar)
    comIdentifiedRes <- comIdentified(datObject$P, estimation$W, comdis, groups)
    disIdentifiedRes <-  disIdentified(datObject$P, estimation$W, comdis)

    return(list(estimationRes         = estimationRes,
                tuckerCongruenceResW  = tuckerCongruenceResW,
                tuckerCongruenceResP  = tuckerCongruenceResP,
                corClassRes           = corClassRes,
                comIdentifiedRes      = comIdentifiedRes,
                disIdentifiedRes      = disIdentifiedRes))
}

# Simulation study skeleton
simulation_study  <-  function(reps, n, p, ncomp, groups, comdis, sparsity,
                               varianceOfComps, error, ridgeSeq, lassoSeq,
                               grouplassoSeq, elitistlassoSeq) {

    datObjectRes <- vector("list", length = reps)
    tuningISandBICRes <- vector("list", length = reps)
    tuningCVRes <- vector("list", length = reps)
    tuningCHullVAFRes <- vector("list", length = reps)
    tuningCHullCVRes <- vector("list", length = reps)
    estimationISRes <- vector("list", length = reps)
    estimationBICRes <- vector("list", length = reps)
    estimationCHullVAFRes <- vector("list", length = reps)
    estimationCHullCVRes <- vector("list", length = reps)
    estimationCVRes <- vector("list", length = reps)
    estimationCV1stErrorRuleRes <- vector("list", length = reps)


    for (i in 1:reps) {
        try({
            comdisspar <- sparsify(comdis, sparsity)
            variances <- makeVariance(varianceOfComps, p = p, error = error)
            datObject <- makeDat(n = n, p = p, ncomp = ncomp, comdisspar, variances = variances)
            Wstart <- varimax(svd(datObject$X)$v[, 1:ncomp])$loadings[, 1:ncomp]

            tuningISandBIC <- indexOfSparseness(dat = datObject$X,
                              groups = groups,
                              ncomp = ncomp,
                              ridge = ridgeSeq,
                              lassoSeq = lassoSeq,
                              grouplassoSeq = grouplassoSeq,
                              elitistlassoSeq = elitistlassoSeq,
                              Wstart = Wstart,
                              nrStarts = 1)

            tuningCV <- CVEigenVector(dat = datObject$X,
                              groups = groups,
                              ncomp = ncomp,
                              positions = NULL,
                              ridgeSeq = ridgeSeq,
                              lassoSeq = lassoSeq,
                              grouplassoSeq = grouplassoSeq,
                              elitistlassoSeq = elitistlassoSeq,
                              nrStarts = 1,
                              Wstart = Wstart)

            tuningCHullVAF <- myCHull(complexity = tuningISandBIC$nNonZeroCoef,
                              criterion = tuningISandBIC$VAF, upperLower = "upper",
                              combs = tuningISandBIC$combs)
            
            tuningCHullCV <- myCHull(complexity = tuningCV$nNonZeroCoef,
                              criterion = tuningCV$MSE, upperLower = "lower",
                              combs = tuningCV$combs)

            datObjectRes[[i]] <- datObject
            tuningISandBICRes[[i]] <- tuningISandBIC
            tuningCVRes[[i]] <- tuningCV
            tuningCHullVAFRes[[i]] <- tuningCHullVAF 
            tuningCHullCVRes[[i]] <- tuningCHullCV 

            estimationISRes[[i]] <- determineQuality(datObject = datObject, comdis = comdis,
                              comdisspar = comdisspar, p = p, ncomp = ncomp, Wstart = Wstart, 
                              groups = groups, tuningSettings = tuningISandBIC$bestTuningSeqsIS)

            estimationBICRes[[i]] <- determineQuality(datObject = datObject, comdis = comdis, 
                              comdisspar = comdisspar, p = p, ncomp = ncomp, Wstart = Wstart,
                              groups = groups, tuningSettings = tuningISandBIC$bestTuningSeqsBIC)

            estimationCHullVAFRes[[i]] <- determineQuality(datObject = datObject, comdis = comdis, 
                              comdisspar = comdisspar, p = p, ncomp = ncomp, Wstart = Wstart, 
                              groups = groups, tuningSettings = tuningCHullVAF$bestTuningSeqsCHull)

            estimationCHullCVRes[[i]] <- determineQuality(datObject = datObject, comdis = comdis,
                              comdisspar = comdisspar, p = p, ncomp = ncomp, 
                              Wstart = Wstart, groups = groups, 
                              tuningSettings = tuningCHullCV$bestTuningSeqsCHull)

            estimationCVRes[[i]] <- determineQuality(datObject = datObject, comdis = comdis, 
                              comdisspar = comdisspar, p = p, ncomp = ncomp, Wstart = Wstart, 
                              groups = groups, tuningSettings = tuningCV$bestTuningSeqsCV)

            estimationCV1stErrorRuleRes[[i]] <- determineQuality(datObject = datObject, 
                             comdis = comdis, comdisspar = comdisspar, p = p, ncomp = ncomp, 
                             Wstart = Wstart, groups = groups, 
                             tuningSettings = tuningCV$bestTuningSeqsCV1stdErrorRule)
        #},
        #    error = function(e){
        #        message(e)
        })

    }
    return(list(datObjectRes = datObjectRes,  
                tuningISandBICRes = tuningISandBICRes,  
                tuningCVRes = tuningCVRes,  
                tuningCHullVAFRes = tuningCHullVAFRes,  
                estimationISRes = estimationISRes,
                estimationBICRes = estimationBICRes ,
                estimationCHullVAFRes = estimationCHullVAFRes,
                estimationCHullCVRes = estimationCHullCVRes ,
                estimationCVRes = estimationCVRes,
                estimationCV1stErrorRuleRes = estimationCV1stErrorRuleRes))
}


##############
# set the conditions

nCond <- list(25, 50, 100)
pCond <- list(50)
ncompCond <- list(3)
sparsityCond <- list(c(0.3, 0.3, 0.3), c(0.8,0.8,0.8))
variancesCond  <- list(c(30, 30, 30))
errorCond <- list(0.05, 0.2)

ridgeSeq <- exp(seq(log(0.00001), log(500), length.out = 10))     
lassoSeq <- exp(seq(log(0.00001), log(500), length.out = 50))     
grouplassoSeq <- 0
elitistlassoSeq <- 0

conditions <- list(nCond, pCond, ncompCond, sparsityCond, variancesCond, errorCond)

conditionCombs <- expand.grid(1:length(nCond), 1:length(pCond), 1:length(ncompCond),
                             1:length(sparsityCond), 1:length(variancesCond), 1:length(errorCond))

nrow(conditionCombs)
conditionCombs

#create duplicates of the conditions so you can run them multiple times
#so the cores get used more optimally
reps <- 1 
x <- 1 
#in total reps*x per conditions
conditionCombsReplicated <- conditionCombs
for (i in 1:(x-1)) {
    conditionCombsReplicated <- rbind(conditionCombsReplicated, conditionCombs)
}

conditionCombsReplicated
(1:nrow(conditionCombsReplicated) %% nrow(conditionCombs)) #these groups belong together


#function to distribute the conditions to the simulation study
distribute_conditions <- function(i, reps, conditionCombs, conditions,
                                  ridgeSeq, lassoSeq, grouplassoSeq, elitistlassoSeq){

    set.seed(i)
    n <- conditions[[1]][[conditionCombs[i, 1]]]
    p <- conditions[[2]][[conditionCombs[i, 2]]]
    ncomp <- conditions[[3]][[conditionCombs[i, 3]]]
    sparsity <- conditions[[4]][[conditionCombs[i, 4]]]
    varianceOfComps <- conditions[[5]][[conditionCombs[i, 5]]]
    error <- conditions[[6]][[conditionCombs[i, 6]]]

    groups <- rep(p/2, 2)
    comdis <- matrix(1, p, ncomp)
    #comdis[1:(p/2), 1]  <- 0
    #comdis[(p/2+1):p, 2]  <- 0

    out <- simulation_study(reps = reps, n = n, p = p, ncomp = ncomp, groups = groups, 
                            comdis = comdis, sparsity = sparsity, 
                            varianceOfComps = varianceOfComps, error = error,
                            ridgeSeq = ridgeSeq, lassoSeq = lassoSeq, grouplassoSeq = grouplassoSeq,
                            elitistlassoSeq = elitistlassoSeq)

    saveRDS(out, file=paste("./simulation_study_modelselection_results/model_selection_result_", i, sep = ""))
    return(out)
}


########### Start the simulation study
# intialise cores
cores <- detectCores() - 1
cores
cl <- makeCluster(cores, outfile = './simstudy_log')
registerDoParallel(cl)

rm(simulation_study_results)

#Start measuring time
startTime <- Sys.time()

#run the simulation study
simulation_study_results <- foreach(i = 1:nrow(conditionCombsReplicated),
                                    .packages = c('testPackage', 'multichull', 'combinat'),
                                    .errorhandling = 'pass') %dopar% {

    source("./modelSelection.R", local = TRUE)
    source("./tucker.R", local = TRUE)
    source("./makeData.R", local = TRUE)

    distribute_conditions(i, reps = reps, conditionCombsReplicated[1:3, ], conditions, ridgeSeq, lassoSeq,
                            grouplassoSeq, elitistlassoSeq)
}

simulation_study_results

stopCluster(cl)
endTime <- Sys.time()
timeTaken  <- endTime - startTime
timeTaken



############################################################################################
############################################################################################
############################################################################################
############################################################################################
############################################################################################
############################################################################################
############################################################################################
############################################################################################
############################################################################################
############################################################################################
############################################################################################
############################################################################################
############################################################################################
############################################################################################
############################################################################################
############################################################################################
############################################################################################
############################################################################################
############################################################################################

############################################################################################
# Simulation study: multi-set!


##############
# set the conditions

nCond <- list(25, 100)
pCond <- list(50)
ncompCond <- list(3)
sparsityCond <- list(c(0.3, 0.3, 0.3), c(0.8,0.8,0.8))
variancesCond  <- list(c(30, 30, 30), c(30, 30, 15))
errorCond <- list(0.05, 0.2)

ridgeSeq <- 0
lassoSeq <- exp(seq(log(0.00001), log(500), length.out = 10))     
grouplassoSeq <- exp(seq(log(0.00001), log(500), length.out = 10))      
elitistlassoSeq <- 0

conditions <- list(nCond, pCond, ncompCond, sparsityCond, variancesCond, errorCond)

conditionCombs <- expand.grid(1:length(nCond), 1:length(pCond), 1:length(ncompCond),
                             1:length(sparsityCond), 1:length(variancesCond), 1:length(errorCond))


nrow(conditionCombs)
conditionCombs

#create duplicates of the conditions so you can run them multiple times
#so the cores get used more optimally
reps <- 5 
x <- 10 
#in total reps*x per conditions
conditionCombsReplicated <- conditionCombs
for (i in 1:(x-1)) {
    conditionCombsReplicated <- rbind(conditionCombsReplicated, conditionCombs)
}

conditionCombsReplicated
(1:nrow(conditionCombsReplicated) %% nrow(conditionCombs)) #these groups belong together


#function to distribute the conditions to the simulation study
distribute_conditions <- function(i, reps, conditionCombs, conditions,
                                  ridgeSeq, lassoSeq, grouplassoSeq, elitistlassoSeq){

    set.seed(i)
    n <- conditions[[1]][[conditionCombs[i, 1]]]
    p <- conditions[[2]][[conditionCombs[i, 2]]]
    ncomp <- conditions[[3]][[conditionCombs[i, 3]]]
    sparsity <- conditions[[4]][[conditionCombs[i, 4]]]
    varianceOfComps <- conditions[[5]][[conditionCombs[i, 5]]]
    error <- conditions[[6]][[conditionCombs[i, 6]]]

    groups <- rep(p/2, 2)
    comdis <- matrix(1, p, ncomp)
    comdis[1:(p/2), 1]  <- 0
    comdis[(p/2+1):p, 2]  <- 0

    out <- simulation_study(reps = reps, n = n, p = p, ncomp = ncomp, groups = groups, 
                            comdis = comdis, sparsity = sparsity, 
                            varianceOfComps = varianceOfComps, error = error,
                            ridgeSeq = ridgeSeq, lassoSeq = lassoSeq, grouplassoSeq = grouplassoSeq,
                            elitistlassoSeq = elitistlassoSeq)

    saveRDS(out, file=paste("./simulation_study_modelselection_multiset_results/model_selection_result_", i, sep = ""))
    return(out)
}


########### Start the simulation study
# intialise cores
cores <- detectCores() - 1
cores
cl <- makeCluster(cores, outfile = './simstudy_log')
registerDoParallel(cl)

rm(simulation_study_results)

#Start measuring time
startTime <- Sys.time()

#run the simulation study
simulation_study_results <- foreach(i = 1:nrow(conditionCombsReplicated),
                                    .packages = c('testPackage', 'multichull'),
                                    .errorhandling = 'pass') %dopar% {

    source("./modelSelection.R", local = TRUE)
    source("./tucker.R", local = TRUE)
    source("./makeData.R", local = TRUE)

    distribute_conditions(i, reps = reps, conditionCombsReplicated, conditions, ridgeSeq, lassoSeq,
                            grouplassoSeq, elitistlassoSeq)
}

stopCluster(cl)
endTime <- Sys.time()
timeTaken  <- endTime - startTime
timeTaken


