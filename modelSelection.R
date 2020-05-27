###############################################################
# Function to perform model selection for sparse SCA
###############################################################

indexOfSparseness <- function(dat, groups, ncomp, positions = NULL,
                              ridgeSeq, lassoSeq, grouplassoSeq,
                              elitistlassoSeq, nrStarts, Wstart){

    I <- nrow(dat)
    p <- ncol(dat)
    total_coef <- sum(groups) * ncomp

    # Function to perform model selection for all combinations of 
    # Given tuning parameters
    combs <- combOfTuningParams(ridgeSeq, lassoSeq, grouplassoSeq,
                       elitistlassoSeq, positions = positions, ncomp = ncomp)

    VarSelect0 <- newAlgoCpp(dat,
                       ridge = rep(0, ncomp),
                       lasso = rep(0, ncomp),
                       #constraints = comdis,
                       constraints = matrix(1, p, ncomp),
                       grouplasso = rep(0, ncomp),
                       elitistlasso = rep(0, ncomp),
                       groups = groups,
                       Q = ncomp,
                       itr = 100000,
                       nStarts = nrStarts,
                       printLoss = FALSE,
                       Wstart = matrix(0,1,1),
                       coorDec = FALSE)

    T_hat0 <- dat %*% VarSelect0$W
    P_hat0 <- VarSelect0$P

    V_s <- sum((T_hat0%*%t(P_hat0))^2) # this is for IS
    V_oo <- sum(dat^2)  # this is for Index of sparseness (IS)
    residual_variance0 <- sum((dat - T_hat0 %*% t(P_hat0))^2)

    IS <- rep(NA, nrow(combs[[1]]))
    VAF <- rep(NA, nrow(combs[[1]]))
    BIC <- rep(NA, nrow(combs[[1]]))
    nNonZeroCoef <- rep(NA, nrow(combs[[1]]))
    ratio_zerocoef <- rep(NA, nrow(combs[[1]]))
    initial_VAF<- rep(NA, nrow(combs[[1]]))

    for(i in 1:nrow(combs[[1]])){
    
        VarSelect <- newAlgoCpp(dat,
                           ridge = combs[[1]][i, ],
                           lasso = combs[[2]][i, ],
                           #constraints = comdis,
                           constraints = matrix(1, p, ncomp),
                           grouplasso = combs[[3]][i, ],
                           elitistlasso = combs[[4]][i, ],
                           groups = groups,
                           Q = ncomp,
                           itr = 100000,
                           nStarts = nrStarts,
                           printLoss = FALSE,
                           Wstart = Wstart,
                           coorDec = FALSE)

        T_hat <- dat %*% VarSelect$W
        P_hat <- VarSelect$P
        W_hat <- VarSelect$W
        nonZeroCoefInW_hat <- sum(W_hat != 0)

        residual_variance <- sum((dat - T_hat %*% t(P_hat))^2)
        V_a <- sum((T_hat %*% t(P_hat))^2)  
        nNonZeroCoef[i] <- nonZeroCoefInW_hat

        IS[i] <- V_a * V_s / V_oo^2 * (sum(W_hat == 0) / total_coef)
        BIC[i]  <- (residual_variance / residual_variance0) + (nonZeroCoefInW_hat * (log(I) / I))

        VAF[i] <- V_a / V_oo
        initial_VAF[i] <- V_s / V_oo
        ratio_zerocoef[i] <- sum(W_hat == 0) / total_coef 
    }
    bestTuningSeqsIS <- lapply(combs, function(x){return(x[which.max(IS), ])})
    bestTuningSeqsBIC <- lapply(combs, function(x){return(x[which.min(BIC), ])})
    return(list(combs = combs, bestTuningSeqsIS = bestTuningSeqsIS,
                bestTuningSeqsBIC = bestTuningSeqsBIC,
                IS = IS, VAF = VAF, BIC = BIC, initial_VAF = initial_VAF, 
                ratio_zerocoef = ratio_zerocoef, nNonZeroCoef = nNonZeroCoef))
}


denester <- function(listobj){
    #Function to de-nest a nested list:
    #if the list contains lists that are not nested
  out <- list()
    for(i in 1:length(listobj)){
        if(is.list(listobj[[i]])){
            out <- c(out, listobj[[i]])
        } else {
            out <- c(out, list(listobj[[i]]))
        }
    }
    return(out)
}

combOfTuningParams <- function(ridgeSeq, lassoSeq, grouplassoSeq, elitistlassoSeq, positions = NULL, ncomp = NULL){
    # A not so elegant function to get all combinations of tuning  parameters.
    # But very flexible and useful,
    # it does not matter that it is horrible, it is only used once
    # per model selection sequence.

    lambdaList <- denester(list(ridgeSeq, lassoSeq, grouplassoSeq, elitistlassoSeq))

    totalNumberOfCombs <- prod(unlist(lapply(lambdaList, length)))
    if(totalNumberOfCombs > 50000){
        stop(paste("Not gonna happen!, total n comb :", totalNumberOfCombs, sep=""))
        return(NA) 
    }

    print(paste("Total number of models: ", totalNumberOfCombs, sep=""))

    if(is.null(positions)){
        if(length(lambdaList) != 4){stop("Given sequences != 4")}
        if(is.null(ncomp)){stop("Given ncomp")}
        combmat <- as.list(expand.grid(lambdaList))
        names(combmat) <- NULL
        outlist <- lapply(combmat, function(x){out <- matrix(NA, length(x), ncomp); out[] <- x; return(out)})
       
        return(outlist)
    }

    nseq <- sum(apply(positions, 1, function(x){length(unique(x))}))
    if(nseq != length(lambdaList)){
        stop(paste("Number of given sequences: ", length(lambdaList),  ". Does not match, the number of sequences specified in positions: ", nseq, sep = ""))
    }

    comblist <- list()
    counter <- 0

    for(j in 1:nrow(positions)){
        for(i in 1:length(unique(positions[j, ]))){
            counter <- counter + 1
            comblist[[length(comblist) + 1]]  <- lambdaList[[counter]]
        }
    }
    combmat <- as.matrix(expand.grid(comblist))
    counter <- 1
    outlist <- list()
    for(j in 1:nrow(positions)){
        out <- matrix(NA, nrow(combmat), ncol(positions))
        identifiers <- unique(positions[j, ])
        for(i in 1:length(unique(positions[j, ]))){
            out[, which(identifiers[i] == positions[j, ])] <- combmat[, counter ] 
            counter <- counter + 1
        }
        outlist[[j]] <- out
    }
    return(outlist)
}


# A function that searches the grid of possible tuning paramters
# It creates a hypercube around the space of interest
# looks for where the best combination of tuning parameters is 
# according to IS and makes a hypercube around that point with half the previous range
# algorithm stops when the range is small, check the function for how small exactly
ISexperimental <- function(dat, groups, ncomp, 
                              ridgeSeq, lassoSeq, grouplassoSeq,
                              elitistlassoSeq, nrStarts, Wstart, maxitr = 10000, 
                              positions = NULL, stepsize, logscale = FALSE){

    p  <- ncol(dat)
    seql <- list(ridgeSeq, lassoSeq, grouplassoSeq, elitistlassoSeq)
    VarSelect0 <- newAlgoCpp(dat,
                       ridge = rep(0, ncomp),
                       lasso = rep(0, ncomp),
                       #constraints = comdis,
                       constraints = matrix(1, p, ncomp),
                       grouplasso = rep(0, ncomp),
                       elitistlasso = rep(0, ncomp),
                       groups = groups,
                       Q = ncomp,
                       itr = 100000,
                       nStarts = nrStarts,
                       printLoss = FALSE,
                       Wstart = matrix(0,1,1),
                       coorDec = FALSE)
    T_hat0 <- dat %*% VarSelect0$W
    P_hat0 <- VarSelect0$P
    V_s <- sum((T_hat0%*%t(P_hat0))^2) # this is for IS
    V_oo <- sum(dat^2)  # this is for Index of sparseness (IS)


    # The function is written like this because you can specify
    # different tuning sequences per pernalty per component.
    # In order for the function to be usable 
    range  <- rapply(seql, function(x){return((max(x)-min(x))/2)}, how="list") 
    best <- rapply(seql, mean, how="list") 
    counter <- 0

    while(mean(unlist(range)) > 0.05){
        for(i in 1:length(seql)){
            if(is.list(seql[[i]])){
                for(j in 1:length(seql[[i]])){
                    seql[[i]][[j]] <- determineRange(best[[i]][[j]], range[[i]][[j]], 
                                                     stepsize, logscale = logscale)
                }
            } else {
                    seql[[i]] <- determineRange(best[[i]], range[[i]], 
                                                stepsize, logscale = logscale)
            }
        }

        range  <- rapply(seql, function(x){return((max(x)-min(x))/4)}, how ="list") 

        combs <- combOfTuningParams(seql[[1]], seql[[2]],
                                    seql[[3]], seql[[4]],
                                    ncomp = ncomp, positions = positions)

        IS <- rep(NA, nrow(combs[[1]]))
        for(i in 1:nrow(combs[[1]])){

            VarSelect <- newAlgoCpp(dat,
                               ridge = combs[[1]][i, ],
                               lasso = combs[[2]][i, ],
                               #constraints = comdis,
                               constraints = matrix(1, p, ncomp),
                               grouplasso = combs[[3]][i, ],
                               elitistlasso = combs[[4]][i, ],
                               groups = groups,
                               Q = ncomp,
                               itr = maxitr,
                               nStarts = nrStarts,
                               printLoss = FALSE,
                               Wstart = Wstart,
                               coorDec = FALSE)

            T_hat <- dat %*% VarSelect$W
            P_hat <- VarSelect$P
            W_hat <- VarSelect$W

            V_a <- sum((T_hat %*% t(P_hat))^2) 
            IS[i] <- V_a * V_s / V_oo^2 * sum(W_hat == 0) / (sum(groups) * ncomp) 
        }

        if(is.null(positions)){
            best <- lapply(combs, function(x){return(x[which.max(IS), 1])})
            if(counter %% 4 == 0){
                oldbest <- best
            } else if (counter %% 4 == 3) {
                for(i in 1:length(best)){
                    if(best[[i]] == oldbest[[i]]){
                        range[[i]]  <- 0
                    }
                }
            }
            counter <- counter + 1
        } else {
            #re-distribute the correct best values over the nested list "Best"
            newBest <- lapply(combs, function(x){return(x[which.max(IS), ])})
            for(i in 1:4){
                identifiers  <- unique(positions[i, ])
                for(j in 1:length(identifiers)){
                    where <- which(identifiers[j] == positions[i, ])[1]
                    best[[i]][[j]] <- newBest[[i]][where]
                }
            }
        }
    }

    bestTuningSeqs <- lapply(combs, function(x){return(x[which.max(IS), ])})
    return(list(bestTuningSeqs=bestTuningSeqs))
}

determineRange <- function(best, range, stepsize, logscale){
    if(range != 0){
        minrange <- ifelse(best - range > 0, 
                           best - range, 0)
        maxrange <- ifelse(best + range > 0,
                           best + range, 0)
        if (logscale == TRUE) {
            return(exp(seq(log(minrange + 0.000001), log(maxrange), length.out = stepsize)))
        } else {
            return(seq(minrange, maxrange, length=stepsize))
        }
    } else {
        return(best)
    }
}





# Eigenvector crossvalidation method
EigenVectorCV2 <- function(dat, groups, ncomp, ridge, lasso, grouplasso, elitistlasso,
                           nrFolds, Wstart, nScale = 0, nrStarts, scaleDat = FALSE){ 

    folds <- rep_len(1:nrFolds, nrow(dat))
    #folds <- sample(folds)
    cvError  <- matrix(NA, nrow(dat), ncol(dat))
    MSEkthFold <- rep(NA, nrFolds)
    p <- ncol(dat)

    for(a in 1:nrFolds){

        # actual split of the data
        fold <- which(folds == a)
        trainDat <- dat[-fold, ]
        testDat <- dat[fold, , drop = FALSE]

        if(scaleDat){
            #Scale training and the testing data
            means <- apply( trainDat, 2, mean )		
            trainDat <- t( t(trainDat) - means)
            sdTrain <- apply(trainDat, 2, function(x) sqrt(sum( x^2, na.rm = TRUE ) / (length(na.omit(x)) - nScale )))
            trainDat <- t( t(trainDat)/sdTrain )
            
            testDat <-  t( t(testDat) - means)
            testDat <- t( t(testDat)/sdTrain ) 
            
            rm(means, sdTrain)		
        } 

        #model estimation
        res <- newAlgoCpp(dat,
                   ridge = ridge,
                   lasso = lasso ,
                   #constraints = comdis,
                   constraints = matrix(1, p, ncomp),
                   grouplasso = grouplasso,
                   elitistlasso = elitistlasso,
                   groups = groups,
                   Q = ncomp,
                   itr = 100000,
                   nStarts = 1,
                   printLoss = FALSE,
                   Wstart = Wstart,
                   coorDec = TRUE)

        #Eigenvector crossvalidation Bro Kiers
        pred <- matrix(NA, nrow(testDat), ncol(testDat))
        print(a)
        for(i in 1:ncol(dat)){
            TMinusVariableJ <- testDat[, -i] %*% res$W[-i, ]
            pred[, i] <- TMinusVariableJ %*% res$P[i, ] 
            }
        cvError[fold, ] <- (testDat - pred)^2
        MSEkthFold[a]  <-  mean(cvError[fold, ]) 
        }

    return(list(cvError = cvError, MSE = mean(MSEkthFold), MSEkthFold = MSEkthFold, 
                stdError = sd(MSEkthFold) / sqrt(nrFolds)))
}



CVEigenVector <- function(dat, groups, ncomp, positions = NULL,
                              ridgeSeq, lassoSeq, grouplassoSeq,
                              elitistlassoSeq, nrStarts, Wstart){

    I <- nrow(dat)
    p <- ncol(dat)
    total_coef <- sum(groups) * ncomp

    # Function to perform model selection for all combinations of 
    # Given tuning parameters
    combs <- combOfTuningParams(ridgeSeq, lassoSeq, grouplassoSeq,
                       elitistlassoSeq, positions = positions, ncomp = ncomp)

    MSE <- rep(NA, nrow(combs[[1]]))
    MSEstdError <- rep(NA, nrow(combs[[1]]))
    nNonZeroCoef <- rep(NA, nrow(combs[[1]]))

    for(i in 1:nrow(combs[[1]])){

        resCV <- EigenVectorCV2(dat,
                       groups = groups,
                       ncomp = ncomp,
                       ridge = combs[[1]][i, ], 
                       lasso = combs[[2]][i, ], 
                       grouplasso = combs[[3]][i, ], 
                       elitistlasso = combs[[4]][i, ], 
                       nrFolds = 10,
                       Wstart = Wstart,
                       nScale = 0,
                       nrStarts = 1,
                       scaleDat = FALSE)

        resforComplexity <- newAlgoCpp(dat,
                       ridge = combs[[1]][i, ], 
                       lasso = combs[[2]][i, ],  
                       #constraints = comdis,
                       constraints = matrix(1, p, ncomp),
                       grouplasso = combs[[3]][i, ],  
                       elitistlasso =combs[[4]][i, ],  
                       groups = groups,
                       Q = ncomp,
                       itr = 100000,
                       nStarts = 1,
                       printLoss = FALSE,
                       Wstart = Wstart,
                       coorDec = FALSE)

        MSE[i] <- resCV$MSE
        MSEstdError[i] <- resCV$stdError 
        nNonZeroCoef[i] <- sum(resforComplexity$W != 0)
    }

    bestModel <- which.min(MSE)
    MSE1stdErrorRule <- MSE[MSE < (MSE[bestModel] + MSEstdError[bestModel])]
    bestModel1stdErrorRule <- which(MSE1stdErrorRule[length(MSE1stdErrorRule)] == MSE)

    bestTuningSeqsCV <- lapply(combs, function(x){return(x[which.min(MSE), ])})
    bestTuningSeqsCV1stdErrorRule <- lapply(combs, function(x){return(x[bestModel1stdErrorRule, ])})
    return(list(combs = combs, bestTuningSeqsCV = bestTuningSeqsCV,
                MSE = MSE, bestModel1stdErrorRule = bestModel1stdErrorRule,
                bestTuningSeqsCV1stdErrorRule = bestTuningSeqsCV1stdErrorRule, 
                MSEstdError = MSEstdError, nNonZeroCoef = nNonZeroCoef ))
}



#wrapper around the CHull function from the multi Chull package
myCHull <- function(complexity, criterion, upperLower, combs) {
    bestTuningSeqsCHull <- list()
    bestModelCHull <- NA

    try({
    CHullData <- data.frame(complexity, criterion)
    CHullRes <- CHull(CHullData, bound = upperLower)

    #sometimes there are more than 1 model on the bounds
    #In that case select the model with the highest fit
    if (nrow(CHullRes$Solution) == 1){
        bestModelCHull <- as.numeric(regmatches(rownames(CHullRes$Solution), regexpr("[0-9]+", rownames(CHullRes$Solution))))
    } else {
        bestfit <- which.max(CHullRes$Solution$fit)
        bestModelCHull <- as.numeric(regmatches(rownames(CHullRes$Solution[bestfit, ]), regexpr("[0-9]+", rownames(CHullRes$Solution[bestfit, ]))))

    }

    bestTuningSeqsCHull <- lapply(combs, function(x){return(x[bestModelCHull, ])})

    })

    return(list(bestTuningSeqsCHull = bestTuningSeqsCHull, bestModelCHull = bestModelCHull))
}


