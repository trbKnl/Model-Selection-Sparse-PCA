################################################################################
# ANALYSIS HERRING DATA
#

rm(list=ls())

library(ggplot2)
library(Rcpp)
source("../modelSelection.R")
sourceCpp("../mmsca.cpp")


load("./Herring.rda")
class(Herring)
names(Herring)
summary(Herring$Herring_ChemPhy)
summary(Herring$Herring_Sensory)

Herring$Herring_Sensory

##################################

Her_ChemPhy <- scale(Herring$Herring_ChemPhy)  
Her_Sensory <- scale(Herring$Herring_Sensory)

Her_ChemPhy  <- Her_ChemPhy / sum(Her_ChemPhy^2)
Her_Sensory  <- Her_Sensory / sum(Her_Sensory^2)

Her_Cat <- cbind(Her_ChemPhy, Her_Sensory)

#values are now very small
min(Her_Cat)
max(Her_Cat)

#scale by a constant
Her_Cat <- Her_Cat*10


## Block 1: ChemPhy
##Inspect singular values to get an idea of the number of components
# Elbow rule: 4 Components
plot(svd(Her_ChemPhy)$d^2)

out <- vector("list", 8)
for (i in 1:8) {
    Q <- i
    out[[i]] <- EigenVectorCV2(dat = Her_ChemPhy ,
                   groups = c(ncol(Her_ChemPhy)),
                   ncomp = Q,
                   ridge = rep(0, Q),
                   lasso = rep(0, Q),
                   grouplasso = rep(0, Q),
                   elitistlasso = rep(0, Q),
                   nrFolds = 10,
                   Wstart = matrix(0, ncol(Her_ChemPhy), Q),
                   nScale = 0,
                   nrStarts = 1,
                   scaleDat = FALSE)
}

MSE <- unlist(lapply(out, FUN = function(x) { x$MSE }))
stdError <- unlist(lapply(out, FUN = function(x) { x$stdError }))
df <- data.frame(MSE, stdError)

#CV says 3 component
ggplot(df, aes(x=1:nrow(df), y=MSE)) +
    geom_point() +
    geom_errorbar(aes(ymin=MSE-stdError, ymax=MSE+stdError)) +
    theme_bw()


## Block 1: Sensory
##Inspect singular values to get an idea of the number of components
# Elbow rule: 3 Components
plot(svd(Her_Sensory)$d^2)

out <- vector("list", 8)
for (i in 1:8) {
    Q <- i
    out[[i]] <- EigenVectorCV2(dat = Her_Sensory,
                   groups = c(ncol(Her_Sensory)),
                   ncomp = Q,
                   ridge = rep(0, Q),
                   lasso = rep(0, Q),
                   grouplasso = rep(0, Q),
                   elitistlasso = rep(0, Q),
                   nrFolds = 10,
                   Wstart = matrix(0, ncol(Her_Sensory), Q),
                   nScale = 0,
                   nrStarts = 1,
                   scaleDat = FALSE)
}

MSE <- unlist(lapply(out, FUN = function(x) { x$MSE }))
stdError <- unlist(lapply(out, FUN = function(x) { x$stdError }))
df <- data.frame(MSE, stdError)

#CV says 3 component
ggplot(df, aes(x=1:nrow(df), y=MSE)) +
    geom_point() +
    geom_errorbar(aes(ymin=MSE-stdError, ymax=MSE+stdError)) +
    theme_bw()



## concatenated Data
##Inspect singular values to get an idea of the number of components
# Elbow rule: 6 Components
plot(svd(Her_Cat)$d^2)

out <- vector("list", 8)
for (i in 1:8) {
    Q <- i
    out[[i]] <- EigenVectorCV2(dat = Her_Cat,
                   groups = c(ncol(Her_Cat)),
                   ncomp = Q,
                   ridge = rep(0, Q),
                   lasso = rep(0, Q),
                   grouplasso = rep(0, Q),
                   elitistlasso = rep(0, Q),
                   nrFolds = 10,
                   Wstart = matrix(0, ncol(Her_Cat), Q),
                   nScale = 0,
                   nrStarts = 1,
                   scaleDat = FALSE)
}

MSE <- unlist(lapply(out, FUN = function(x) { x$MSE }))
stdError <- unlist(lapply(out, FUN = function(x) { x$stdError }))
df <- data.frame(MSE, stdError)

#CV says 6 component
#CV 1std Err rule 4 Components
ggplot(df, aes(x=1:nrow(df), y=MSE)) +
    geom_point() +
    geom_errorbar(aes(ymin=MSE-stdError, ymax=MSE+stdError)) +
    theme_bw()

###########################################################################
# Tuning of the parameters
###########################################################################
# We take 6 components 3 for the first block 3 for the second block


ridgeSeq <- exp(seq(log(0.00001), log(10), length.out = 25))     
lassoSeq <- exp(seq(log(0.00001), log(500), length.out = 25))     
grouplassoSeq <- exp(seq(log(0.00001), log(500), length.out = 25))   
elitistlassoSeq <- 0
Q <- 6 
ncol(Her_ChemPhy)
ncol(Her_Sensory)

groups <- c(ncol(Her_ChemPhy), ncol(Her_Sensory))

combs <- combOfTuningParams(ridgeSeq, lassoSeq, grouplassoSeq, elitistlassoSeq, ncomp = Q)

nmodels <- nrow(combs[[1]])

out <- vector("list", nmodels)
for (i in 1:nmodels) {
    print(i)
    print(combs[[1]][i, ])
    print(combs[[2]][i, ])
    print(combs[[3]][i, ])
    print(combs[[4]][i, ])
    out[[i]] <- EigenVectorCV2(dat = Her_Cat ,
               groups = groups,
               ncomp = Q,
               ridge = combs[[1]][i, ], 
               lasso = combs[[2]][i, ],
               grouplasso = combs[[3]][i, ],
               elitistlasso = combs[[4]][i, ],
               nrFolds = 10,
               Wstart = matrix(0, ncol(Her_Cat), Q),
               nScale = 0,
               nrStarts = 1,
               scaleDat = FALSE)
}

best <- which.min(MSE)
candidates <- stdError[best] + MSE[best] > MSE


candidates
models <- 1:nmodels
candidates_models <- models[candidates]
length(candidates_models)


#Select the most sparse model based given the elegible candidate models 

ncandidates_models <- length(candidates_models)
sparsity <- rep(NA, ncandidates_models)

for (i in 1:ncandidates_models) {
    results <- newAlgoCpp(X = Her_Cat,
               ridge = combs[[1]][candidates_models[i], ],
               lasso = combs[[2]][candidates_models[i], ],
               constraints = matrix(1, ncol(Her_Cat), Q),
               grouplasso = combs[[3]][candidates_models[i], ],
               elitistlasso = combs[[4]][candidates_models[i], ],
               groups = groups,
               Q = Q,
               itr = 100000,
               Wstart = matrix(0, ncol(Her_Cat), Q),
               nStarts = 1,
               printLoss = FALSE,
               coorDec = TRUE)
    sparsity[i] <- sum(results$W == 0)
}

mostsparse <- which.max(sparsity)
bestmodel <- candidates_models[mostsparse]

saveRDS(bestmodel, "./bestmodelnumber.rds")
bestmodel <- readRDS("./bestmodelnumber.rds")

combs[[1]][bestmodel, ]
combs[[2]][bestmodel, ]
combs[[3]][bestmodel, ]
combs[[4]][bestmodel, ]

results <- newAlgoCpp(X = Her_Cat,
       ridge = combs[[1]][bestmodel, ],
       lasso = combs[[2]][bestmodel, ],
       constraints = matrix(1, ncol(Her_Cat), Q),
       grouplasso = combs[[3]][bestmodel, ],
       elitistlasso =combs[[4]][bestmodel, ],
       groups = groups,
       Q = Q,
       itr = 100000,
       Wstart = matrix(0, ncol(Her_Cat), Q),
       nStarts = 1,
       printLoss = FALSE,
       coorDec = TRUE)


rownames(results$W) <-  c(colnames(Herring$Herring_ChemPhy), colnames(Herring$Herring_Sensory))

library(xtable)

results$W

xtable(results$W, digits = 3)



