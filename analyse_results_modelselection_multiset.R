########### Make the multi set plots

rm(list=ls())
# script to make the plots
# This scripts makes the tables 1 and 2, Figure 4 (Tucker congruences) and 5 (SRD plot)
# 

source('./tucker.R')
##########
# These were the conditions used
# set the conditions

nCond <- list(25, 100)
pCond <- list(50)
ncompCond <- list(3)
sparsityCond <- list(c(0.3, 0.3, 0.3), c(0.8,0.8,0.8))
variancesCond  <- list(c(30, 30, 30))
errorCond <- list(0.05, 0.2)

ridgeSeq <- exp(seq(log(0.00001), log(500), length.out = 10))     
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
grouping <- (0:(nrow(conditionCombsReplicated)-1) %% nrow(conditionCombs)) + 1 #these groups belong together

#test <- readRDS("./simulation_study_modelselection_multiset_results/model_selection_result_1")
#names(test)
results_to_analyse <- vector("list", length = nrow(conditionCombs))

for (a in 1:nrow(conditionCombs)) {
    res <- matrix(NA, reps, 36)
    files <- 1:nrow(conditionCombsReplicated)  
    files <- files[grouping == a]
    resList <- vector("list", length = x)

    for (i in 1:x) {
        fileName <- paste("~/simulation_study_modelselection_multiset_results/model_selection_result_", files[i], sep = "")
        #fileName <- paste("./simulation_study_modelselection_multiset_results/model_selection_result_", files[i], sep = "")
        loadedFile <- readRDS(fileName)

        for (j in 1:reps) {
            res[j, 1] <- loadedFile$estimationISRes[[j]]$tuckerCongruenceResW
            res[j, 2] <- loadedFile$estimationBICRes[[j]]$tuckerCongruenceResW
            res[j, 3] <- loadedFile$estimationCVRes[[j]]$tuckerCongruenceResW
            res[j, 4] <- loadedFile$estimationCV1stErrorRuleRes[[j]]$tuckerCongruenceResW
            res[j, 5] <- loadedFile$estimationCHullCVRes[[j]]$tuckerCongruenceResW
            res[j, 6] <- loadedFile$estimationCHullVAFRes[[j]]$tuckerCongruenceResW

            res[j, 7] <- mean(loadedFile$estimationISRes[[j]]$corClassRes$zeroCoefFound)
            res[j, 8] <- mean(loadedFile$estimationBICRes[[j]]$corClassRes$zeroCoefFound)
            res[j, 9] <- mean(loadedFile$estimationCVRes[[j]]$corClassRes$zeroCoefFound)
            res[j, 10] <- mean(loadedFile$estimationCV1stErrorRuleRes[[j]]$corClassRes$zeroCoefFound)
            res[j, 11] <- mean(loadedFile$estimationCHullCVRes[[j]]$corClassRes$zeroCoefFound)
            res[j, 12] <- mean(loadedFile$estimationCHullVAFRes[[j]]$corClassRes$zeroCoefFound)

            res[j, 13] <- mean(loadedFile$estimationISRes[[j]]$corClassRes$nonZeroCoefFound)
            res[j, 14] <- mean(loadedFile$estimationBICRes[[j]]$corClassRes$nonZeroCoefFound)
            res[j, 15] <- mean(loadedFile$estimationCVRes[[j]]$corClassRes$nonZeroCoefFound)
            res[j, 16] <- mean(loadedFile$estimationCV1stErrorRuleRes[[j]]$corClassRes$nonZeroCoefFound)
            res[j, 17] <- mean(loadedFile$estimationCHullCVRes[[j]]$corClassRes$nonZeroCoefFound)
            res[j, 18] <- mean(loadedFile$estimationCHullVAFRes[[j]]$corClassRes$nonZeroCoefFound)

            res[j, 19] <- sum(loadedFile$estimationISRes[[j]]$comIdentifiedRes, na.rm = TRUE)
            res[j, 20] <- sum(loadedFile$estimationBICRes[[j]]$comIdentifiedRes, na.rm = TRUE)
            res[j, 21] <- sum(loadedFile$estimationCVRes[[j]]$comIdentifiedRes, na.rm = TRUE)
            res[j, 22] <- sum(loadedFile$estimationCV1stErrorRuleRes[[j]]$comIdentifiedRes, na.rm = TRUE)
            res[j, 23] <- sum(loadedFile$estimationCHullCVRes[[j]]$comIdentifiedRes, na.rm = TRUE)
            res[j, 24] <- sum(loadedFile$estimationCHullVAFRes[[j]]$comIdentifiedRes, na.rm = TRUE)

            res[j, 25] <- sum(loadedFile$estimationISRes[[j]]$disIdentifiedRes, na.rm = TRUE)
            res[j, 26] <- sum(loadedFile$estimationBICRes[[j]]$disIdentifiedRes, na.rm = TRUE)
            res[j, 27] <- sum(loadedFile$estimationCVRes[[j]]$disIdentifiedRes, na.rm = TRUE)
            res[j, 28] <- sum(loadedFile$estimationCV1stErrorRuleRes[[j]]$disIdentifiedRes, na.rm = TRUE)
            res[j, 29] <- sum(loadedFile$estimationCHullCVRes[[j]]$disIdentifiedRes, na.rm = TRUE)
            res[j, 30] <- sum(loadedFile$estimationCHullVAFRes[[j]]$disIdentifiedRes, na.rm = TRUE)

            res[j, 31] <- conditionCombs[a, 1]
            res[j, 32] <- conditionCombs[a, 2]
            res[j, 33] <- conditionCombs[a, 3]
            res[j, 34] <- conditionCombs[a, 4]
            res[j, 35] <- conditionCombs[a, 5]
            res[j, 36] <- conditionCombs[a, 6]

        }
        resList[[i]] <- res
    }
    resultsForCondition <-  do.call(rbind, resList)
    results_to_analyse[[a]] <- cbind(resultsForCondition, a)
}

results_to_analyse

results_to_analyse <- do.call(rbind, results_to_analyse)
results_to_analyse <- as.data.frame(results_to_analyse)

colnames(results_to_analyse)  <- c("tucker_IS", "tucker_BIC", "tucker_10foldCV", "tucker_10foldCV1stdError", "tucker_CHullMSE", "tucker_CHullVAF", "zeroCoefFound_IS", "zeroCoefFound_BIC", "zeroCoefFound_10foldCV",
"zeroCoefFound_10foldCV1stdError", "zeroCoefFound_CHullMSE", "zeroCoefFound_CHullVAF",
"nonZeroCoefFound_IS", "nonZeroCoefFound_BIC", "nonZeroCoefFound_10foldCV",
"nonZeroCoefFound_10foldCV1stdError", "nonZeroCoefFound_CHullMSE", "nonZeroCoefFound_CHullVAF", "commonIden_IS", "commonIden_BIC", "commonIden_10foldCV", "commonIden_10foldCV1stdError", "commonIden_CHullMSE", "commonIden_CHullVAF", "distinctiveIden_IS", "distinctiveIden_BIC", "distinctiveIden_10foldCV", "distinctiveIden_10foldCV1stdError", "distinctiveIden_CHullMSE", "distinctiveIden_CHullVAF", "n", "p", "ncomp", "sparsity", "variances", "error", "conditionNumber")

results_to_analyse <- reshape(results_to_analyse, varying = 1:30, sep = "_", timevar = "method", direction = "long")

results_to_analyse$conditionNumber <- as.factor(results_to_analyse$conditionNumber)
results_to_analyse$error <- as.factor(results_to_analyse$error)
results_to_analyse$n <- as.factor(results_to_analyse$n)
results_to_analyse$sparsity <- as.factor(results_to_analyse$sparsity)
results_to_analyse$variances <- as.factor(results_to_analyse$variances)


unique(results_to_analyse$method)
results_to_analyse$method[which(results_to_analyse$method == "10foldCV")] <- "10-fold CV"
results_to_analyse$method[which(results_to_analyse$method == "10foldCV1stdError")] <- "10-fold CV 1 SE rule"
results_to_analyse$method[which(results_to_analyse$method == "CHullMSE")] <- "CHull MSE"
results_to_analyse$method[which(results_to_analyse$method == "CHullVAF")] <- "CHull VAF"

levels(results_to_analyse$sparsity) <- c("Sparsity 30%", "Sparsity 80%")
#levels(results_to_analyse$variances) <- c("Variances 30 30 30", "Variances 30 30 15")
levels(results_to_analyse$variances) <- c("Variances 30 30 30")
levels(results_to_analyse$error) <- c("Error 5%", "Error 20%")
levels(results_to_analyse$n) <- c("25", "100")

library(ggplot2)
library(xtable)

##### Make a table of whether the common components are found
results_to_analyse1 <- results_to_analyse[results_to_analyse$commonIden == 1, ]
tab1 <- ((ftable(results_to_analyse1$method, results_to_analyse1$error,   results_to_analyse1$n,
                 results_to_analyse1$commonIden, results_to_analyse1$sparsity, row.vars = c(1, 4)) / 50) * 100)

tab1

xtableFtable(tab1)

results_to_analyse1 <- results_to_analyse
results_to_analyse1 <- results_to_analyse1[results_to_analyse1$distinctiveIden == 2, ]

##### Make a table of whether the distinctive components are found
tab2  <- (ftable(results_to_analyse1$method, results_to_analyse1$error,   results_to_analyse1$n, results_to_analyse1$distinctiveIden == 2, results_to_analyse1$sparsity, row.vars = c(1, 4)) / 50) * 100
tab2
xtableFtable(tab2)


########################### sum rank differences plot

# create scores 100% is the top scores, lower is better so 100 - scores
# The reference benchmark for each method is 0, which means 100% correct, methods start out equally good
scores <- 100 - cbind(rep(100, nrow(tab1)), as.matrix(tab1),  as.matrix(tab2))
scores
ranks <- apply(scores, 2, FUN = function(x) { rank(x, ties.method = "min")}  )
ranks
rankDiffs <- apply(ranks, 2, FUN = function(x) { abs(ranks[, 1] - x) } )
rankDiffs
SRDscores <- apply(ranks, 1, sum)
SRDscores

# approximate theoretical distribution of SRD scores in a simulation based on
# the rankings being equal chance
set.seed(1)
reps <- 10000
rankPermutations <- do.call(cbind, combinat::permn(6))
res <- matrix(NA, nrow(scores), reps) 
for (i in 1:reps) {
    ranksTheor <- cbind(rep(1, nrow(scores)),
                        rankPermutations[, sample(1:ncol(rankPermutations), ncol(scores) -1, replace = TRUE)])
    rankDiffsTheor <- apply(ranksTheor, 2, FUN = function(x) { abs(ranksTheor[, 1] - x) } )
    SRDscoresTheor <- apply(ranksTheor, 1, sum )
    res[, i] <- SRDscoresTheor
}

res <- c(res)
plot(density(res))
mean(res)
var(res)

plot(dnorm(0:100, mean = mean(res), sd = sd(res)), type = 'l')
for (i in 1:length(SRDscores)) {
    abline(v = SRDscores[i], col = i)
}


df <- data.frame(SRDdistr = 100*pnorm(0:100, mean = mean(res), sd = sd(res)))
firstQuantileSRD <- quantile(res, c(0.05, 0.975))[1]
#secondQuantileSRD <- quantile(res, c(0.025, 0.975))[2]
medianSRD <- median(res)

library(ggplot2)

plotSRD  <- ggplot(df, aes(x=0:100, y=SRDdistr)) +
    scale_x_continuous(breaks = seq(0, 100, by = 10)) +
    scale_y_continuous(breaks = seq(0, 100, by = 10)) +
    xlab("Observed SRD scores") +
    ylab("Observed SRD scores") +
    geom_vline(xintercept = firstQuantileSRD, col = "grey", linetype = "dashed") +
    #geom_vline(xintercept = secondQuantileSRD, col = "grey", linetype = "dashed") +
    geom_vline(xintercept = medianSRD, col = "grey", linetype = "dashed") +
    geom_line(col="grey") +
    geom_segment(aes(x=SRDscores[1], y = 0, xend = SRDscores[1], yend = SRDscores[1])) +
    geom_segment(aes(x=SRDscores[2], y = 0, xend = SRDscores[2], yend = SRDscores[2])) +
    geom_segment(aes(x=SRDscores[3], y = 0, xend = SRDscores[3], yend = SRDscores[3])) +
    geom_segment(aes(x=SRDscores[4], y = 0, xend = SRDscores[4], yend = SRDscores[4])) +
    geom_segment(aes(x=SRDscores[5], y = 0, xend = SRDscores[5], yend = SRDscores[5])) +
    geom_segment(aes(x=SRDscores[6], y = 0, xend = SRDscores[6], yend = SRDscores[6])) +
    annotate(geom="text", x=SRDscores[1] - 0.5, y = SRDscores[1], hjust = 1, label = "10-fold CV") +
    annotate(geom="text", x=SRDscores[2] - 0.5, y = SRDscores[2], hjust = 1, label = "10-fold CV 1 SE rule") +
    annotate(geom="text", x=SRDscores[3] - 0.5, y = SRDscores[3], hjust = 1, label = "BIC") +
    annotate(geom="text", x=SRDscores[4] - 0.5, y = SRDscores[4], hjust = 1, label = "Chull MSE") +
    annotate(geom="text", x=SRDscores[5] - 0.5, y = SRDscores[5], hjust = 1, label = "Chull VAF") +
    annotate(geom="text", x=SRDscores[6] - 0.5, y = SRDscores[6], hjust = 1, label = "IS") +
    theme_bw()

plotSRD

#ggsave("/home/turbo/work/mmsca/main/main_v1/plots/SRD.pdf", plot = plotSRD, height = 5, width = 10, unit='in', dpi = 300)

########################################################################################################################

#make plot
plot <- ggplot(results_to_analyse1, aes(y = tucker, x = method))
plot + geom_boxplot(aes(fill = n), outlier.size=0.3, lwd=0.2) +
    facet_grid(rows = vars(error), cols = vars(sparsity)) +
    geom_hline(yintercept = 0.85, col = "black", linetype = 2) +
    guides(fill=guide_legend(title="I")) +
    xlab("Method") + 
    scale_fill_manual(values=c("#484848", "#808080", "#D3D3D3")) +
    ylab("Tucker congruence") +
    coord_fixed(ratio = 26/5) +
    theme_bw() %+replace% theme(axis.text.x = element_text(angle = 60, hjust = 0.5, size=10),
                axis.text.y = element_text(size=10),
                strip.text = element_text(size=10),
                legend.text = element_text(size=10),
                legend.title = element_text(size=10),
                legend.key.size = unit(0.5, "cm"))


#ggsave("~/work/mmsca/mmscaDescAlgo/plots/multi_set_tucker.pdf", height = 10, width = 10, unit='in', dpi = 300)

#crop whitespace with pdfcrop 
#the margin option adds extra margin to the cropping, this is needed because else the legend is cropped off
#system("pdfcrop --margin 5 ~/work/mmsca/mmscaDescAlgo/plots/multi_set_tucker.pdf ~/work/mmsca/mmscaDescAlgo/plots/multi_set_tucker.pdf")

## weights found plot
#results_to_analyse2 <- results_to_analyse 
#rownames(results_to_analyse2) <- NULL
#colnames(results_to_analyse2)[10] <- "percentage_zeroCoefFound" 
#colnames(results_to_analyse2)[11] <- "percentage_nonZeroCoefFound"
#results_to_analyse2$id <- 1:3600
#levels(results_to_analyse2$n) <- c("I = 25", "I = 50", "I = 100")
#
#results_to_analyse2 <- reshape(results_to_analyse2, sep ="_", direction = "long", varying = 10:11, timevar="found") 
#results_to_analyse2
#results_to_analyse2$found
#results_to_analyse2$found <- as.factor(results_to_analyse2$found)
#levels(results_to_analyse2$found)
#levels(results_to_analyse2$found) <- c("Non-zero weights", "Zero weights")
#
#plot <- ggplot(results_to_analyse2, aes(y = percentage, x = method ))
#plot + geom_boxplot(aes(fill = found), outlier.size =0.3, lwd=0.2) +
#    facet_grid(error + sparsity ~ n) +
#    guides(fill=guide_legend(title="Weights")) +
#    xlab("Method") + 
#    ylab("Proportion of correctly identified weights") +
#    coord_fixed(ratio = 26/5) +
#    scale_fill_manual(values=c("#808080", "#D3D3D3")) +
#    theme_bw() %+replace% theme(axis.text.x = element_text(angle = 60, hjust = 0.5, size=10),
#                axis.text.y = element_text(size=10),
#                strip.text = element_text(size=10),
#                legend.text = element_text(size=10),
#                legend.title = element_text(size=10),
#                legend.key.size = unit(0.5, "cm"))


#crop whitespace with pdfcrop


#################################################################################


reps <- 5 
x <- 10 
#in total reps*x per conditions
conditionCombsReplicated <- conditionCombs
for (i in 1:(x-1)) {
    conditionCombsReplicated <- rbind(conditionCombsReplicated, conditionCombs)
}

conditionCombsReplicated
grouping <- (0:(nrow(conditionCombsReplicated)-1) %% nrow(conditionCombs)) + 1 #these groups belong together
test <- readRDS("./simulation_study_modelselection_results/model_selection_result_1")
Tmat <- test$datObjectRes[[1]]$X %*% test$datObjectRes[[1]]$P[, 1:3]
Tmat_hat <- test$datObjectRes[[1]]$X %*% test$estimationISRes[[1]]$estimationRes$W
tuckerCongruence(Tmat, Tmat_hat)

results_to_analyse <- vector("list", length = nrow(conditionCombs))

for (a in 1:nrow(conditionCombs)) {
    res <- matrix(NA, reps, 12)
    files <- 1:nrow(conditionCombsReplicated)  
    files <- files[grouping == a]
    resList <- vector("list", length = x)
    for (i in 1:x) {
        #fileName <- paste("./simulation_study_modelselection_results/model_selection_result_", files[i], sep = "")
        fileName <- paste("~/simulation_study_modelselection_multiset_results/model_selection_result_", files[i], sep = "")
        loadedFile <- readRDS(fileName)

        for (j in 1:reps) {
            X <-    loadedFile$datObjectRes[[j]]$X
            Tmat <- X %*% loadedFile$datObjectRes[[j]]$P[, 1:3]

            print(i, j)
            res[j, 1] <- tuckerCongruence(Tmat, X %*% loadedFile$estimationISRes[[j]]$estimationRes$W)
            res[j, 2] <- tuckerCongruence(Tmat, X %*% loadedFile$estimationBICRes[[j]]$estimationRes$W)
            res[j, 3] <- tuckerCongruence(Tmat, X %*% loadedFile$estimationCVRes[[j]]$estimationRes$W)
            res[j, 4] <- tuckerCongruence(Tmat, X %*% loadedFile$estimationCV1stErrorRuleRes[[j]]$estimationRes$W)
            res[j, 5] <- tuckerCongruence(Tmat, X %*% loadedFile$estimationCHullCVRes[[j]]$estimationRes$W)
            res[j, 6] <- tuckerCongruence(Tmat, X %*% loadedFile$estimationCHullVAFRes[[j]]$estimationRes$W)

            res[j, 7]  <- conditionCombs[a, 1]
            res[j, 8]  <- conditionCombs[a, 2]
            res[j, 9]  <- conditionCombs[a, 3]
            res[j, 10] <- conditionCombs[a, 4]
            res[j, 11] <- conditionCombs[a, 5]
            res[j, 12] <- conditionCombs[a, 6]

        }
        resList[[i]] <- res
    }
    resultsForCondition <-  do.call(rbind, resList)
    results_to_analyse[[a]] <- cbind(resultsForCondition, a)
}

results_to_analyse
results_to_analyse <- do.call(rbind, results_to_analyse)

apply(results_to_analyse, 2, median)

results_to_analyse < 100


