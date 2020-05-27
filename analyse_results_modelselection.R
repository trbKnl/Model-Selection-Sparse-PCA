##################################
# script to make the plots for PCA model selection
# This script makes Figures 1, 2 and 3

rm(list=ls())
source('./tucker.R')

#some testing
test <- readRDS("~/simulation_study_modelselection_results/model_selection_result_1")
names(test)
test$tuningISandBIC
test$tuningCVRes
test$tuningCHullVAFRes
test$estimationISRes

##########
# These were the conditions used
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
reps <- 5 
x <- 10 
#in total reps*x per conditions
conditionCombsReplicated <- conditionCombs
for (i in 1:(x-1)) {
    conditionCombsReplicated <- rbind(conditionCombsReplicated, conditionCombs)
}

conditionCombsReplicated
grouping <- (0:(nrow(conditionCombsReplicated)-1) %% nrow(conditionCombs)) + 1 #these groups belong together


results_to_analyse <- vector("list", length = nrow(conditionCombs))

for (a in 1:nrow(conditionCombs)) {
    res <- matrix(NA, reps, 24)
    files <- 1:nrow(conditionCombsReplicated)  
    files <- files[grouping == a]
    resList <- vector("list", length = x)

    for (i in 1:x) {
        #fileName <- paste("./simulation_study_modelselection_results/model_selection_result_", files[i], sep = "")
        fileName <- paste("~/simulation_study_modelselection_results/model_selection_result_", files[i], sep = "")
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

            res[j, 19] <- conditionCombs[a, 1]
            res[j, 20] <- conditionCombs[a, 2]
            res[j, 21] <- conditionCombs[a, 3]
            res[j, 22] <- conditionCombs[a, 4]
            res[j, 23] <- conditionCombs[a, 5]
            res[j, 24] <- conditionCombs[a, 6]

        }
        resList[[i]] <- res
    }
    resultsForCondition <-  do.call(rbind, resList)
    results_to_analyse[[a]] <- cbind(resultsForCondition, a)
}

results_to_analyse

results_to_analyse <- do.call(rbind, results_to_analyse)
results_to_analyse <- as.data.frame(results_to_analyse)
results_to_analyse 

colnames(results_to_analyse)  <- c("tucker_IS", "tucker_BIC", "tucker_10foldCV", "tucker_10foldCV1stdError", "tucker_CHullMSE", "tucker_CHullVAF", "zeroCoefFound_IS", "zeroCoefFound_BIC", "zeroCoefFound_10foldCV",
"zeroCoefFound_10foldCV1stdError", "zeroCoefFound_CHullMSE", "zeroCoefFound_CHullVAF",
"nonZeroCoefFound_IS", "nonZeroCoefFound_BIC", "nonZeroCoefFound_10foldCV",
"nonZeroCoefFound_10foldCV1stdError", "nonZeroCoefFound_CHullMSE", "nonZeroCoefFound_CHullVAF", 
"n", "p", "ncomp", "sparsity", "variances", "error", "conditionNumber")

results_to_analyse <- reshape(results_to_analyse, varying = 1:18, sep = "_", timevar = "method", direction = "long")


unique(results_to_analyse$method)
results_to_analyse$method[which(results_to_analyse$method == "10foldCV")] <- "10-fold CV"
results_to_analyse$method[which(results_to_analyse$method == "10foldCV1stdError")] <- "10-fold CV 1 SE rule"
results_to_analyse$method[which(results_to_analyse$method == "CHullMSE")] <- "CHull MSE"
results_to_analyse$method[which(results_to_analyse$method == "CHullVAF")] <- "CHull VAF"

results_to_analyse$conditionNumber <- as.factor(results_to_analyse$conditionNumber)

results_to_analyse$error <- as.factor(results_to_analyse$error)
results_to_analyse$n <- as.factor(results_to_analyse$n)
results_to_analyse$sparsity <- as.factor(results_to_analyse$sparsity)

levels(results_to_analyse$sparsity) <- c("Sparsity 30%", "Sparsity 80%")
levels(results_to_analyse$error) <- c("Error 5%", "Error 20%")
levels(results_to_analyse$n) <- c("25", "50", "100")

library(ggplot2)

# scripts to make the boxplots

# tucker plot
plot <- ggplot(results_to_analyse, aes(y = tucker, x = method))
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

#ggsave("~/work/mmsca/mmscaDescAlgo/plots/single_set_tucker.pdf", height = 10, width = 10, unit='in', dpi = 300)

#crop whitespace with pdfcrop
#system("pdfcrop --margin 5 ~/work/mmsca/mmscaDescAlgo/plots/single_set_tucker.pdf ~/work/mmsca/mmscaDescAlgo/plots/single_set_tucker.pdf")

################################################## weights found plot
results_to_analyse2 <- results_to_analyse 
rownames(results_to_analyse2) <- NULL
colnames(results_to_analyse2)[10] <- "percentage_zeroCoefFound" 
colnames(results_to_analyse2)[11] <- "percentage_nonZeroCoefFound"
results_to_analyse2$id <- 1:3600
levels(results_to_analyse2$n) <- c("I = 25", "I = 50", "I = 100")

results_to_analyse2 <- reshape(results_to_analyse2, sep ="_", direction = "long", varying = 10:11, timevar="found") 
results_to_analyse2
results_to_analyse2$found
results_to_analyse2$found <- as.factor(results_to_analyse2$found)
levels(results_to_analyse2$found)
levels(results_to_analyse2$found) <- c("Non-zero weights", "Zero weights")

plot <- ggplot(results_to_analyse2, aes(y = percentage, x = method ))
plot + geom_boxplot(aes(fill = found), outlier.size =0.3, lwd=0.2) +
    facet_grid(error + sparsity ~ n) +
    guides(fill=guide_legend(title="Weights")) +
    xlab("Method") + 
    ylab("Proportion of correctly identified weights") +
    coord_fixed(ratio = 26/5) +
    scale_fill_manual(values=c("#808080", "#D3D3D3")) +
    theme_bw() %+replace% theme(axis.text.x = element_text(angle = 60, hjust = 0.5, size=10),
                axis.text.y = element_text(size=10),
                strip.text = element_text(size=10),
                legend.text = element_text(size=10),
                legend.title = element_text(size=10),
                legend.key.size = unit(0.5, "cm"))

ggsave("~/work/mmsca/mmscaDescAlgo/plots/single_set_weightsFound.pdf", height = 10, width = 10, unit='in', dpi = 300)

#crop whitespace with pdfcrop
system("pdfcrop --margin 5 ~/work/mmsca/mmscaDescAlgo/plots/single_set_weightsFound.pdf ~/work/mmsca/mmscaDescAlgo/plots/single_set_weightsFound.pdf")


#################################################################################
#################################################################################
#################################################################################
#################################################################################
#################################################################################
#################################################################################
#################################################################################
#################################################################################
#################################################################################
#################################################################################
#################################################################################
#################################################################################
#################################################################################
#################################################################################
#################################################################################
#################################################################################
# The tucker congruence coefficients between T and That

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
        fileName <- paste("~/simulation_study_modelselection_results/model_selection_result_", files[i], sep = "")
        loadedFile <- readRDS(fileName)

        for (j in 1:reps) {
            X <-    loadedFile$datObjectRes[[j]]$X
            Tmat <- X %*% loadedFile$datObjectRes[[j]]$P[, 1:3]

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
results_to_analyse <- as.data.frame(results_to_analyse)
results_to_analyse

colnames(results_to_analyse)  <- c("tucker_IS", "tucker_BIC", "tucker_10foldCV", "tucker_10foldCV1stdError", "tucker_CHullMSE", "tucker_CHullVAF", "n", "p", "ncomp", "sparsity", "variances", "error", "conditionNumber")

results_to_analyse <- reshape(results_to_analyse, varying = 1:6, sep = "_", timevar = "method", direction = "long")


unique(results_to_analyse$method)
results_to_analyse$method[which(results_to_analyse$method == "10foldCV")] <- "10-fold CV"
results_to_analyse$method[which(results_to_analyse$method == "10foldCV1stdError")] <- "10-fold CV 1 SE rule"
results_to_analyse$method[which(results_to_analyse$method == "CHullMSE")] <- "CHull MSE"
results_to_analyse$method[which(results_to_analyse$method == "CHullVAF")] <- "CHull VAF"

results_to_analyse$conditionNumber <- as.factor(results_to_analyse$conditionNumber)

results_to_analyse$error <- as.factor(results_to_analyse$error)
results_to_analyse$n <- as.factor(results_to_analyse$n)
results_to_analyse$sparsity <- as.factor(results_to_analyse$sparsity)

levels(results_to_analyse$sparsity) <- c("Sparsity 30%", "Sparsity 80%")
levels(results_to_analyse$error) <- c("Error 5%", "Error 20%")
levels(results_to_analyse$n) <- c("25", "50", "100")


plot <- ggplot(results_to_analyse, aes(y = tucker, x = method))
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


#ggsave("/home/turbo/work/mmsca/main/main_v1/plots/single_set_T_That_tucker.pdf", height = 10, width = 10, unit='in', dpi = 300)

###crop whitespace with pdfcrop
#system("pdfcrop --margin 5 /home/turbo/work/mmsca/main/main_v1/plots/single_set_T_That_tucker.pdf /home/turbo/work/mmsca/main/main_v1/plots/single_set_T_That_tucker.pdf")




