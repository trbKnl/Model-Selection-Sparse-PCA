## Data repository of paper titled: Model selection techniques for sparse weight based PCA

### Abstract of the paper

Many studies make use of multiple types of data that are collected for the same set of samples, resulting in so-called multi-block data (for example multi-omics studies). A popular analysis framework is sparse principal component analysis (PCA) of the concatenated data. The sparseness in the component weights of these models is usually induced by penalties. A crucial factor in the use of such penalized methods, is a proper tuning of the regularization parameters used to give more or less weight to the penalties. In this paper we examine several model selection procedures to tune these regularization parameters for sparse PCA. The model selection procedures include cross-validation, BIC, Index of sparseness, and the convex Hull procedure. Furthermore, to account for the multi-block structure, we present a sparse PCA algorithm with a group LASSO penalty added to it, to either select or cancel out blocks of data in an automated way. Also the tuning of the group LASSO parameter is studied for the proposed model selection procedures. We conclude that when the component weights are to be interpreted, cross-validation with the one standard error rule is preferred; alternatively, if the interest lies in obtaining component scores using a very limited set of variables, the convex Hull, BIC, and index of sparseness are all suitable.

Authors: Niek C. de Schipper, Katrijn Van Deun

#### Description of files:

- analyse_results_modelselection_multiset.R: Creates Plots and Figures for the multi-set part of the paper
- analyse_results_modelselection.R: Creates Figures for the single set part of the paper
- data_example: Folder containing the data analysis of the included example 
- makeData.R: Function to generate a dataset given some parameters
- mmsca.cpp: Contains the c++ version of the algorithm described in the paper
- mmscaCppPackage: Contains the files and scripts to wrap mmsca.cpp in a package, needed if you want to run the simulation over multiple cores
- modelSelection.R: Contains function to perform model selection
- simulation_study_modelselection.R: Script to perform the simulation study
- tucker.R: Function to calculate the tucker congruences between columns of two matrices 

#### Not included (send an email if you want them)
- simulation_study_modelselection_multiset_results: Folder including the results from the simulation study
- simulation_study_modelselection_results: Folder including the results from the simulation study

