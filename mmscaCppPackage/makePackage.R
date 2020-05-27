####################################################
# Create a package for mmsca cpp
#####################################################

library(testPackage)

require(Rcpp)
require(RcppArmadillo)
require(devtools)

# This works:
# This function create a skeleton package 
# it works with RcppArmadillo
# copy the .cpp file manually in the directory
RcppArmadillo.package.skeleton(name = "testPackage", list = character(), 
     environment = .GlobalEnv, path = ".", force = FALSE, 
     code_files = character(), example_code = TRUE)

#This is needed so the cpp functions will be exported
#when loading the package
Rcpp::compileAttributes(pkgdir="./testPackage")

#build the package with build from devtools
#install the package from the commandline with:
build(pkg="./testPackage")
#install.packages("./testPackage_1.0.tar.gz")
library(testPackage)



