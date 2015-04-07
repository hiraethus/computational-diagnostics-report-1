source("http://bioconductor.org/biocLite.R")
biocLite() # install all the main biocLite packages
biocLite("limma")
biocLite("GEOquery")
biocLite("affy")
biocLite("affyPLM")

library("limma")
library("GEOquery")
library("affy")
library("affyPLM")

# retrieve geneset
G.ExpressionSet <- GEOquery::getGEO("GSE24249", GSEMatrix=TRUE)[[1]] # this is 51.8MB - will take some time to download

# already normalised otherwise would use
# affyPLM::normalize.ExpressionSet.quantiles(G.ExpressionSet)

# design matrix
design <- cbind(case = c(1,1,1,0,0,0),
                control = c(0,0,0,1,1,1))

fit <- limma::lmFit(G.ExpressionSet, design)
fit <- limma::eBayes(fit)

result <- topTable(fit, coef=2, n=length(featureNames(G.ExpressionSet)), adjust="fdr")
