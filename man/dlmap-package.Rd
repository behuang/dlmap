\name{dlmap-package}
\alias{dlmap-package}
\alias{dlmap-pkg}
\docType{package}
\title{
DLMapping for QTL detection
}
\description{
QTL mapping in a mixed model framework with separate detection and 
localization stages. The former detects the number of QTL on each chromosome 
based on the genetic variation due to the grouped markers on the chromosome, 
while the latter stage uses this information to determine the most likely QTL 
positions. The mixed model can accommodate general fixed and random effects, 
including spatial effects in field trials and random pedigree effects. 
}
\details{
\tabular{ll}{
Package: \tab dlmap\cr
Type: \tab Package\cr
Version: \tab 1.13.1\cr
Date: \tab 2013-12-11\cr
License: \tab GPL 2\cr
}

The function \code{\link{dlcross}} is a constructor function for \code{dlcross}
objects to be input to \code{\link{dlmap}}. 
It will read in files in several different formats, including the \pkg{qtl}
\code{cross} format and two new formats to 
accommodate association mapping populations and designs with extensive 
phenotypic data. 

The primary function is \code{\link{dlmap}}, which performs the iterative 
algorithm to detect and position QTL on all chromosomes with significant 
genetic variation. This can accomodate sophisticated mixed models for phenotypic
variation in addition to the genetic modeling. 

Because ASReml-R is proprietary, we provide the option of using the \pkg{nlme} package to fit mixed models via the \code{algorithm} argument. However, there
are some features which are only accessible when the \code{algorithm} used is \code{asreml}.

The vignette included in this package gives more background on the methodology, input file structure, and examples of how to use each of the important functions
in the package. 

}
\author{
Emma Huang and Andrew George \cr 
Maintainer: Emma Huang <Emma.Huang@csiro.au>
}
\references{Huang, B.E. and George, A.W. 2009. Look before you leap: A new approach to QTL mapping. TAG 119:899-911

B. Emma Huang, Rohan Shah, Andrew W. George (2012). dlmap: An R Package for Mixed Model QTL and Association Analysis. Journal of Statistical Software 50(6): 1-22. URL http://www.jstatsoft.org/v50/i06/.
}
\seealso{
\code{\link[ibdreg]{ibdreg}}, \code{\link[nlme]{lme}}
}
