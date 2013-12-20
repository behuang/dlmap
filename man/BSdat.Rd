\name{BSdat}
\alias{BSdat}
\docType{data}
\title{Simulated data for a backcross}
\description{
Dataset simulated according to Broman and Speed (2002). 
}
\usage{data(BSdat)}
\format{
An object of class \code{cross}. List with two components:
\itemize{
\item \code{geno}    is a list with elements corresponding to chromosomes. 
\code{names(geno)} contains the names of the chromosomes. There are two 
components for each chromosome: \code{data}, a matrix whose rows are 
individuals and whose columns are markers, and \code{map}, a vector of 
marker positions in cM. There is no missing data, and genotypes are coded as 1=AA, 2=AB
\item \code{pheno}   is a data frame of size (250 x 2) containing the trait data. The first trait is generated from a random normal distribution with variance 1.0 and mean determined by the QTL genotypes as described below. The second trait
is an ID for each individual
}
}
\details{
The data was generated for a sample size of 250 from a map with 9 chromosomes.
Each chromosome had length 100 cM and contained 11 equally spaced markers 
(spaced 10 cM apart). The background phenotypic variation was 1.0 and there 
was no missing data. The QTL were located as follows:
\itemize{
\item Chr 1: 2 QTL located at 30 and 70 cM, both with effect size of 0.76
\item Chr 2: 2 QTL located at 30 and 70 cM with effect size of 0.76 and -0.76
\item Chr 3: 1 QTL located at 50 cM with effect size of 0.76
\item Chr 4: 1 QTL located at 30 cM with effect size of 0.76
\item Chr 5: 1 QTL located at 0 cM with effect size of 0.76
}
}
\references{
Broman, KW and Speed TP. 2002. A model selection approach for the identification of quantitative trait loci in experimental crosses. JRSS-B 64:641-656.
}
\examples{
data(BSdat)
library(qtl) 

# Summary of chromosomes and markers
nchr(BSdat)
nmar(BSdat)

# linkage map of data
plot.map(BSdat)

# interval mapping
BSgp <- calc.genoprob(BSdat, step=2)
BSim <- scanone(BSgp)

# composite interval mapping
BScim <- cim(BSgp, n.marcov=5, method="hk")
 
# LOD profile from CIM
plot(BScim)

# LOD threshold for 5 cofactors from paper
abline(h=3.56)
}
\keyword{datasets}
