\name{BSphe}
\alias{BSphe}
\docType{data}
\title{Simulated phenotypic data}
\description{
Phenotypic data simulated with BSdat which has multiple observations per genotype
}
\usage{data(BSphe)}
\format{
A data frame with three columns:
\describe{
\item{ID:}{the ID for each individual. These are the same 250 
individuals appearing in BSdat, each with 4 replicates}
\item{Block:}{a factor with 4 levels indicating the replicate}
\item{phenotype:}{the phenotype value. It is generated from the same
QTL effects as in BSdat, with the additional block effect and random error}
}
}
\references{
Broman, KW and Speed TP. 2002. A model selection approach for the identification of quantitative trait loci in experimental crosses. JRSS-B 64:641-656.
}
\examples{
data(BSphe)

boxplot(BSphe$phenotype~BSphe$Block)
}
\keyword{datasets}
