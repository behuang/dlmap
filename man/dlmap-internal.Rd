\name{dlmap-internal}
\alias{dldetect}
\alias{dllocalize}
\alias{dlmapdet}
\alias{dlmaploc}
\alias{dltest}
\alias{calcpos}
\alias{calc.genoprob2}
\alias{waldtest.asreml}
\alias{cintern}
\alias{dlcross.cross}
\alias{dlcross.dlmap}
\alias{dlcross.other}
\title{ Internal dlmap functions}
\description{
  Internal dlmap functions. These are generally not to be called by the user. 
}
\usage{

dlcross.cross(genobj, pheobj, idname, step, fixpos, estmap)
dlcross.dlmap(genobj, pheobj, mapobj, idname, genfile, mapfile, phefile, type, step, fixpos, estmap, ...)
dlcross.other(genobj, pheobj, mapobj, idname, genfile, mapfile, phefile)

dldetect(input, algorithm, filestem, ...)
dllocalize(input, algorithm, QTLperChr, ...)
dlmapdet(input, algorithm, s.chr, chrSet, prevLoc = NULL, ...)
dlmaploc(input, algorithm, s.chr, chrSet, prevLoc = NULL, ...)
dltest(input, algorithm, chrSet, prevLoc = NULL, ...)

calcpos(cross, step, fixpos)
calc.genoprob2(cross, pos, error.prob=1e-04, map.function=c("haldane", "kosambi", "c-f", "morgan"))

waldtest.asreml(object, cc)
cintern(cc, tau, vrb, sigma2)
}
\details{
The detect, localize, test and mapping functions form the basis of the 
DLMapping algorithm. The argument \code{algorithm} allows for versions using
either \command{lme} or \command{asreml} to fit mixed models.\cr 

The merge function is used to combine the genotype and environmental data,
and also imputes missing values in the genotype data according to the Viterbi
algorithm. See \code{\link[qtl]{fill.geno}} for more details. Note that 
individuals with no phenotypic response data are omitted in the merged 
dataset. Individuals with phenotypic but no genotypic data are retained, but
the genotypes are not imputed. 

calc.genoprob2 is a modification of a function from the R/qtl package to allow
for calculation of genotype probabilities at more general positions than a fixed
step size. 
}
\references{Huang, B.E. and George, A.W. 2009. Look before you leap: A New Approach to QTL Mapping. TAG 119:899-911

Broman et al. (2003) R/qtl: QTL mapping in experimental crosses. Bioinformatics 19:889-890

B. Emma Huang, Rohan Shah, Andrew W. George (2012). dlmap: An R Package for Mixed Model QTL and Association Analysis. Journal of Statistical Software 50(6): 1-22. URL http://www.jstatsoft.org/v50/i06/. 
}
\author{ Emma Huang and Andrew George}
\keyword{internal}
