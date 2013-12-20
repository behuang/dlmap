\name{data access}
\alias{ngen}
\alias{nphen}
\alias{nmrk}
\alias{ngen.dlcross}
\alias{ngen.dlmap}
\alias{nphen.dlcross}
\alias{nphen.dlmap}
\alias{nmrk.dlcross}
\alias{nmrk.dlmap}
\title{ Data summaries of dlcross and dlmap objects}
\description{
Access the number of unique genotyped individuals; unique phenotyped individuals; and number of markers on each map chromosome}
\usage{
ngen(object, \dots)

\method{ngen}{dlmap}(object, \dots)

\method{ngen}{dlcross}(object, \dots)

nphen(object, \dots)

\method{nphen}{dlmap}(object, \dots)

\method{nphen}{dlcross}(object, \dots)

nmrk(object, \dots)

\method{nmrk}{dlmap}(object, \dots)

\method{nmrk}{dlcross}(object, \dots)
}
\arguments{
  \item{object}{ Object of type dlcross or dlmap}
  \item{\dots}{ Additional arguments}
}
\value{ 
ngen returns the number of unique genotyped individuals. 
nphen returns the number of unique phenotyped individuals - generally greater than or equal to ngen because of replicates. 
nmrk returns a vector indicating the number of markers on each chromosome. 
}
\references{ Huang, BE and George, AW. 2009. Look before you leap: A new approach to QTL mapping. TAG 119:899-911

B. Emma Huang, Rohan Shah, Andrew W. George (2012). dlmap: An R Package for Mixed Model QTL and Association Analysis. Journal of Statistical Software 50(6): 1-22. URL http://www.jstatsoft.org/v50/i06/.}
\author{ Emma Huang and Andrew George }
\examples{
# load dataset
data(BSdat)
data(BSphe2)

\dontrun{
# convert data to dlmap format
dl.in1 <- dlcross(format="rqtl", genobj=BSdat, idname="ID", fixpos=1)

ngen(dl.in1)
nphen(dl.in1)
nmrk(dl.in1)

# convert data with separate phenotypic trait file
dl.in2 <- dlcross(format="rqtl", genobj=BSdat, pheobj=BSphe2, idname="ID", step=5)

ngen(dl.in2)
nphen(dl.in2)
nmrk(dl.in2)
}
}
