\name{dlcross}
\alias{dlcross}
\alias{plot.dlcross}
\alias{summary.dlcross}
\title{ Constructor, summary and plotting functions for dlcross format}
\description{
Reads in objects in \code{cross} format and files or objects in \code{dlmap} format and converts them to dlmap input format. }
\usage{
dlcross(format = c("rqtl", "dlmap", "other"), genobj, pheobj, mapobj, idname="ID", genfile, phefile, mapfile, type, step=0, fixpos=0, estmap=TRUE, ...)
\method{plot}{dlcross}(x, chr, pheno.col, \dots)
\method{summary}{dlcross}(object, \dots)
}
\arguments{
 \item{format}{ See documentation for \code{\link[qtl]{read.cross}}. Also supports the 
input of an object of class \code{cross} ("rqtl" format), "dlmap" format 
(described below) and continuous data for association mapping populations 
("other" format)}
  \item{genobj}{ if format="rqtl", object of class \code{cross}; 
if format="dlmap" or "other" data frame containing genotypes}
  \item{pheobj}{ Data frame or matrix containing supplementary 
phenotypic or environmental data}
  \item{mapobj}{ if format="other" or "dlmap", 2/3-columned data-frame 
containing marker names, chromosomes and marker positions}
  \item{idname}{ Unique identifier variable name; will be used to match 
phenotypic and marker data }
  \item{genfile}{ if format="dlmap" see below}
  \item{mapfile}{ if format="dlmap" see below}
  \item{phefile}{ if format="dlmap" see below}
  \item{type}{Type of experimental cross; only necessary if format="dlmap" or 
"other"; if not input, will be assumed that data is from unrelated individuals}
  \item{step}{ Step size for localization stage, i.e. if \code{step=2}, grid 
of positions spaced 2 cM apart are considered for QTL locations. If 
\code{step=0} (default) positions are only located at markers.}
  \item{fixpos}{ Alternative to specifying a step size - if \code{fixpos=2}, 
2 evenly spaced positions between each marker are considered as QTL locations. 
If \code{fixpos=0} (default) positions are only located at markers. }
  \item{estmap}{Flag for whether to re-estimate the linkage map. Cannot be done 
for format="other"}
  \item{\dots}{ additional arguments}  
  \item{x}{input to plot function, object of class \code{dlcross}}
  \item{chr}{character string naming the subset of chromosomes to plot; if numeric, a vector of chromosome indices}
  \item{pheno.col}{vector of phenotypes for which to show barplots/histograms}
  \item{object}{input to summary function, object of class \code{dlcross}}
}
\details{
The main function constructs \code{dlcross} objects for input into 
\code{\link{dlmap}}. The \code{format} 
argument allows for compatibility with the \code{cross} format as supported
by \code{\link[qtl]{read.cross}} in R/qtl. 
In addition, format="dlmap" or "other" takes three files or 
objects in the following form. 

\itemize{
\item \code{genfile}: First row contains an identifier variable and marker names. 
Following rows contain genotype values for each individual
\item \code{phefile}: First row contains names of phenotypic and environmental traits. One of these traits must be an identifier for each individual. 
Following rows contain values of the traits for each individual. Note: This file may contain more observations than the file containing the marker data
\item \code{mapfile}: Contains marker names, chromosome groupings and positions
}

If a single set of trait values is available for each genotype, then 
phenotypic data will be input through the arguments \code{genobj} or \code{phefile} (depending on the file format). 
The argument \code{pheobj} allows for input of phenotypic data on replicates
or additional individuals which are not necessarily genotyped. 

Choosing format="other" allows for association mapping populations to be 
analyzed, and in this case the data can be input as with format="dlmap", 
but a genetic linkage map is not required. Hence the mapfile only needs the
first two columns of marker names and chromosome groupings. 

The plot function plots diagnostics summarizing the data. If the type of cross is not "other", will plot the genetic map for the cross. Will also plot a histogram or barplot of the first few phenotypic variables. Will only plot a maximum of three phenotypic variables; pheno.col can be 
used to select which are plotted. 

The summary function outputs a summary of data stored in dlcross object, including number of genotypes, number of phenotypes, number of phenotypic variables, number of chromosomes and markers per chromosome. For experimental crosses, based off of \code{cross} 
object summary; similar output for association mapping populations. 

} 
\value{ Object with class "dlcross" which can be input to dlmap. 
Contains the following elements:
\item{dfMerged}{Data frame to be used in dlMapping analysis}
\item{map}{Original genetic map}
\item{nphe}{Number of phenotypic traits}
\item{loc}{A flag for whether to run the localization stage}
\item{idname}{the ID name input to the function}
\item{mapp}{If format is not "other", genetic map augmented by imputed genotypes at grid of positions}
\item{genCross}{if format is not "other", rqtl cross object. CAUTION: if 
there are replicates, i.e., more phenotypic data than genotypic, this object
will not contain all of the phenotypic data for the sample}
}
\seealso{\code{\link[qtl]{plot.cross}}}
\references{ Huang, BE and George, AW. 2009. Look before you leap: A new approach to QTL mapping. TAG 119:899-911

B. Emma Huang, Rohan Shah, Andrew W. George (2012). dlmap: An R Package for Mixed Model QTL and Association Analysis. Journal of Statistical Software 50(6): 1-22. URL http://www.jstatsoft.org/v50/i06/.
}
\author{ Emma Huang and Andrew George }
\examples{
\dontrun{
# load dataset
data(BSdat)
data(BSphe)

# convert data to dlmap format
dl.in1 <- dlcross(format="rqtl", genobj=BSdat, idname="ID", fixpos=1)

# convert data with separate phenotypic trait file
dl.in2 <- dlcross(format="rqtl", genobj=BSdat, pheobj=BSphe, idname="ID", step=5)

plot(dl.in2)
summary(dl.in2)
}
}
