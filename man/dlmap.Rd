\name{dlmap}
\alias{dlmap}
\alias{plot.dlmap}
\alias{summary.dlmap}
\alias{profileplot}
\alias{profileplot.dlmap}
\title{ Perform DLMapping}
\description{
Fits the iterative algorithm for DLMapping. Reads in data, performs detection and localization stages and outputs summary of selected QTL effects. 
}
\usage{
dlmap(object, phename, baseModel, algorithm=c("asreml", "lme"), fixed = NULL, 
random = NULL, rcov = NULL, sparse = NULL, pedigree, 
seed = 1, maxit=60, n.perm = 0, multtest=c("holm", "bon"), 
alpha = 0.05, filestem = "dl", ...)

\method{plot}{dlmap}(x, chr, max.dist, qcol="light blue", mcol="red", pcol="purple", marker.names=FALSE, \dots)

profileplot(object, \dots)
\method{profileplot}{dlmap}(object, chr, marker.names=TRUE, QTLpos=TRUE, pch=20, \dots)

\method{summary}{dlmap}(object, \dots)
}
\arguments{
  \item{object}{For input to dlmap, dlcross object created by create.dlcross; for plotting functions, object of class \code{dlmap}}
  \item{phename}{ Response variable name}
  \item{algorithm}{Indicates whether to fit mixed models using asreml function or lme function (using packages asreml or nlme)}
  \item{fixed}{ A formula object specifying the fixed effects part of the 
base model, with the terms, separated by + operators, on the right of 
a ~ operator. There is no left side to the ~ expression. If no fixed effect 
is specified, the model defaults to ~1, i.e. intercept only.} 
  \item{random}{ A formula object specifying the random effects part of the 
base model, with the terms, separated by + operators, on the right of 
a ~ operator. See asreml documentation for more detail. }
  \item{rcov}{ A formula object specifying the error structure of the model, with the terms, separated by + operators, on the right of a ~ operator. See asreml documentation for more detail. }
  \item{sparse}{ A formula object specifying the fixed effects to be 
absorbed, with the terms, separated by + operators, on the right of 
a ~ operator. See asreml documentation for more detail. }
  \item{baseModel}{ An alternative to specifying \code{fixed}, \code{random}, 
\code{sparse}, and \code{rcov} separately. If a base model has already been 
fit in asreml-R for the phenotypic variation, this can be input directly}
  \item{pedigree}{ Either a pedigree object consisting of three columns or a kinship matrix with number of rows and columns equivalent to the number of unique genotypes. For the pedigree, the first 
column is the individual ID, then the mother's ID and the father's ID. 
The name of the ID variable in the first column must match the \code{idname} variable. Rownames and colnames for the kinship matrix must match the \code{idname} variable.  }
  \item{seed}{ Random number seed. Default=1}
  \item{n.perm}{ Number of permutations used to get adjusted p-values at 
each iteration of detection. If \code{n.perm=0} (default) the Holm 
correction is used. }
  \item{multtest}{ Correction used for multiple testing. If n.perm>0 will use 
  permutation, but otherwise can choose between Holm and Bonferroni}
  \item{alpha}{ Significance level for testing }
  \item{filestem}{ Stem to add to names of any files generated in DL Mapping 
process. Default="dl"}
  \item{maxit}{ Maximum number of iterations to attempt for convergence of \code{lme}}
  \item{\dots}{additional arguments}
  \item{x}{object of class \code{dlmap}}
 \item{chr}{character string naming the subset of chromosomes to plot; or, if numeric, a vector of chromosome indices}
 \item{max.dist}{a numerical value in cM determining the distance the genetic map should be subsetted by}
 \item{qcol}{colour of intervals surrounding QTL (see \code{par} for colour options)}
 \item{mcol}{colour of QTL flanking markers (see \code{par} for colour options)}
 \item{pcol}{colour of QTL positions (see \code{par} for colour options)}
  \item{marker.names}{logical value. For profileplot, if \code{TRUE} then marker names are plotted along the top of the profileplot. Defaults to \code{TRUE}. For plot function, if \code{TRUE} then flanking marker names are highlighted. Defaults to \code{FALSE}}
 \item{QTLpos}{logical value. if \code{TRUE} then QTL positions are indicated with vertical lines in profileplot. Defaults to \code{TRUE}}
 \item{pch}{Character to be used for points in plotting; default is solid circle}
}
\details{
There are two versions of the main function, which use different engines to 
fit the linear mixed models which form the framework of the algorithm.
Which is used depends on the value of the argument \code{algorithm}. 
\code{algorithm="asreml"} provides a much more general implementation of the
DLMapping algorithm and is the preferred method of analysis. 
\code{algorithm="lme"} is more restricted in its capabilities, in 
that it cannot model random effects or covariance structure, cannot handle 
more than 200 markers, and only allows for a single phenotypic observation
per genotype. Also, permutation has not been implemented for this function
because it is very slow. However, this version will fit the basic 
algorithm and is useful should a license for ASReml not be available. 

In studies where the number of genetic markers is much larger than the number 
of phenotyped individuals, we reduce the dimension of the analysis to the number
of genetic lines used in the analysis multiplied by the number of chromosomes 
in the genetic map. This is done in a similar
manner to \pkg{wgaim}, with thanks to Julian Taylor and Ari Verbyla for 
the suggestion. The transformation of the genetic data
reduces the time for computational analysis for high-dimensional data and is 
particularly useful in association analysis. 

This version of wgaim allows high dimensional marker information to be analysed. A simple transformation of the collated high dimensional marker set shows that it may be reduced to the number of genetic lines used in the analysis. This transformation is internal to the wgaim.asreml call and users can now expect a considerably large acceleration in the performance of wgaim.



In the asreml version, there are two options for specifying the model for 
phenotypic variation. 
The individual model components can either be input directly as they would be 
in an ASReml call, or a previous model (\code{baseModel}) output from ASReml 
can be input and the components will be retrieved from it. The latter 
formulation may be useful if prior phenotypic modelling has taken place. Note
that in either case, variables appearing in the rcov statement must be 
ordered appropriately in the dataset. For example, if 
\command{rcov=~ar1(Column):ar1(Row)} the data must be sorted as 
\emph{Row} within \emph{Column}. 

Missing values in \command{asreml} are replaced with zeros, so it is important
to centre the covariate in question. This is done for all genotypes when 
\code{algorithm="asreml"}. Thus individuals with phenotypic but not
genotypic data, which play important roles in field trials, may be included
safely. When \code{algorithm="lme"} these individuals cannot be included, so the
default behavior is to omit observations with missing values. 

It is recommended that \code{n.perm} be set to 0 for initial exploratory 
analysis, as the permutation analysis may be lengthy. The Holm 
correction is used to adjust for the number of chromosomes under 
consideration at each detection stage. While this is a conservative 
measure it seems to perform well in practice.

Two files are output with names set by the argument \code{filestem}, which
has a default value of "dl". The 
file "filestem.trace" contains ASReml licensing and likelihood convergence
output which otherwise would be dumped to the screen and possibly obscure
other messages. Errors, warnings and other messages will still appear on the 
screen. Some warnings which appear may be passed through from an ASReml call
and output on exit. These may generally be ignored. This file is not created
if \code{algorithm="lme"} is used.  

The file "filestem.det.log" is a record of iterations in the detection stage. 
For each iteration the REMLRT testing for genetic variation on each chromosome
is output, along with adjusted p-values, genomewide threshold and markers 
selected as fixed effects. The p-values are corrected for the
number of chromosomes tested either by the Holm correction or by 
permutation. If the number of permutations (\code{n.perm}) is greater than 
0, then for the Xth iteration an additional file "filestem.permX" will be 
created which contains the test statistics for the permuted datasets. 
See the accompanying vignette for an example of how to interpret the ".det.log" file. 

If the type of cross is not "other", the plotting function plots the genetic linkage map for a selection of chromosomes. Indicates marker
locations, marker names, and detected QTL positions and associated flanking 
markers obtained from a \code{dlmap} fit. This function relies upon \code{link.map.cross}, which was written by Julian
Taylor for the \code{wgaim} package. It is built upon here by adding QTL 
regions and estimated positions to the map.

The function \code{plot.dlmap} provides a neat visual display of 
chromosomes. If no QTL are detected, only the linkage map will be plotted; 
otherwise detected QTL will be placed at their estimated positions 
and the intervals around them (and flanking markers) will be highlighted. 
If a subset of chromosomes are plotted and detected QTLs exist outside that 
subset a warning will be given that QTLs have been omitted from the display. 

The arguments \code{mcol}, \code{qcol} and \code{pcol} have been added for
personal colour highlighting the flanking markers, QTL regions and QTL 
positions respectively. The procedure may also be given the usual \code{col}
argument which will be passed on to the other markers. 

In order to ensure that all marker names are displayed 
without vertical overlap, the default value of the "cex" 
parameter passed to "text" should be manipulated. 
For large maps with many chromosomes, marker names and adjacent chromosomes
will overlap horizontally. In this case it is suggested that the user 
horizontally maximize the plotting window to remove overlap, or subset the
chromosomes displayed. 

The profileplot function plots the Wald statistic profile for each chromosome with detected QTL on first
interval mapping scan. Indicates marker locations, marker names, and detected 
QTL positions obtained from a \code{dlmap} fit. It provides a neat visual display of 
the Wald profile for chromosomes with detected QTL. If no QTL are detected, 
nothing will be plotted. Otherwise, the Wald profile will be plotted by cM
position of points on the interval mapping grid. Marker names will be 
displayed at the appropriate positions along the top of the plot. Vertical
lines will mark the position of detected QTL.

The summary function outputs a summary of a dlmap object and detected QTL.  It primarily prints the summary table
computed from \code{dlmap}. This includes the 
chromosome QTL are detected on, estimated positions, flanking markers, QTL
effects and standard deviations, Z-ratio and p-value. 
}
\value{
  \item{Summary}{Table with one row per QTL detected, columns for which 
chromosome the QTL is on, its position (cM), flanking markers, additive 
(dominant) effect and standard deviation, Z-ratio and p-value.}
  \item{no.qtl }{ Total number of QTL detected on all chromosomes}
  \item{final.model }{Object of class \code{asreml} for final model 
containing all terms in the base model, as well as effects for every QTL 
detected at the appropriate locations. No random effects for markers are fit}
  \item{profile}{If QTL are detected on \emph{C} chromosomes, this is a list 
with \emph{C} elements, each a matrix with 2 rows and a column for each 
position on the chromosome. The first row contains the cM position; the second 
row contains the Wald statistic for the model fit in the localization stage}
  \item{input}{Original dlcross input to the analysis}
}
\references{ Huang, BE and George, AW. 2009. Look before you leap: A new approach to QTL mapping. TAG 119:899-911

B. Emma Huang, Rohan Shah, Andrew W. George (2012). dlmap: An R Package for Mixed Model QTL and Association Analysis. Journal of Statistical Software 50(6): 1-22. URL http://www.jstatsoft.org/v50/i06/.
}
\author{ Emma Huang and Andrew George; Julian Taylor }
\seealso{ \code{\link[dlmap]{dlcross}} }
\examples{
\dontrun{
data(BSdat)
data(BSphe)

# Convert cross object to DL Mapping format
dl.in1 <- dlcross(format="rqtl", genobj=BSdat, idname="ID", fixpos=1)

# Analyze data
BSdl <- dlmap(object=dl.in1, algorithm="lme", phename="phenotype", filestem="BS")

plot(BSdl)

# With additional phenotypic data
dl.in2 <- dlcross(format="rqtl", genobj=BSdat, pheobj=BSphe, idname="ID", step=5)
BSph <- dlmap(object=dl.in2, algorithm="asreml", phename="phenotype", env=TRUE, random=~Block)

profileplot(BSph)
summary(BSph)
}
}
