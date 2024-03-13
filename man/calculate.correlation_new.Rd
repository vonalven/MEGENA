\name{calculate.correlation_new}
\alias{calculate.correlation_new}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
correlation calculation
}
\description{
correlation analysis with FDR calculation. This is an adapted version of the original calculate.correlation function, that is orders of magnitude faster for the correlation calculation without permutations. The permutation-based approach to compute correlation pValues was redesigned to determine significance by genes instead of determining a global FDR threshold. Please notice that for the permutation approach, the more permutations, the better the estimation of the empirical pValue.
}
\usage{
calculate.correlation(datExpr,doPerm = 100,doPar = FALSE,num.cores = 8,method = "pearson",
FDR.cutoff = 0.05,n.increment = 100,is.signed = FALSE,
output.permFDR = TRUE,output.corTable = TRUE,saveto = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{datExpr}{gene expression data matrix}
  \item{doPerm}{Number of permutations to perform. If \code{doPerm = FALSE} (or \code{doPerm = 0}), calculates BH FDR p-values instead of permutation based FDR.}
  \item{doPar}{TRUE/FALSE logical variable to choose parallelization. Parallelization is utilized when BH FDR p-values are calculated for all pairs.}
  \item{num.cores}{number of cores to use in parallelization.}
  \item{method}{correlation method to be passed to \code{cor} for \code{method} argument.}
  \item{FDR.cutoff}{FDR threshold to output final results of significant correlations.}
  \item{n.increment}{When permutation is utilized, 0 <= |rho| <= 1 is broken down into n.increment to map each |rho| cutoff to respective FDR.}
  \item{is.signed}{TRUE/FALSE to indicate using signed/unsigned correlation.}
  \item{output.permFDR,output.corTable}{TRUE/FALSE to choose to output permutation indices and FDR table.}
  \item{saveto}{folder to output results.}
}
\details{
If \code{doPar = TRUE}, then num.cores are registered for PCP. 
}
\value{
 output is three column edgelist data.frame, third column being the weight.
}
\author{
Won-Min Song, adapted by Giacomo von Alvensleben
}
\examples{
# test simplest case of planar network (a 3-clique).
data(Sample_Expression)
calculate.correlation(datExpr[1:100,],doPerm = 5)
}