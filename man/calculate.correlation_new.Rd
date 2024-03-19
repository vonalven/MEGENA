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
calculate.correlation_new(datExpr,
                          doPerm            = 100,
                          num.cores         = 8,
                          method            = "pearson",
                          use.obs           = "pairwise.complete.obs",
                          FDR.cutoff        = 0.05,
                          direction         = "absolute",
                          rho.thresh.step   = 0.0025,
                          output.permFDR    = TRUE,
                          output.corTable   = TRUE,
                          saveto            = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{datExpr}{Gene expression data matrix (rows = probes/genes, columns = samples).}
  \item{doPerm}{Number of permutations to perform. If \code{doPerm = FALSE} (or \code{doPerm = 0}), calculates BH FDR p-values instead of permutation based FDR.}
  \item{num.cores}{Number of cores to use in parallelization.}
  \item{method}{Correlation method to be passed to the \code{method} argument of the \code{cor} function.}
  \item{use.obs}{Method to handle missing values when computing correlation, to be passed to \code{use} argument of the \code{cor} function.}
  \item{FDR.cutoff}{FDR threshold to subset only significant correlations returned in the output list object.}
  \item{direction}{If direction = "absolute", absolute correlation values are considered. If direction = "positive" or "negative", signed correlation is considered, and significance is calculated for positive/negative values.}
  \item{rho.thresh.step}{Used in the permutation-approach to compute correlation. This is the step used to define a vector or correlation thresholds at which FDR empirical p-values are estimated.}
  \item{output.permFDR,output.corTable}{TRUE/FALSE to choose to output permutation FDR and correlation tables.}
  \item{saveto}{folder to output results.}
}
\details{
  Using a permutation approach to compute global empirical FDR pValue for a set of correlation thresholds has the advantage that pValues are computed without prior assumptions on the data distribution.
  This has however the disadvantage that a global threshold is determined. A permutation approach per-edge (or gene-gene pair) would be ideal but given the huge amount of edges to test this is very computational expensive.
  Some additional notes on the permutation approach:

  The permutation approach operates under the assumptions that:
  \itemize{
    \item The observed "real" correlation are a combination of False Positive (FP) and True Positive (FP), i.e. FP + TP to optimize computations we split this computations in N permutations.
    \item The permuted correlations are FP.
    \item The computed FDR range must be [0, 1]. For low correlation thresholds, often an empirical FDR > 1 is observed (this is because at low correlation there is no linear trend in the data that are noisy). In this cases, FRSs are set to 1. The FDR also relies on the assumption that higher correlation thresholds should be associated with lower FDR pValues. To ensure monotonicity, a correction of the empirical FDR is set in place, where if the FDR of the previous sorted correlation threshold step is lower than the next one (for a higher correlation), it is assigned to the previous one.
    \item Combining this, a global False Discovery Rate (FDR) pValue is computed for each orrelation threhsold os step rho.thresh.step as FDR = FP / (FP + TP). Ideally, a similar permutation approach could be used to compute an empirical pValue for each edge (gene-gene correlation pair). This was the original adapted implementation here, available in a previous commit of this forked MEGENA branch, but it is too computationally expensive.
  }

  When does it make sense to use absolute, positive or negative directions?
  \itemize{
    \item absolute: to include all the genes in the correlation. In the PFN network, edges are added iteratively to an initially empty graph, starting with the most important edges (highest absolute correlation value), to impose planarity. This will define cluster with highly interconnected genes independently on the correlation sign.
    \item positive: to build networks with positively-correlated gene pairs (edges) only. This can be useful if we want to identify modules of strongly co-expressed (positively correlated) genes only.
    \item negative: to build networks with negatively-correlated gene pairs (edges) only. Biologically, this doesn't make a lot of sense as those modules would only include genes that are anti-correlated among themselves. For direction = "negative", we suggest to use \code{d.func = function(x) 1 - abs(x)} in function \code{do.MEGENA}.
  }
}
\value{
 output a list with significant gene-gene correlation pairs (edges) and FDR table.
}
\author{
Won-Min Song, adapted by Giacomo von Alvensleben
}
\examples{
# test simplest case of planar network (a 3-clique).
data(Sample_Expression)
calculate.correlation_new(datExpr[1:100, ], doPerm = 5)
}