\name{do.MEGENA}
\alias{do.MEGENA}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
MEGENA clustering + MHA
}
\description{
multiscale clustering analysis (MCA) and multiscale hub analysis (MHA) pipeline
}
\usage{
do.MEGENA(g,
          do.hubAnalysis = TRUE,
          mod.pval       = 0.05,
          hub.pval       = 0.05,
          remove.unsig   = TRUE,
          min.size       = 10,
          max.size       = 2500,
          d.func         = function(x) {1-x},
          doPar          = FALSE,
          num.cores      = 4,
          n.perm         = 100,
          singleton.size = 3,
          save.output    = FALSE)
}

\arguments{
  \item{g}{igraph object of PFN.}
  \item{do.hubAnalysis}{TRUE/FALSE indicating to perform multiscale hub analysis (MHA) in downstream. Default is TRUE.}
  \item{mod.pval}{cluster significance p-value threshold w.r.t random planar networks}
  \item{hub.pval}{hub significance p-value threshold w.r.t random planar networks}
  \item{remove.unsig}{TRUE/FALSE indicating to remove insignificant clusters in MHA.}
  \item{min.size}{minimum cluster size}
  \item{max.size}{maximum cluster size}
  \item{d.func}{Distance function. See details below.}
  \item{doPar}{TRUE/FALSE indicating parallelization usage}
  \item{num.cores}{number of cores to use in parallelization.}
  \item{n.perm}{number of permutations to calculate hub significance p-values/cluster significance p-values.}
  \item{singleton.size}{Minimum module size to regard as non-singleton module. Default is 3.}
  \item{save.output}{TRUE/FALSE to save outputs from each step of analysis}
}
\details{
Performs MCA and MHA by taking PFN as input. Returns a list object containing clustering outputs, hub analysis outputs, and node summary table. 

Recommended d.func (see explanations in the next paragraph):
\itemize{
    \item for unsigned correlation: \code{d.func = function(x) 1-x}
    \item for signed correlation: \code{d.func = function(x) 1-abs(x)}
  }

As explained in the reference publication, the MCA algorithm incorporates three distinct criteria to identify locally coherent clusters while maintaining a globally optimal partition:
\enumerate{
    \item First, shortest path distances (SPD) are utilized to optimize within-cluster compactness.
    \item Second, local path index (LPI) is used to optimize local clustering structure. The LPI is a measure of the importance of nodes in a network, based on the number of walks of length 3 between pairs of nodes. These paths represent local connectivity patterns within the network.
    \item Third, overall modularity (Q) is employed to identify optimal partition.
  }

The edges weight (= correlation values) are used to compute SPDs for each pairwise combination of the nodes of a sub-network defined by connected vertexes (significant gene-gene correlation pairs), resulting in a symmetric matrix of SPD. The shortest path between two vertices in a graph is the path with the minimum sum of edge weights connecting them.
To compute meaningful SPDs, a measure of distance is used instead of simply the edges weight (= gene-gene correlation). The default d.func distance function is d = 1 - x, where x is an edge weight. By this way, if 2 genes are highly correlated (rho ~1), d = 1 - ~1 = ~0 and the distance is minimal. 
When using signed correlations, the default distance function is not appropriate anymore. in fact, if we have some edges with negative correlation, using d = 1 - x would produce distances larger than one, inflating the increase of SPD between 2 genes (computed as sum of edge distance metrics).
Let's imagine we have an edge with correlation -0.9. In this case d = 1 - (-0.9) = 1.9, indicating very far away genes. By using d = 1 - abs(x), the problem is solved: d = 1 - abs(-0.9) = 0.1, making those 2 genes very close and resulting in a correct SBD not deflating the compactness of sub-networks. 

Note: the edges weight are only used in this SPD-computation step, we should therefore not worry in other parts of the code about the sign of the correlation. 

Note: since the SPD is optimized and not for example the concordance of sub-network correlation sign, we can have modules with co-positively and negatively correlated genes. This simply means that those genes are very correlated between them, independently on the directionality of the correlation.
}
\value{
A series of output files are written in wkdir. Major outputs are,
\item{module.output}{outputs from MCA}
\item{hub.output}{outputs from MHA}
\item{node.summary}{node table summarizing clustering results.}
}
\author{
Won-Min Song, documentation adapted by Giacomo von Alvensleben
}
\examples{
\dontrun{
rm(list = ls())
data(Sample_Expression)
ijw <- calculate.correlation(datExpr[1:100,],doPerm = 2)
el <- calculate.PFN(ijw[,1:3])
g <- graph.data.frame(el,directed = FALSE)
MEGENA.output <- do.MEGENA(g = g,remove.unsig = FALSE,doPar = FALSE,n.perm = 10)
}
}