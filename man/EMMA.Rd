\name{EMMA}
\alias{EMMA}
\title{Function for end-member modelling analysis.}
\description{A multivariate data set (m samples composed of n variables) is decomposed
by eigenspace analysis and modelled with a given number of end-members 
(q). Several steps of scaling, transformation, normalisation, eigen 
space decomposition, factor rotation, data modelling and evaluation are
performed.}
\usage{EMMA(X, q, lw, c, Vqn, EM.ID, classunits, ID, rotation = "Varimax", 
    plot = FALSE, legend, ..., pm = FALSE)}
\arguments{
  \item{X}{Numeric matrix with m samples (rows) and n variables (columns).}
  \item{q}{Numeric scalar with number of end-members to be modelled.}
  \item{lw}{Numeric scalar with the weight tranformation limit, i.e. quantiles, 
cf. Klovan & Imbrie (1971); default is 0.}
  \item{c}{Numeric scalar specifying the constant sum scaling parameter, e.g. 1, 
100, 1000; default is 100.}
  \item{Vqn}{Numeric matrix specifying optional unscaled user-defined end-member
loadings. If provided, these are used instead of model-derived ones.}
  \item{EM.ID}{Character vector with end-member names. If present, these will be set
as row-names of the output data set and used in the legend text if
a legend is enabled.}
  \item{classunits}{Numeric vector, optional class units (e.g. micrometers or phi-units) of 
the same length as columns of X.}
  \item{ID}{Numeric or character vector, optional sample IDs of the same
length as rows of X.}
  \item{rotation}{Character scalar, rotation type, default is "Varimax" (cf. Dietze et 
al., 2012). One out of the rotations provided in GPArotation is 
possible (cf. \code{\link{rotations}}).}
  \item{plot}{Logical scalar, optional graphical output of the results, default is 
FALSE. If set to TRUE, end-member loadings and end-member scores are 
plotted.}
  \item{legend}{Character scalar, legend position (cf. \code{\link{legend}}). If
omitted, no legend will be plotted, default is no legend.}
  \item{\dots}{Additional arguments passed to the plot function. Since the function 
returns two plots some additional graphical parameters must be
specified as vector with the first element for the first plot
and the second element for the second plot. Legend contents and 
colours apply to both plots. If colours are specified, \code{colour} 
should be used instead of \code{col}. See example section for 
further advice.}
  \item{pm}{Logical scalar to enable pm.}
}
\details{The function values \code{$loadings} and \code{$scores} are redundant.
They are essentially the same as \code{$Vqsn} and \code{$Mqs}. However,
they are included for user convenience.}
\value{A list with numeric matrix objects.
\item{loadings}{Normalised rescaled end-member loadings.}
\item{scores}{Rescaled end-member scores.}
\item{Vqn}{Normalised end-member loadings.}
\item{Vqsn}{Normalised rescaled end-member loadings.}
\item{Mqs}{Rescaled end-member scores.}
\item{Xm}{Modelled data.}
\item{modes}{Mode class of end-member loadings.}
\item{Mqs.var}{Explained variance of end-members}
\item{Em}{Absolute row-wise model error.}
\item{En}{Absolute column-wise model error.}
\item{Rm}{Row-wise (sample-wise) explained variance.}
\item{Rn}{Column-wise (variable-wise) explained variance.}
\item{ol}{Number of overlapping end-members.}}
\references{Dietze E, Hartmann K, Diekmann B, IJmker J, Lehmkuhl F, Opitz S, 
Stauch G, Wuennemann B, Borchers A. 2012. An end-member algorithm for 
deciphering modern detrital processes from lake sediments of Lake Donggi 
Cona, NE Tibetan Plateau, China. Sedimentary Geology 243-244: 169-180.\cr
Klovan JE, Imbrie J. 1971. An Algorithm and FORTRAN-IV Program for 
Large-Scale Q-Mode Factor Analysis and Calculation of Factor Scores. 
Mathematical Geology 3: 61-77.
Miesch AT. 1976. Q-Mode factor analysis of geochemical and petrologic  
data matrices with constant row sums. U.S. Geological Survey 
Professsional Papers 574.}
\author{Michael Dietze, Elisabeth Dietze}



\seealso{\code{\link{test.parameters}}, \code{\link{rotations}},
\code{\link{eigen}}, \code{\link{nnls}}}
\examples{
## load example data and set phi-vector
data(X.artificial, envir = environment())
phi <- seq(from = 1, to = 10, length.out = ncol(X.artificial))

## perform EMMA with 5 end-members
EM <- EMMA(X = X.artificial, q = 5, lw = 0.05, c = 100, plot = TRUE)

## perform EMMA with 4 end-members and more graphical settings
EM <- EMMA(X = X.artificial, q = 4, lw = 0.05, c = 100, 
           plot = TRUE,
           EM.ID = c("EM 1", "EM 2", "EM 3", "EM 4"),
           classunits = phi,
           xlab = c(expression(paste("Class [", phi, "]")), "Sample ID"),
           legend = "topleft", 
           cex = 0.7,
           colour = colors()[c(441, 496, 499, 506)])
}

\keyword{EMMA}
