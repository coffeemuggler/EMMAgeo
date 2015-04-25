\name{test.parameters}
\alias{test.parameters}
\title{Function to evaluate influence of model parameters.}
\description{All possible combinations of number of end-members and weight 
transformation limits are used to perform EMMA. The function returns 
matrices of absolute and relative measures of individual model performance.}
\usage{test.parameters(X, q, lw, c, rotation = "Varimax", plot, legend, 
    ..., pm = FALSE)}
\arguments{
  \item{X}{Numeric matrix with m samples (rows) and n variables (columns).}
  \item{q}{Numeric vector of length two, specifying the minimum and maximum number 
of end-members to be modelled.}
  \item{lw}{Numeric vector specifying the weight tranformation limit, i.e. 
quantile; default is 0.}
  \item{c}{Numeric scalar specifying the constant sum scaling parameter, e.g. 1, 
100, 1000; default is 0.}
  \item{rotation}{Character scalar, rotation type, default is "Varimax" (cf. Dietze et 
al., 2012). One out of the rotations provided in GPArotation is 
possible (cf. \code{\link{rotations}}).}
  \item{plot}{Character scalar, optional graphical output of the results. Specify
which tested parameter will be plotted: "mEm" (mean absolute row-wise 
error), "mEn" (mean absolute column-wise error), "mRm" (mean relative 
row-wise error), "mRn" (mean relative column-wise error), "mRt" (mean
relative total error), "ol" (number of overlapping end-members). All
plots except "ol" are colour-coded bitmaps of q, lw and the 
specified test parameter and line-plots the specified parameter vs. q.}
  \item{legend}{Character scalar, specifying legend position (cf. \code{\link{legend}}). 
If omitted, no legend will be plotted, default is no legend.}
  \item{\dots}{Additional arguments passed to the plot function. Since the function 
returns two plots (except for plot option "ol"), additional graphical
parameters must be specified as vector with the first element for the
first plot and the second element for the second plot. If graphical 
parameters are natively vectors (e.g. a sequence of colours), they must
be specified as matrices with each vector as a row. A legend can only
be added to the second plot. Colours only apply to the second plot as
well. If colours are specified, \code{colour} should be used instead of 
\code{col}. See example section for further advice.}
  \item{pm}{Logical scalar to enable pm.}
}
\details{The mean total explained variance mRt may be used to define a maximum 
number of meaningful end-members for subsequent modelling, e.g. as the 
number of end-members, which reaches the first local mRt maximum.\cr\cr
Overlapping is defined as one end-member having its mode within the
"area" of any other end-member, which is genetically not explainable.}
\value{A list with result objects
\item{mEm}{Absolute row-wise model error.}
\item{mEn}{Absolute column-wise model error.}
\item{mRm}{Mean row-wise explained variance.}
\item{mRn}{Mean column-wise explained variance.}
\item{mRt}{Mean total explained variance.}
\item{ol}{Number of overlapping end-member loadings.}
\item{q.max}{Maximum number of meaningful end-members.}}
\references{Dietze E, Hartmann K, Diekmann B, IJmker J, Lehmkuhl F, Opitz S, 
Stauch G, Wuennemann B, Borchers A. 2012. An end-member algorithm for 
deciphering modern detrital processes from lake sediments of Lake Donggi 
Cona, NE Tibetan Plateau, China. Sedimentary Geology 243-244: 169-180.}
\author{Michael Dietze, Elisabeth Dietze}



\seealso{\code{\link{EMMA}}}
\examples{
## load example data set
data(X.artificial, envir = environment())

## truncate the data set for faster computation
X.trunc <- X.artificial[1:20,]

## define test parameters
q <- 2:8 # number of end-members
lw <- seq(from = 0, to = 0.3, by = 0.1)

## test parameter influence and plot mean total explained variance
TP <- test.parameters(X = X.trunc, q = q, lw = lw, plot = "mRt",
                      legend = "bottomright", cex = 0.7,
                      colour = rgb((1:7) / 7, 0.9, 0.2, 1))

## show maximum number of end-members
TP$q.max
}

\keyword{EMMA}
