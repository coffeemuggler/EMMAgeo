\name{test.L}
\alias{test.L}
\title{Function to evaluate the initial cumulative explained variance.}
\description{This function performs eigenspace decomposition using the 
weight-transformed matrix W to determine the minimum number of end-members.
Depending on the number of provided weight transformation limits (lw) a 
single vector or a matrix is returned.}
\usage{test.L(X, lw, c, r.min = 0.95, plot = FALSE, legend, ..., pm = FALSE)}
\arguments{
  \item{X}{Numeric matrix with m samples (rows) and n variables (columns).}
  \item{lw}{Numeric vector specifying the weight tranformation limits, i.e. 
quantiles; default is 0.}
  \item{c}{Numeric scalar specifying the constant sum scaling parameter, e.g. 1, 
100, 1000; default is 100.}
  \item{r.min}{Numeric scalar, minimum value of explained variance to be reached by 
the end-members included, default is 0.95.}
  \item{plot}{Logical scalar, optional graphical output of the results, default is 
FALSE.,}
  \item{legend}{Character scalar, specify legend position (cf. \code{\link{legend}}). If
omitted, no legend will be plotted, default is no legend.}
  \item{\dots}{Additional arguments passed to the plot function. Use \code{colour} instead
of \code{col} to create different colours.}
  \item{pm}{Logical scalar to enable pm.}
}
\details{The results may be used to define a minimum number of end-members for
subsequent modelling steps, e.g. by using the Kaiser criterion, which 
demands a minimum number of eigenvalues to reach a squared R of 0.95.}
\value{A list with objects
\item{L}{Vector or matrix of cumulative explained variance.}
\item{q.min}{Vector with number of factors that passed r.min.
}}
\references{Dietze E, Hartmann K, Diekmann B, IJmker J, Lehmkuhl F, Opitz S, 
Stauch G, Wuennemann B, Borchers A. 2012. An end-member algorithm for 
deciphering modern detrital processes from lake sediments of Lake Donggi 
Cona, NE Tibetan Plateau, China. Sedimentary Geology 243-244: 169-180.}
\author{Michael Dietze, Elisabeth Dietze}




\examples{
## load example data set
data(X.artificial, envir = environment())

## create sequence of weight transformation limits
lw <- seq(from = 0, to = 0.2, 0.02)

## perform the test and show q.min
L <- test.L(X = X.artificial, lw = lw, c = 100, plot = TRUE)
L$q.min

## a visualisation with more plot parameters
L <- test.L(X = X.artificial, lw = lw, c = 100, plot = TRUE, 
            ylim = c(0.5, 1), xlim = c(1, 7), 
            legend = "bottomright", cex = 0.7)

## another visualisation, a close-up
plot(1:7, L$L[1,1:7], type = "l", 
     xlab = "q", ylab = "Explained variance")
for(i in 2:7) {lines(1:7, L$L[i,1:7], col = i)}
}

\keyword{EMMA}
