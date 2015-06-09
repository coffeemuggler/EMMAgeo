\name{test.lw.max}
\alias{test.lw.max}
\title{Function to find maximum possible lw value.}
\description{This function approximates the highest possible value for lw in a 
nested loop.}
\usage{test.lw.max(X, n = 10)}
\arguments{
  \item{X}{Numeric matrix with m samples (rows) and n variables (columns).}
  \item{n}{Numeric scalar, number of loop runs and values per loop.}
}

\value{A numeric scalar, maximal possible lw value.}
\references{Dietze E, Hartmann K, Diekmann B, IJmker J, Lehmkuhl F, Opitz S, 
Stauch G, Wuennemann B, Borchers A. 2012. An end-member algorithm for 
deciphering modern detrital processes from lake sediments of Lake Donggi 
Cona, NE Tibetan Plateau, China. Sedimentary Geology 243-244: 169-180. \cr
Klovan JE, Imbrie J. 1971. An Algorithm and FORTRAN-IV Program for 
Large-Scale Q-Mode Factor Analysis and Calculation of Factor Scores. 
Mathematical Geology 3: 61-77.}
\author{Michael Dietze, Elisabeth Dietze}



\seealso{\code{\link{EMMA}}, \code{\link{check.data}}, 
\code{\link{test.parameters}}}
\examples{
## load example data set
data(X.artificial, envir = environment())

## create weight transformation limits vector
lw <- seq(from = 0, to = 0.6, by = 0.02)

## test the vector
test.lw(X = X.artificial, lw = lw)
}

\keyword{EMMA}
