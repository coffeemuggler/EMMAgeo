\name{interpolate.classes}
\alias{interpolate.classes}
\title{Function to interpolate classes.}
\description{This function interpolates grain-size classes, either to higher or to lower
resolution.}
\usage{interpolate.classes(X, boundaries.in, boundaries.out, method = "natural", 
    fixed.start = TRUE)}
\arguments{
  \item{X}{Numeric matrix with m samples (rows) and n variables (columns).}
  \item{boundaries.in}{Numeric vector with class boundaries of the input data.}
  \item{boundaries.out}{Numeric vector with class boundaries of the output data.}
  \item{method}{Logical scalar, interpolation method, one out of "linear" (linear 
interpolation), "fmm" (cubic spline), "natural" (natural spline), 
"periodic" (periodic spline). Default is \code{"natural"}.}
  \item{fixed.start}{Logical scalar, specifying if the outer boundaries should be set 
to the same values as in the original matrix, default is \code{TRUE}.
This may become necessary to avoid interpolation errors, see example.}
}

\value{Numeric matrix with interpolated class values.}
\references{Dietze E, Hartmann K, Diekmann B, IJmker J, Lehmkuhl F, Opitz S, 
Stauch G, Wuennemann B, Borchers A. 2012. An end-member algorithm for 
deciphering modern detrital processes from lake sediments of Lake Donggi 
Cona, NE Tibetan Plateau, China. Sedimentary Geology 243-244: 169-180.}
\author{Michael Dietze, Elisabeth Dietze}



\seealso{\code{\link{EMMA}}, \code{\link{approx}}, \code{\link{spline}}}
\examples{
## load example data
data(X.artificial, envir = environment())
classes.in <- seq(from = 1, to = 10, length.out = ncol(X.artificial))
  
## Example 1 - decrease the class numbers
## define number of output classes
classes.out <- seq(1, 10, length.out = 20)

## interpolate the data set
Y <- interpolate.classes(X = X.artificial, 
                         boundaries.in = classes.in, 
                         boundaries.out = classes.out,
                         method = "linear")

## show original vs. interpolation for first 10 samples
plot(NA, xlim = c(1, 10), ylim = c(0, 11))
for(i in 1:10) {
  lines(classes.in, X.artificial[i,] * 20 + i)
  lines(classes.out, Y[i,] * 20 + i, col = 2)
}

## Example 2 - increase the class numbers
## define number of output classes
classes.out <- seq(1, 10, length.out = 200)

## interpolate the data set
Y <- interpolate.classes(X = X.artificial, 
                         boundaries.in = classes.in, 
                         boundaries.out = classes.out)

## show original vs. interpolation for first 10 samples
plot(NA, xlim = c(1, 10), ylim = c(0, 11))
for(i in 1:10) {
  lines(classes.in, X.artificial[i,] * 20 + i)
  lines(classes.out, Y[i,] * 20 + i, col = 2)
}
}

\keyword{EMMA}
