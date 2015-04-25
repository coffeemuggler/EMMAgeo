\name{convert.units}
\alias{convert.units}
\title{Function to convert between phi and micrometers.}
\description{The function converts values from the phi-scale to the 
micrometer-scale and vice versa.}
\usage{convert.units(phi, mu)}
\arguments{
  \item{phi}{Numeric vector with grain-size class values in phi to be converted.}
  \item{mu}{Numeric vector with grain-size class values in micrometres to be converted.}
}

\value{Numeric vector with converted grain-size class values.}

\author{Michael Dietze, Elisabeth Dietze}



\seealso{\code{\link{interpolate.classes}}}
\examples{
## generate phi-values
phi <- -2:5

## convert and show phi to mu
mu  <- convert.units(phi = phi)
mu

## convert and show phi to mu
convert.units(mu = mu)
}

\keyword{EMMA}
