\name{create.sample}
\alias{create.sample}
\title{Function to create grain-size distributions}
\description{The function creates discrete samples of grain-size distributions from
continuous particle size measurements, e.g. obtained by SEM-imagery
(Francus, 1998), by moving a window of user-defined size with specified
step-sizes over the diameter locations.}
\usage{create.sample(data, window.width, step.size, class.limits)}
\arguments{
  \item{data}{Numeric matrix or data frame with two columns: depth (y-coordinates, in 
metric units) and diameter of the particles (either metric units or 
phi-scale).}
  \item{window.width}{Numeric scalar, vertical window width. All particles whose y-coordinates
fall into this moving window are sampled to build the grain-size 
distribution for this window.}
  \item{step.size}{Numeric scalar, amount of vertical distance between each sample, i.e.
moving window mid point.}
  \item{class.limits}{Numeric vector, specifying the class limits of the grain-size ditribution
(either in metric units or phi-scale). If missing, the limits will be set
automatically based on the diameter range of the data, with 20 classes, 
i.e. 21 limits.}
}

\value{A list with function output
\item{sample}{Numeric matrix, sample grain-size distributions.}
\item{n}{Numeric vector, number of speciment per sample.
}}
\references{Francus P. 1998. An image-analysis technique to measure
grain-size variation in thin sections of soft clastic sediments.
Sedimentary Geology 121(3-4): 289-298.}
\author{Michael Dietze, Elisabeth Dietze}



\seealso{\code{\link{interpolate.classes}}}
\examples{
## load example data set and create grain-size class vector
data(SEM.data, envir = environment())
phi <- seq(from = 3, 9, by = 0.5)

## create samples each 100 micrometers with a 500 micrometers wide window
sample <- create.sample(data = SEM.data,
                        window.width = 200,
                        step.size = 10,
                        class.limits = phi)
X <- sample$sample
n <- sample$n

## plot grain-size distributions colorised by sample density
plot(NA, xlim = range(phi), ylim = range(X), 
     main = "Grain-size distributions",
     xlab = "Grain-size class", ylab = "Relative amount [\%]")
for(i in 1:nrow(X)) {lines(phi, X[i,], col = rgb(0, 0, 1, n[i] / max(n)))}
}

\keyword{EMMA}
