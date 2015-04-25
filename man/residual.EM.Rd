\name{residual.EM}
\alias{residual.EM}
\title{Function to evaluate residual end-member loading.}
\description{This function calculates an optional residual end-member loading. It uses 
the modelled end-member loadings as input and evaluates the root of 1 minus 
the sum of all squared loadings to analyse the remaining variance, e.g. 
if not all (robust) EMs are included (cf. Dietze et al., 2012). Negative 
values are set to zero.}
\usage{residual.EM(Vqn)}
\arguments{
  \item{Vqn}{Numeric matrix with m robust end-member loadings.}
}

\value{Numeric vector with residual end-member loading.}
\references{Dietze E, Hartmann K, Diekmann B, IJmker J, Lehmkuhl F, Opitz S, 
Stauch G, Wuennemann B, Borchers A. 2012. An end-member algorithm for 
deciphering modern detrital processes from lake sediments of Lake Donggi 
Cona, NE Tibetan Plateau, China. Sedimentary Geology 243-244: 169-180.}
\author{Michael Dietze, Elisabeth Dietze}



\seealso{\code{\link{EMMA}}, \code{\link{robust.EM}}}
\examples{
## Some preparing steps to retrieve only robust end-members
## load example data
data(X.artificial, envir = environment())

## truncate the data set for faster computation
X.trunc <- X.artificial[1:20,]

## define limits data set
limits = cbind(c(11, 31, 60, 78), 
               c(13, 33, 62, 80))

## define end-member numbers to test
q  <- 4:8
## define weight transformation limits
lw <- seq(from = 0, to = 0.1, by = 0.05)
## perform the robustness test
TR  <- test.robustness(X.trunc, q, lw, c = 100)

## extract robust end-members
REM <- robust.EM(Vqsn = TR$Vqsn, limits = limits, Vqn = TR$Vqn)

## Calculation of residual end-member loading
## define mean robust end-member loadings
Vqn.rob <- REM$Vqn.mean
## perform residual end-member loading calculation
Vqn.res <- residual.EM(Vqn.rob)

# Visualisation of the result
plot(NA, xlim = c(1, 80), ylim = c(0, 1))
for(i in 1:4) {lines(Vqn.rob[i,])}
lines(Vqn.res, col = 2)
}

\keyword{EMMA}
