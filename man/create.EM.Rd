% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/create.EM.R
\name{create.EM}
\alias{create.EM}
\title{Create grain-size-distributions.}
\usage{
create.EM(p1, p2, s, boundaries, n)
}
\arguments{
\item{p1}{\code{Numeric} vector, means of normal distributions, i.e. mode
positions.}

\item{p2}{\code{Numeric} vector, standard deviations of normal 
distributions, i.e. mode width.}

\item{s}{\code{Numeric} vector, relative proportions of each mode, i.e.
relative mode height.}

\item{boundaries}{\code{Numeric} vector of length two with class boundaries 
(i.e. \code{c(lower boundary, upper boundary)}).}

\item{n}{\code{Numeric} scalar with number of classes, i.e. resolution of 
the end-member.}
}
\value{
\code{Numeric} vector with normalised end-member loadings, 
consisting of the mixed normal distributions according to the input 
parameters.
}
\description{
This function allows creating artificial grain-size end-members. One such
"artificial end-member loading" may be composed of one or more superimposed
normal distributions.
}
\details{
When building a data set of many artificial end member loadings, these 
should all have the same \code{boundaries} and \code{n}. The function 
builds composites of individual normal distributions. Each distribution is 
scaled according to \code{s}. Finally the distribution is scaled to 100 \%.
}
\examples{

## set lower and upper class boundary, number of classes and class units
boundaries <- c(0, 11)
n <- 40
phi <- seq(from = boundaries[1], 
           to = boundaries[2], 
           length.out = n)

## create two artificial end-member loadings
EMa.1 <- create.EM(p1 = c(2, 5), p2 = c(1, 0.8), s = c(0.7, 0.3), 
                   boundaries = boundaries, n = n)
EMa.2 <- create.EM(p1 = c(4, 7), p2 = c(1.1, 1.4), s = c(0.5, 0.5),
                   boundaries = boundaries, n = n)

## plot the two artificial end-member loadings
plot(phi, EMa.1, type = "l")
lines(phi, EMa.2, col = "red")

}
\seealso{
\code{mix.EM}
}
\author{
Michael Dietze, Elisabeth Dietze
}
\keyword{EMMA}
