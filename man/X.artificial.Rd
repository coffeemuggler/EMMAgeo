\name{X.artificial}
\alias{X.artificial}
\docType{data}
\title{example data}
\description{Artificial data set created by randomly mixed normal distributions}
\usage{data(X.artificial)}
\format{
  The format is:
 num [1:80, 1:80] 0.00032 0.00057 0.00037 0.00029 ...
}
\details{The dataset is the result of five randomly mixed normal distributions,
         forming 80 samples, each represented by 80 classes. Some white noise
         is added.}
\examples{
## load example data set
data(X.artificial)

## plot first 10 samples stacked in one line plot
plot(NA, xlim = c(1, 80), ylim = c(1, 14))
for(i in 1:10) {lines(1:80, X.artificial[i,] * 50 + i)}
}
\keyword{datasets}
