% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/test.robustness.R
\name{test.robustness}
\alias{test.robustness}
\title{Test model robustness.}
\usage{
test.robustness(
  X,
  q,
  l,
  P,
  c,
  classunits,
  ID,
  rotation = "Varimax",
  ol.rej,
  mRt.rej,
  plot = FALSE,
  ...
)
}
\arguments{
\item{X}{Numeric matrix with m samples (rows) and n variables (columns).}

\item{q}{Numeric vector with number of end-members to be modelled.}

\item{l}{Numeric vector specifying the weight tranformation limits, i.e.
quantiles; default is 0.}

\item{P}{Numeric matrix, optional alternative input parameters for q and l,
either of the form m:3 with m variations in the columns q.min, q.max, l or
of the form m:2 with m variations in the columns q, l.}

\item{c}{Numeric scalar specifying the constant sum scaling parameter, e.g.
1, 100, 1000; default is 100.}

\item{classunits}{Numeric vector, optional class units (e.g. phi classes or
micrometers) of the same length as columns of X.}

\item{ID}{Numeric or character vector, optional sample IDs of the same
length as columns of X.}

\item{rotation}{Character scalar, rotation type, default is "Varimax" (cf.
Dietze et al., 2012). One out of the rotations provided in GPArotation is
possible.}

\item{ol.rej}{Numeric scalar, optional rejection threshold for overlapping
criterion.  All model runs with overlapping end-members greater than the
specified integer will be removed.}

\item{mRt.rej}{Numeric scalar, optional rejection threshold for mean total
explained variance criterion. All modelled end-members below the specified
value will be removed.}

\item{plot}{Logical scalar, optional graphical output of the results,
default is FALSE. If set to TRUE, end-member loadings and end-member scores
are plotted.}

\item{\dots}{Additional arguments passed to the plot function (see details).}
}
\value{
A list with objects \item{q}{Vector with q.} \item{l}{Vector with
l.} \item{modes}{Vector with mode class.} \item{mRt}{Vector with mean total
explained variance.} \item{ol}{Vector with n overlapping end-members.}
\item{loadings}{Matrix with normalised rescaled end-member loadings.}
\item{Vqsn}{Matrix with rescaled end-member loadings.} \item{Vqn}{Matrix
with normalised factor loadings.}
}
\description{
This function takes a definition of weight transformation 
limits and corresponding minimum and maximum numbers of end-members to 
model all end-member scenarios in accordance with these parameters. Based 
on the output the user can decide on robust end-members.
}
\details{
The function value \code{$loadings} is redundant but was added for user
convenience.\cr Since the
function returns two plots, additional graphical parameters must be
specified as vector with the first element for the first plot and the second
element for the second plot. If graphical parameters are natively vectors
(e.g. a sequence of colours), they must be specified as matrices with each
vector as a row. If colours are specified, \code{colour} should be used
instead of \code{col}. \code{ylim} can only be modified for the first plot.
See example section for further advice.
}
\examples{

## load example data set
data(example_X)

## Example 1 - perform the most simple test
q  <- 4:7
l <- seq(from = 0, to = 0.1, by = 0.05)

M1  <- test.robustness(X = X, q = q, l = l, 
                       ol.rej = 1, mRt.rej = 0.8, 
                       plot = TRUE,
                       colour = c(4, 7),
                       xlab = c(expression(paste("Grain size (", phi, ")", 
                                                 sep = "")), 
                                expression(paste("Grain size (", phi, ")", 
                                                 sep = ""))))

## Example 2 -  perform the test without rejection criteria and plots
P  <- cbind(rep(q[1], length(l)),
            rep(q[3], length(l)),
            l)
M2  <- test.robustness(X = X, P = P)

## Plot 1 - end-member loadings which do not overlap and yielded mRt > 0.80.
plot(M2$Vqsn[1,], type = "l", ylim = c(0, max(M2$Vqsn, na.rm = TRUE)),
     main = "End-member loadings")
  for (i in 2:nrow(M2$Vqsn)) lines(M2$Vqsn[i,])

# Plot 2 - histogram of mode positions
hist(M2$modes,
     breaks = 1:ncol(X), 
     main = "Mode positions",
     xlab = "Class")

# Plot 3 - positions of modelled end-member modes by number of end-members
# Note how scatter in end-member position decreases for the "correct" number 
# of modelled end-members (6) and an appropriate weight limit (ca. 0.1).
ii <- order(M2$q, M2$modes)
modes <- t(rbind(M2$modes, M2$q))[ii,]
plot(modes[,1],
     seq(1, nrow(modes)), 
     main = "Model overview",
     xlab = "Class", 
     ylab = "EM number in model run", 
     pch = as.character(modes[,2]), 
     cex = 0.7)

# Illustrate mode positions as stem-and-leave-plot, useful as a simple
# check, which mode maxima are consistently fall into which grain-size 
# class (useful to define "limits" in robust.EM).
stem(M2$modes, scale = 2)


}
\references{
Dietze E, Hartmann K, Diekmann B, IJmker J, Lehmkuhl F, Opitz S,
Stauch G, Wuennemann B, Borchers A. 2012. An end-member algorithm for
deciphering modern detrital processes from lake sediments of Lake Donggi
Cona, NE Tibetan Plateau, China. Sedimentary Geology 243-244: 169-180.
}
\author{
Michael Dietze, Elisabeth Dietze
}
\keyword{EMMA}
