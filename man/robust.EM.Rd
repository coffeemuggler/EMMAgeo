\name{robust.EM}
\alias{robust.EM}
\title{Function to extract robust end-members.}
\description{This function takes a matrix with end-member loadings and extracts those
whose modes fall into specified limits. The function returns a list with 
all passing end-member loadings and scores, along with their respective 
coumn-wise (variable-wise) measures of centrality and dispersion.}
\usage{robust.EM(Vqsn, limits, quantiles, Vqn, classunits, ID, plot = FALSE, 
    legend, ..., pm = FALSE)}
\arguments{
  \item{Vqsn}{Numeric matrix with m samples (rows) and n variables (columns).}
  \item{limits}{Numeric matrix with two columns that contain the boundaries of mode classes
for each end-member. The first column contains the lower, the second column
the upper limit. If \code{classunits} are provided, the limits are assumed 
to relate to these units, if omitted column-numbers of \code{Vqsn} 
are used.}
  \item{quantiles}{Optional numeric vector of length two with the quantiles to be
evaluated for the robust end-member loadings; default is 
\code{c(0.25, 0.75)}.}
  \item{Vqn}{Numeric matrix with optional normalised factor loadings. If present, the 
same factor loadings as the respectively selected end-member loadings 
are returned.}
  \item{classunits}{Numeric vector, optional class units (e.g. phi classes or micrometers) of 
the same length as columns of X.}
  \item{ID}{Numeric or character vector, optional sample IDs of the same
length as columns of X.}
  \item{plot}{Logical scalar, optional graphical output of the results, default is FALSE.
If set to TRUE, selected end-member loadings are plotted in different 
colours, according to the specified classes. All end-member loadings are 
plotted in pale colour, means and standard deviations are plotted above 
in thicker lines. To plot median and quantile range instead of mean and 
standard deviation, add \code{median = TRUE} as further plot parameter. 
See examples section for further advice.}
  \item{legend}{Character scalar, specifing legend position (cf. \code{\link{legend}}). If
omitted, no legend will be plotted, default is no legend.}
  \item{\dots}{Additional arguments passed to the plot function. Use \code{colour} instead
of \code{col} to create different colours.}
  \item{pm}{Logical scalar to enable pm.}
}

\value{A list object containing:
\item{Vqsn.data}{A list with Vqsn values.}
\item{Vqsn.mean}{A matrix with Vqsn means.}
\item{Vqsn.median}{A matrix with Vqsn medians.}
\item{Vqsn.sd}{A matrix with Vqsn standard deviations.}
\item{Vqsn.qt1}{A matrix with Vqsn quantiles 1.}
\item{Vqsn.qt2}{A matrix with Vqsn quantiles 2.}
\item{Vqn.data}{A list with Vqn values.}
\item{Vqn.mean}{A matrix with Vqn means.}
\item{Vqn.median}{A matrix with Vqn medians.}
\item{Vqn.sd}{A matrix with Vqn standard deviations.}
\item{Vqn.qt1}{A matrix with Vqn quantiles 1.}
\item{Vqn.qt2}{A matrix with Vqn quantiles 2.}}
\references{Dietze E, Hartmann K, Diekmann B, IJmker J, Lehmkuhl F, Opitz S, 
Stauch G, Wuennemann B, Borchers A. 2012. An end-member algorithm for 
deciphering modern detrital processes from lake sediments of Lake Donggi 
Cona, NE Tibetan Plateau, China. Sedimentary Geology 243-244: 169-180.}
\author{Michael Dietze, Elisabeth Dietze}



\seealso{\code{\link{EMMA}}, 
\code{\link{test.robustness}}, 
\code{\link{define.limits}}}
\examples{
## load example data, i.e. here TR
data(TR, envir = environment())

## define end-member limits
limits = cbind(c(11, 31, 60, 78), 
               c(13, 33, 62, 80))

## extract robust end-members with limits matrix
REM <- robust.EM(Vqsn = TR$Vqsn, limits = limits,
                 plot = TRUE,
                 legend = "topleft", 
                 cex = 0.7, 
                 colour = c("orange", "navyblue", "springgreen4", "red4"),
                 median = TRUE)
}

\keyword{EMMA}
