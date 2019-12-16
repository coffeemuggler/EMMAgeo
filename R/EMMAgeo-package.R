

#' End-member modelling algorithm and supporting functions for unmixing 
#' grain-size distributions and further compositional data.
#' 
#' EMMAgeo provides a set of functions for end-member modelling analysis
#' (EMMA) of grain-size data and other cases of compositional data. EMMA 
#' describes a multivariate data set of m samples, each comprising n parameters 
#' (e.g. grain-size classes), as a linear combination of end-member loadings 
#' (the underlying distributions) and end-member scores (the contribution of
#' each loading to each sample).\cr EMMA can be run in two principal ways, 
#' a deterministic and a robust, including modelling the uncertainties. The 
#' deterministic way can be accessed simply with the function \code{EMMA()}. 
#' For the robust way there are two protocols that need to be respected. There 
#' is a compact protocol, which is mainly automated but needs adjustments by 
#' the user, and there is an extended protocol, which allows access to all 
#' parameterisation steps of robust EMMA.\cr The package contains further 
#' auxiliary functions to check and prepare input data, test parameters and
#' use a graphic user interface for deterministic EMMA. The package also 
#' contains an example data set, comprising meaured grain-size distributions 
#' of real world sediment end-members.
#' 
#' \tabular{ll}{ Package: \tab EMMAgeo\cr Type: \tab Package\cr Version: \tab
#' 0.9.7\cr Date: \tab 2019-05-10\cr License: \tab GPL-3\cr }
#' 
#' @name EMMAgeo-package
#' @aliases EMMAgeo
#' @docType package
#' @author Michael Dietze, Elisabeth Dietze
#' @keywords package
#' @import GPArotation
#' @import shiny
#' @importFrom grDevices adjustcolor col2rgb hsv rainbow rgb2hsv
#' @importFrom graphics abline axis barplot contour hist image 
#' layout legend lines locator mtext par plot points rug 
#' text polygon box segments
#' @importFrom stats median approx complete.cases cor density na.omit 
#' quantile rnorm runif sd spline var
#' @importFrom utils setTxtProgressBar txtProgressBar
#' @importFrom limSolve nnls
#' @importFrom caTools runmean
#' @importFrom matrixStats rowMins rowMaxs colVars
NULL

#' example data
#' 
#' Robust end-members, a list with output of the function robust.EM()
#' 
#' The dataset is the result of the function robust.EM() of the R-package
#' EMMAgeo.
#' 
#' @name EMrob
#' @docType data
#' @format The format is: List of 12 $ Vqsn.data :List of 4 ..$ : num [1:15,
#' 1:80] 0.18929 0.184 0.18304 0.00698 0.02033 ...
#' @keywords datasets
#' @examples
#' 
#' ## load example data set
#' data(example_EMrob)
#' 
NULL





#' example data
#' 
#' A list with output of the function test.robustness()
#' 
#' The dataset is the result of the function test.robustness() of the R-package
#' EMMAgeo.
#' 
#' @name EMpot
#' @docType data
#' @format The format is: List of 8 $ q : num [1:90] 4 4 4 4 4 4 4 4 4 4 ...  $
#' lw : num [1:90] 0 0 0 0 0.05 0.05 0.05 0.05 0.1 0.1 ...  $ modes : num
#' [1:90] 12 32 61 80 12 32 61 80 12 32 ...
#' @keywords datasets
#' @examples
#' 
#' ## load example data set
#' data(example_EMpot)
#' 
NULL





#' example data
#' 
#' Synthetic data set created by randomly mixed natural end-members
#' 
#' The dataset is the result of four mixed natural end-members.
#' 
#' @name X
#' @docType data
#' @format num [1:100, 1:116] 0.000899 0.000516 0.00136 0.000989 0.00102 ...
#' @keywords datasets
#' @examples
#' 
#' ## load example data set
#' data(example_X)
#' 
#' ## extract grain-size classes
#' s <- as.numeric(colnames(X))
#' 
#' ## plot first 10 samples stacked in one line plot
#' plot(NA, 
#'      xlim = c(1, ncol(X)), 
#'      ylim = c(1, 20))
#'      
#' for(i in 1:10) {
#'   lines(x = s, 
#'         y = X[i,] + i)
#' }
#' 
#' ## plot grain-size map
#' image(x = s, 
#'       z = t(X), 
#'       log = "x", 
#'       col = rainbow(n = 250))
#' 
NULL




