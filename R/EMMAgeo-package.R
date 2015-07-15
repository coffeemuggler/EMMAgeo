

#' End-member modelling algorithm and supporting functions for grain-size
#' analysis
#' 
#' This package provides a set of functions for end-member modelling
#' analysis of grain-size data (EMMAgeo).
#' 
#' \tabular{ll}{ Package: \tab EMMAgeo\cr Type: \tab Package\cr Version: \tab
#' 0.9.2\cr Date: \tab 2015-06-07\cr License: \tab GPL-3\cr }
#' 
#' @name EMMAgeo-package
#' @aliases EMMAgeo
#' @docType package
#' @author Michael Dietze, Elisabeth Dietze
#' @keywords package
#' @import GPArotation
#' @import limSolve
#' @import shape
#' @import shiny
NULL

#' example data
#' 
#' Robust end-members, a list with output of the function robust.EM()
#' 
#' The dataset is the result of the function robust.EM() of the R-package
#' EMMAgeo.
#' 
#' @name REM
#' @docType data
#' @format The format is: List of 12 $ Vqsn.data :List of 4 ..$ : num [1:15,
#' 1:80] 0.18929 0.184 0.18304 0.00698 0.02033 ...
#' @keywords datasets
#' @examples
#' 
#' ## load example data set
#' data(examples)
#' 
NULL





#' example data
#' 
#' Depth and diameter of SEM-imaged grains
#' 
#' The dataset is the result of a classified and digitised SEM image and
#' contains two columns: depth and diameter, both in micrometers.
#' 
#' @name SEM.data
#' @docType data
#' @format The format is: num [1:1237, 1:2] 5.45 2.59 2.77 4.43 4.85 ...
#' @keywords datasets
#' @examples
#' 
#' ## load example data set
#' data(examples)
#' 
NULL





#' example data
#' 
#' A list with output of the function test.robustness()
#' 
#' The dataset is the result of the function test.robustness() of the R-package
#' EMMAgeo.
#' 
#' @name TR
#' @docType data
#' @format The format is: List of 8 $ q : num [1:90] 4 4 4 4 4 4 4 4 4 4 ...  $
#' lw : num [1:90] 0 0 0 0 0.05 0.05 0.05 0.05 0.1 0.1 ...  $ modes : num
#' [1:90] 12 32 61 80 12 32 61 80 12 32 ...
#' @keywords datasets
#' @examples
#' 
#' ## load example data set
#' data(examples)
#' 
NULL





#' example data
#' 
#' Artificial data set created by randomly mixed normal distributions
#' 
#' The dataset is the result of five randomly mixed normal distributions,
#' forming 80 samples, each represented by 80 classes. Some white noise is
#' added.
#' 
#' @name X.artificial
#' @docType data
#' @format The format is: num [1:80, 1:80] 0.00032 0.00057 0.00037 0.00029 ...
#' @keywords datasets
#' @examples
#' 
#' ## load example data set
#' data(X.artificial)
#' 
#' ## plot first 10 samples stacked in one line plot
#' plot(NA, xlim = c(1, 80), ylim = c(1, 14))
#' for(i in 1:10) {lines(1:80, X.artificial[i,] * 50 + i)}
#' 
NULL





#' example data
#' 
#' This file is simply used to securely locate the package path
#' 
#' No details needed
#' 
#' @name X.to.find.EMMAgeo.path.MDT
#' @docType data
#' @format The format is: List of 12 $ Vqsn.data :List of 4 ..$ : num [1:15,
#' 1:80] 0.18929 0.184 0.18304 0.00698 0.02033 ...
#' @keywords datasets
#' @examples
#' 
#' ## load example data set
#' data(examples)
#' 
NULL



