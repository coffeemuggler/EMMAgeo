#' Calculate the initial cumulative explained variance of factors.
#' 
#' This function performs eigenspace decomposition using the weight-transformed
#' matrix W to determine the explained variance with increasing number of 
#' factors. Depending on the number of provided weight transformation limits 
#' (\code{l}) a vector or a matrix is returned.
#' 
#' The results may be used to define a minimum number of end-members for
#' subsequent modelling steps, e.g. by using the Kaiser criterion, which
#' demands a minimum number of eigenvalues to reach a squared R of 0.95.
#' 
#' @param X \code{Numeric} matrix, input data set with m samples (rows) 
#' and n variables (columns).
#' 
#' @param l \code{Numeric} vector, weight tranformation limits, i.e.
#' quantiles; default is 0.
#' 
#' @param c \code{Numeric} scalar, constant sum scaling parameter, e.g.
#' 1, 100, 1000; default is 100.
#' 
#' @param r.min \code{Numeric} scalar, minimum value of explained variance to
#' be reached by the end-members included, default is 0.95.
#' 
#' @param plot \code{Logical} scalar, optional graphical output of the results,
#' default is FALSE.
#' 
#' @param legend \code{Character} scalar, specify legend position (cf.
#' \code{\link{legend}}). If omitted, no legend will be plotted, default is no
#' legend.
#' 
#' @param \dots Additional arguments passed to the plot function. Use
#' \code{colour} instead of \code{col} to create different colours.
#' 
#' @return \code{List} with objects \item{L}{Vector or matrix of cumulative
#' explained variance.} \item{q.min}{Vector with number of factors that passed
#' r.min. }
#' 
#' @author Michael Dietze, Elisabeth Dietze
#' @references Dietze E, Hartmann K, Diekmann B, IJmker J, Lehmkuhl F, Opitz S,
#' Stauch G, Wuennemann B, Borchers A. 2012. An end-member algorithm for
#' deciphering modern detrital processes from lake sediments of Lake Donggi
#' Cona, NE Tibetan Plateau, China. Sedimentary Geology 243-244: 169-180.
#' @keywords EMMA
#' @examples
#' 
#' ## load example data set
#' data(example_X)
#' 
#' ## create sequence of weight transformation limits
#' l <- seq(from = 0, to = 0.2, 0.02)
#' 
#' ## perform the test and show q.min
#' L <- test.factors(X = X, l = l, c = 100, plot = TRUE)
#' L$q.min
#' 
#' ## a visualisation with more plot parameters
#' L <- test.factors(X = X, l = l, c = 100, plot = TRUE, 
#'                   ylim = c(0.5, 1), xlim = c(1, 7), 
#'                   legend = "bottomright", cex = 0.7)
#' 
#' ## another visualisation, a close-up
#' plot(1:7, L$L[1,1:7], type = "l", 
#'      xlab = "q", ylab = "Explained variance")
#' for(i in 2:7) {lines(1:7, L$L[i,1:7], col = i)}
#' 
#' @export test.factors
test.factors <- function(
  X,
  l,
  c,
  r.min = 0.95,
  plot = FALSE,
  legend,
  ...
){
  
  ## check for l vs. lw
  if("lw" %in% names(list(...))) {
    stop('Parameter "lw" is depreciated. Use "l" instead.')
  }
  
  ## check/set default values
  if(missing(l) == TRUE) {l <- 0}
  if(missing(c) == TRUE) {c <- 100}
  
  ## rescale X to constant sum
  X <- X / apply(X, 1, sum) * c
  
  ## create result matrix
  L <- matrix(nrow = length(l), ncol = ncol(X))

  ## loop through all l values
  for (i in 1:length(l)) {

    ## calculate weight limit quantiles column-wise
    ls <- sapply(X = 1:ncol(X), FUN = function(j) {
      quantile(x = X[,j], probs = c(l[i], 1 - l[i]), type = 5)})

    ## perform weight-transformation
    W <- t((t(X) - ls[1,]) / (ls[2,] - ls[1,]))

    ## create similarity matrix as outer product
    A <- t(W) %*% W

    ## perform eigenspace decomposition
    EIG <- eigen(A)

    ## assign raw eigenvectors V and eigenvalues L
    V <- EIG$vectors[,order(seq(ncol(A), 1, -1))]
    L.raw <- EIG$values[order(seq(ncol(A), 1, -1))]

    ## calculate cumulative sums of eigenvalues
    L[i,] <- cumsum(sort(L.raw / sum(L.raw), decreasing = TRUE))
  }
  
  ## define result vector of minumum number of end-members
  q.min <- rep(NA, length(l))
  
  ## find minimum number of end-members to pass threshold value r.min
  for (i in 1:length(l)) {q.min[i] <- seq(1, ncol(X))[L[i,] >= r.min][1]}
  
  ## optional plot of the results (explained variance vs. number of factors)
  if(plot == TRUE) {
    ## adjust plot margins
    par(oma = c(0, 1, 0, 0))
    ## read additional arguments list and check/set default values
    extraArgs <- list(...)

    main <- if("main" %in% names(extraArgs)) {
      extraArgs$main
    } else {
      "Variance explained by factors"
    }
    
    xlab <- if("xlab" %in% names(extraArgs)) {
      extraArgs$xlab
    } else {
      "Number of factors (q)"
    }
    
    ylab <- if("ylab" %in% names(extraArgs)) {
      extraArgs$ylab
    } else {
      expression(paste("Explained variance (", R^2, ")", sep = ""))
    }
    
    xlim <- if("xlim" %in% names(extraArgs)) {
      extraArgs$xlim
    } else {
      c(1, ncol(L))
    }
    
    ylim <- if("ylim" %in% names(extraArgs)) {
      extraArgs$ylim
    } else {
      range(L)
    }
    
    legend.text <- if("legend" %in% names(extraArgs)) {
      extraArgs$legend
    } else {
        round(l, 3)
    }
    
    legend.cex <- if("cex" %in% names(extraArgs)) {
      extraArgs$cex
    } else {
      1
    }
    
    legend.lty <- if("lty" %in% names(extraArgs)) {
      extraArgs$lty
    } else {
        1
    }
    
    legend.title <- if("title" %in% names(extraArgs)) {
      extraArgs$title
    } else {
      "l"
    }
    
    col <- if("col" %in% names(extraArgs)) {
      extraArgs$col
    } else {
      seq(1, length(l))
    }
        
    ## plot data
    plot(NA, 
         main = main,
         xlab = xlab, 
         ylab = ylab,
         xlim = xlim,
         ylim = ylim)
    
    for(i in 1:length(l)) {
      
      lines(x = 1:ncol(L), y = L[i,], col = i)
      points(x = 1:ncol(L), y = L[i,], col = i)
    }
    
    ## optionally add legend
    if(missing(legend) == FALSE) {
      legend.position <- legend
      legend(x = legend.position,
             legend = legend.text,
             cex = legend.cex,
             lty = legend.lty,
             title = legend.title,
             bty = "n",
             col = if("col" %in% names(extraArgs)) {
               extraArgs$col
             } else {
               seq(1, length(l))
             })
    }
  }
  
  ## return output
  return(list(L = L,
              q.min = q.min))
}
