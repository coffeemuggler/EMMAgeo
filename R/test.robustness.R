#' Function to test model robustness.
#' 
#' All possible combinations of number of end-members and weight transformation
#' limits are used to perform EMMA. The resulting loadings are written to an
#' output matrix and mode positions (i.e. class with maximum loading) of all
#' loadings are evaluated and returned.
#' 
#' The function value \code{$loadings} is redundant but was added for user
#' convenience.
#' 
#' @param X Numeric matrix with m samples (rows) and n variables (columns).
#' @param q Numeric vector with number of end-members to be modelled.
#' @param lw Numeric vector specifying the weight tranformation limits, i.e.
#' quantiles; default is 0.
#' @param P Numeric matrix, optional alternative input parameters for q and lw,
#' either of the form m:3 with m variations in the columns q.min, q.max, lw or
#' of the form m:2 with m variations in the columns q, lw.
#' @param c Numeric scalar specifying the constant sum scaling parameter, e.g.
#' 1, 100, 1000; default is 100.
#' @param classunits Numeric vector, optional class units (e.g. phi classes or
#' micrometers) of the same length as columns of X.
#' @param ID Numeric or character vector, optional sample IDs of the same
#' length as columns of X.
#' @param rotation Character scalar, rotation type, default is "Varimax" (cf.
#' Dietze et al., 2012). One out of the rotations provided in GPArotation is
#' possible (cf. \code{\link{rotations}}).
#' @param ol.rej Numeric scalar, optional rejection threshold for overlapping
#' criterion.  All model runs with overlapping end-members greater than the
#' specified integer will be removed.
#' @param mRt.rej Numeric scalar, optional rejection threshold for mean total
#' explained variance criterion. All modelled end-members below the specified
#' value will be removed.
#' @param plot Logical scalar, optional graphical output of the results,
#' default is FALSE. If set to TRUE, end-member loadings and end-member scores
#' are plotted.
#' @param \dots Additional arguments passed to the plot function. Since the
#' function returns two plots, additional graphical parameters must be
#' specified as vector with the first element for the first plot and the second
#' element for the second plot. If graphical parameters are natively vectors
#' (e.g. a sequence of colours), they must be specified as matrices with each
#' vector as a row. If colours are specified, \code{colour} should be used
#' instead of \code{col}. \code{ylim} can only be modified for the first plot.
#' See example section for further advice.
#' @param pm Logical scalar to enable pm.
#' @return A list with objects \item{q}{Vector with q.} \item{lw}{Vector with
#' lw.} \item{modes}{Vector with mode class.} \item{mRt}{Vector with mean total
#' explained variance.} \item{ol}{Vector with n overlapping end-members.}
#' \item{loadings}{Matrix with normalised rescaled end-member loadings.}
#' \item{Vqsn}{Matrix with rescaled end-member loadings.} \item{Vqn}{Matrix
#' with normalised factor loadings.}
#' @author Michael Dietze, Elisabeth Dietze
#' @references Dietze E, Hartmann K, Diekmann B, IJmker J, Lehmkuhl F, Opitz S,
#' Stauch G, Wuennemann B, Borchers A. 2012. An end-member algorithm for
#' deciphering modern detrital processes from lake sediments of Lake Donggi
#' Cona, NE Tibetan Plateau, China. Sedimentary Geology 243-244: 169-180. \cr
#' @keywords EMMA
#' @examples
#' 
#' ## load example data set
#' data(X.artificial, envir = environment())
#' 
#' ## Example 1 - perform the most simple test
#' q  <- 4:7
#' lw <- seq(from = 0, to = 0.1, by = 0.02)
#' 
#' M1  <- test.robustness(X = X.artificial, q = q, lw = lw, 
#'                        ol.rej = 1, mRt.rej = 0.8, 
#'                        plot = TRUE,
#'                        colour = c(4, 7),
#'                        xlab = c(expression(paste("Grain size (", phi, ")", 
#'                                                  sep = "")), 
#'                                 expression(paste("Grain size (", phi, ")", 
#'                                                  sep = ""))))
#' 
#' ## Example 2 -  perform the test without rejection criteria and plots
#' P  <- cbind(rep(q[1], length(lw)),
#'             rep(q[3], length(lw)),
#'             lw)
#' M2  <- test.robustness(X = X.artificial, P = P)
#' 
#' ## Plot 1 - end-member loadings which do not overlap and yielded mRt > 0.80.
#' plot(M2$Vqsn[1,], type = "l", ylim = c(0, max(M2$Vqsn, na.rm = TRUE)),
#'      main = "End-member loadings")
#'   for (i in 2:nrow(M2$Vqsn)) lines(M2$Vqsn[i,])
#' 
#' # Plot 2 - histogram of mode positions
#' hist(M2$modes,
#'      breaks = 1:ncol(X.artificial), 
#'      main = "Mode positions",
#'      xlab = "Class")
#' 
#' # Plot 3 - positions of modelled end-member modes by number of end-members
#' # Note how scatter in end-member position decreases for the "correct" number 
#' # of modelled end-members (6) and an appropriate weight limit (ca. 0.1).
#' ii <- order(M2$q, M2$modes)
#' modes <- t(rbind(M2$modes, M2$q))[ii,]
#' plot(modes[,1],
#'      seq(1, nrow(modes)), 
#'      main = "Model overview",
#'      xlab = "Class", 
#'      ylab = "EM number in model run", 
#'      pch = as.character(modes[,2]), 
#'      cex = 0.7)
#' 
#' # Illustrate mode positions as stem-and-leave-plot, useful as a simple
#' # check, which mode maxima are consistently fall into which grain-size 
#' # class (useful to define "limits" in robust.EM).
#' stem(M2$modes, scale = 2)
#' 
#' 
#' @export test.robustness
test.robustness <-
structure(function # Function to test model robustness.
### All possible combinations of number of end-members and weight 
### transformation limits are used to perform EMMA. The resulting loadings are 
### written to an output matrix and mode positions (i.e. class with maximum 
### loading) of all loadings are evaluated and returned.
(X,
### Numeric matrix with m samples (rows) and n variables (columns).
q, 
### Numeric vector with number of end-members to be modelled.
lw,
### Numeric vector specifying the weight tranformation limits, i.e. 
### quantiles; default is 0.
P,
### Numeric matrix, optional alternative input parameters for q and lw, either
### of the form m:3 with m variations in the columns q.min, q.max, lw or of 
### the form m:2 with m variations in the columns q, lw.
c,
### Numeric scalar specifying the constant sum scaling parameter, e.g. 1, 
### 100, 1000; default is 100.
classunits,
### Numeric vector, optional class units (e.g. phi classes or micrometers) of 
### the same length as columns of X.
ID,
### Numeric or character vector, optional sample IDs of the same
### length as columns of X.
rotation = "Varimax",
### Character scalar, rotation type, default is "Varimax" (cf. Dietze et 
### al., 2012). One out of the rotations provided in GPArotation is 
### possible (cf. \code{\link{rotations}}).
ol.rej,
### Numeric scalar, optional rejection threshold for overlapping criterion. 
### All model runs with overlapping end-members greater than the specified 
### integer will be removed.
mRt.rej,
### Numeric scalar, optional rejection threshold for mean total explained 
### variance criterion. All modelled end-members below the specified value 
### will be removed.
plot = FALSE,
### Logical scalar, optional graphical output of the results, default is 
### FALSE. If set to TRUE, end-member loadings and end-member scores are 
### plotted.
...,
### Additional arguments passed to the plot function. Since the function 
### returns two plots, additional graphical parameters must be specified as
### vector with the first element for the first plot and the second element
### for the second plot. If graphical parameters are natively vectors (e.g.
### a sequence of colours), they must be specified as matrices with each 
### vector as a row. If colours are specified, \code{colour} 
### should be used instead of \code{col}. \code{ylim} can only be modified
### for the first plot. See example section for further advice.
pm = FALSE
### Logical scalar to enable pm.
){
  ## check/set class units vector and test for consistency
  if(missing(classunits) == TRUE) {classunits <- 1:ncol(X)}
  if(ncol(X) != length(classunits)) {stop(
    "Units vector is not of same length as variables.")}
  
  ## check/set constant sum value
  if(missing(c) == TRUE) {c <- 100}
  
  ## check/set ID vector and test for consistency
  if(missing(ID) == TRUE) ID <- 1:nrow(X)
  if(nrow(X) != length(ID)) stop(
    "ID vector is not of same length as variables.")
  
  ## create vectors with test values of q and lw
  if(missing(P) == TRUE) {
  ## option 1 - no input matrix P given
    ## create test vectors q and lw
    q.t  <- rep(q, each = length(lw))
    lw.t <- rep(seq(min(lw), max(lw), length.out = length(lw)), 
              length(q))
  } else if(ncol(P) == 3) {
  
  ## option 2 - P with q.min, q.max and lw
    ## create help and dummy variables
    N <- sum(P[,2] - P[,1]) # total numerb of q 
    q.t  <- NA # dummy vector q.t
    lw.t <- NA # dummy vector lw.t
    ## attach q and lw series to dummy vectors 
    for(i in 1:nrow(P)) {
      q.n  <- seq(P[i,1], P[i,2])
      lw.n <- rep(P[i,3], length(q.n))
      q.t  <- c(q.t, q.n)
      lw.t <- c(lw.t, lw.n)
    }
    ## remove dummy vector values
    q.t  <- q.t[2:length(q.t)]
    lw.t <- lw.t[2:length(lw.t)]
  } else if(ncol(P) == 2) {
    
  ## option 3 - P with q and lw
    ## assgin values to test vectors q and lw
    q.t  <- P[,1]
    lw.t <- P[,2]
  }
  
  ## create result matrices
  data.t <- matrix(nrow = sum(q.t), 
                   ncol = 5) # metadata (q, lw, modes, mRt, ol)
  Vqn.t  <- matrix(nrow = sum(q.t), 
                   ncol = ncol(X)) # end-member loadings
  Vqsn.t <- matrix(nrow = sum(q.t), 
                   ncol = ncol(X)) # normalised end-member loadings
  
  ## set counter variables
  ni <- 1
  nj <- q.t[1]
  i.pb <- 1
  
  ## loop through all parameter combinations
  for (i in 1:length(q.t)) {
    ## perform EMMA
    EM <- EMMA(X, 
               q = q.t[i], 
               lw = lw.t[i], 
               classunits = classunits, 
               ID = ID, 
               c = c, 
               rotation = rotation)

    ## assign result values to matrices
    Vqn.t[ni:nj,]     <- EM$Vqn
    Vqsn.t[ni:nj,]    <- EM$Vqsn
    data.t[ni:nj,]    <- c(rep(q.t[i], q.t[i]), # q
                           rep(lw.t[i], q.t[i]), # lw
                           rep(NA, q.t[i]), # dummy peak position
                           rep(mean(c(EM$Rm, EM$Rn), 
                                    na.rm = TRUE), q.t[i]), # mRt
                           rep(EM$ol, q.t[i])) # ol
    
    ## update counter variables
    ni <- ifelse(ni < sum(q.t), nj + 1, ni)
    nj <- ifelse(ni < sum(q.t), ni + q.t[i+1] - 1, ni)
  }
  
  ## determine mode position
  for(i in 1:nrow(Vqsn.t)) {
    data.t[i,3] <- classunits[Vqsn.t[i,1:ncol(
      Vqsn.t)] == max(Vqsn.t[i,1:ncol(Vqsn.t)])]
  }
  
  ## optionally remove all data sets that failed rejection criterion ol.rej
  if(missing(ol.rej) == FALSE) {
    ## identify rows that passed criterion
    ID     <- data.t[,5] < ol.rej
    ## keep data rows that passed criterion
    Vqsn.t <- Vqsn.t[ID,]
    Vqn.t  <- Vqn.t[ID,]
    data.t <- data.t[ID,]
    
    ## stop if result is NULL
    if(is.matrix(Vqsn.t) == FALSE) stop("No output passed threshold ol.rej.")
  }
  
  ## optionally remove all data sets that failed rejection criterion mRt.rej
  if(missing(mRt.rej) == FALSE) {
    ## identify rows that passed criterion
    ID     <- data.t[,4] > mRt.rej
    ## keep data rows that passed criterion
    Vqsn.t <- Vqsn.t[ID,]
    Vqn.t  <- Vqn.t[ID,]
    data.t <- data.t[ID,]
    ## stop if result is NULL
    if(is.matrix(Vqsn.t) == FALSE) stop("No output passed threshold mRt.rej.")
  }
  
  ## optionally plot resulting end-member loadings and a histogram
  if(plot == TRUE) {
    ## read out additional arguments list
    extraArgs <- list(...)
    colour <- if("colour" %in% names(extraArgs)) {extraArgs$colour} else
      {c("black", "grey")}
    main <- if("main" %in% names(extraArgs)) {extraArgs$main} else
      {c(expression(paste("Loadings (", V[qsn], ")", sep = "")), 
         "Mode positions")}
    xlab <- if("xlab" %in% names(extraArgs)) {extraArgs$xlab} else
      {c("Classes", "Classes")}
    ylab <- if("ylab" %in% names(extraArgs)) {extraArgs$ylab} else
      {c("Amount, relative", "Amount, relative")}
    ylim <- if("ylim" %in% names(extraArgs)) {extraArgs$ylim} else
      {c(0, max(Vqsn.t, na.rm = TRUE))}
    
    ## setup plot area
    par(mfrow = c(1, 2),
        oma = c(0, 1, 0, 0))
    
    ## plot end-member loadings
    plot(classunits,
         Vqsn.t[1,], type = "l", 
         main = main[1],
         xlab = xlab[1],
         ylab = ylab[1],
         ylim = ylim,
         col = colour[1])
    for (i in 2:nrow(Vqsn.t)) lines(classunits, Vqsn.t[i,], col = colour[1])
    
    ## plot histogram of mode positions
    hist(data.t[,3], 
         breaks = classunits, 
         main = main[2],
         xlab = xlab[2],
         ylab = ylab[2],
         axes = FALSE,
         col = colour[2])
    axis(side = 1)
    rug(data.t[,3])
    
    ## reset format of the plot area
    par(mfrow = c(1, 1),
        oma = c(0, 0, 0, 0))
  }
  
  ## optionally add pm
  if(pm == TRUE) {pm <- check.data(matrix(runif(4), ncol = 2), 
                                   5, 0.01, 100, invisible = FALSE)}
  
  ##value<< A list with objects
  list(q        = data.t[,1], ##<< Vector with q.
       lw       = data.t[,2], ##<< Vector with lw.
       modes    = data.t[,3], ##<< Vector with mode class.
       mRt      = data.t[,4], ##<< Vector with mean total explained variance.
       ol       = data.t[,5], ##<< Vector with n overlapping end-members.
       loadings = Vqsn.t,     ##<< Matrix with normalised rescaled end-member loadings.
       Vqsn     = Vqsn.t,     ##<< Matrix with rescaled end-member loadings.
       Vqn      = Vqn.t)      ##<< Matrix with normalised factor loadings.
  ##end<<

  ##details<<
  ## The function value \code{$loadings} is redundant but was added for
  ## user convenience.
  
  ##references<<
  ## Dietze E, Hartmann K, Diekmann B, IJmker J, Lehmkuhl F, Opitz S, 
  ## Stauch G, Wuennemann B, Borchers A. 2012. An end-member algorithm for 
  ## deciphering modern detrital processes from lake sediments of Lake Donggi 
  ## Cona, NE Tibetan Plateau, China. Sedimentary Geology 243-244: 169-180. \cr
  
  ##keyword<<
  ## EMMA
}, ex = function(){
  ## load example data set
  data(X.artificial, envir = environment())
  
  ## Example 1 - perform the most simple test
  q  <- 4:7
  lw <- seq(from = 0, to = 0.1, by = 0.02)
  
  M1  <- test.robustness(X = X.artificial, q = q, lw = lw, 
                         ol.rej = 1, mRt.rej = 0.8, 
                         plot = TRUE,
                         colour = c(4, 7),
                         xlab = c(expression(paste("Grain size (", phi, ")", 
                                                   sep = "")), 
                                  expression(paste("Grain size (", phi, ")", 
                                                   sep = ""))))
  
  ## Example 2 -  perform the test without rejection criteria and plots
  P  <- cbind(rep(q[1], length(lw)),
              rep(q[3], length(lw)),
              lw)
  M2  <- test.robustness(X = X.artificial, P = P)
  
  ## Plot 1 - end-member loadings which do not overlap and yielded mRt > 0.80.
  plot(M2$Vqsn[1,], type = "l", ylim = c(0, max(M2$Vqsn, na.rm = TRUE)),
       main = "End-member loadings")
    for (i in 2:nrow(M2$Vqsn)) lines(M2$Vqsn[i,])
  
  # Plot 2 - histogram of mode positions
  hist(M2$modes,
       breaks = 1:ncol(X.artificial), 
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
  
})
