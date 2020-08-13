#' Infer lower and upper mode position limits to define robust end-members.
#' 
#' This function identifies the lower and upper limits within which robust 
#' end-members have clustered mode positions. It uses a kernel density estimate
#' of the mode positions of all input end-member loadings, clips it at a 
#' user-defined minimum density and returns the resulting rising and falling 
#' shoulders of the kde peaks as limits.
#' 
#' Note that the threshold above which a mode cluster is identified is an 
#' arbitrary, user-defined value and probably needs to be adjusted iteratively
#' to get reasonable results. The default value may or may not be adequate! 
#' 
#' 
#' @param loadings \code{R} object, output of function \code{model.EM}.
#' 
#' @param classunits \code{Numeric} vector, optional class units 
#' (e.g. micrometers or phi-units) of the same length as columns of \code{X}.
#' 
#' @param bw \code{Numeric} scalar, bandwidth of the kernel, moved over the 
#' data set. If omitted, the default value of 1 % of the number of classes is
#' used.
#' 
#' @param threshold \code{Numeric} scalar, threshold quantile which is used to 
#' identify mode clusters. Only kde densities above this values are kept and 
#' used to derieve mode cluster limits.
#' 
#' @return \code{Numeric} matrix with lower and upper mode limits.
#' 
#' @author Michael Dietze, Elisabeth Dietze
#' @seealso \code{\link{EMMA}}, \code{\link{model.EM}}
#' @keywords EMMA
#' @examples
#' 
#' ## load example data set
#' data(example_EMpot)
#' 
#' ## infer mode cluster limits
#' limits <- get.limits(loadings = EMpot)
#' 
#' @export get.limits
get.limits <- function(
  loadings,
  classunits,
  bw,
  threshold = 0.7
) {
  
  ## fill mode vector
  if(class(loadings)[1] == "EMMAgeo_empot") {
  
    loadings_mode <- loadings$modes  
  } else {
    
    loadings_mode <- numeric(nrow(loadings))

    for(i in 1:length(loadings_mode)) {
      loadings_mode[i] <- seq(from = 1, 
                              to = ncol(loadings))[
                                loadings[i,] == max(loadings[i,], 
                                                    na.rm = TRUE)]
    }
  }
  
  ## check/set classunits
  if(missing(classunits) == TRUE) {
    
    classunits <- seq(from = 1, 
                      to = ncol(loadings$loadings))
  }
  
  ## check/set bw
  if(missing(bw) == TRUE) {
    
    bw <- (max(loadings_mode) - min(loadings_mode)) / 100
  }
  
  ## create kde of modes
  kde <- density(x = loadings_mode,
                 bw = bw)
  
  ## keep kde parts above threshold value and convert to limits
  kde.ok <- kde$y >= quantile(x = kde$y, 
                              probs = threshold)
  
  kde.limits.1 <- diff(x = kde.ok) == 1
  kde.limits.2 <- diff(x = kde.ok) == -1
  
  ## stop if no meaningful limits can be produced
  if(sum(kde.limits.1) != sum(kde.limits.1)) {
    
    stop("No coherent limits can be found! Set limits manually.")
  }
  
  ## create limits matrix
  limits <- cbind(kde$x[kde.limits.1], 
                  kde$x[kde.limits.2])
  
  ## sort limits row-wise
  limits <- t(apply(X = limits, 
                    MARGIN = 1, 
                    FUN = sort))
  
  ## round limits to integer values
  limits <- round(x = limits, 
                  digits = 0)
  
  ## convert limit classes to class units
  limits <- cbind(classunits[limits[,1]],
                  classunits[limits[,2]])
  
  ## print threshold value
  print(paste("Threshold in KDE density values is", 
              round(x = quantile(x = kde$y, 
                                 probs = threshold),
                    digits = 6)))
  
  ## return result
  return(limits)
}
