#' Generate a vector of weight transformation values from l.min to l.max.
#' 
#' This function generates a sequence of weight transformation values that
#' range from l_min (by default zero) to l_max (by default 95 \% of the 
#' maximum possible value). It uses the function \code{test.l.max()}.
#' 
#' 
#' @param X \code{Numeric} matrix, input data set with m samples (rows) 
#' and n variables (columns).
#' 
#' @param n \code{Numeric} scalar, length of the output vector (by default 10).
#' 
#' @param max \code{Numeric} scalar, fraction of the maximum value 
#' (by default 0.95).
#' 
#' @param min \code{Numeric} scalar, minimum value (by default zero).
#' 
#' @return \code{Numeric} vector of class \code{"EMMAgeo_l"}, weight 
#' transformation values.
#' 
#' @author Michael Dietze, Elisabeth Dietze
#' @seealso \code{\link{test.l.max}}
#' @keywords EMMA
#' @examples
#' 
#' ## load example data set
#' data(example_X)
#' 
#' ## infer l-vector
#' l <- get.l(X = X, 
#'            n = 5, 
#'            max = 0.8, 
#'            min = 0.02)
#' 
#' @export get.l
get.l <- function(
  X, 
  n = 10, 
  max = 0.95,
  min = 0
){
  
  ## estimate maximum possible l-value
  l.max <- test.l.max(X = X, 
                      n = 10)
  
  ## calculate potentially smaller max value
  l.max <- l.max * max
  
  ## generate output sequence
  l <- seq(from = min, 
           to = l.max, 
           length.out = n)
  
  ## set class of l
  class(l) <- "EMMAgeo_l"
  
  ## return output  
  return(l)
}