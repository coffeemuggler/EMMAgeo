#' Generate a parameter matrix with q.min and q.max values for robust EMMA.
#' 
#' This function uses the input data matrix \code{X} and a vector of weight
#' transformation limits to generate a matrix of minimum and maximum likely 
#' numbers of end-members to be used to model and extract robust end-members.
#' 
#' The parameter \code{q.min} should be at least 2 because otherwise the 
#' entire dataset would consist of one end-member and there would be no
#' variability at all. The parameter \code{q.max} is set to 10 by default,
#' based on practical issues. In natural systems, there are only rarely 
#' occasions when such a high number of sediment transport regimes may be
#' preserved in and can be resolved from sedimentary deposits. The parameter
#' \code{l} should be a vector between the minimum possible (zero) and maximum
#' possible value (by definition the median, 0.5, but usually a lower value).
#' Whensubmitting only a scalar, the variability can be only due to the range
#' of possible endmembers (between \code{q.min} and \code{q.max}). If the 
#' parameter \code{correct.output} is enabled, this can decrease the number of 
#' valid values for \code{l}, i.e. the number of rows of the output matrix 
#' may no longer be the same as the length of the input vector of \code{l}. In 
#' such a case the vector \code{l} must be replaced by the rownames of the 
#' output matrix (\code{l <- as.numeric(rownames(get.q()))}).
#' 
#' 
#' @param X \code{Numeric} matrix, input data set with m samples (rows) 
#' and n variables (columns).
#' 
#' @param l \code{Numeric} vector, weight transformation limits, default 
#' is zero.
#' 
#' @param q.min \code{Numeric} scalar, minimum number of end-members to use, 
#' default is 2.
#' 
#' @param q.max \code{Numeric} scalar, maximum number of end-members to use, 
#' default is 10.
#' 
#' @param criteria.min \code{Numeric} scalar, minimum value of explained 
#' variance reached to be a valid model realisation, default is 0.5.
#' 
#' @param criteria.max \code{Character} or numeric scalar, either keyword 
#' \code{"local_max"} to use first local maximum or any numeric value of 
#' explained variance, default is \code{"local_max"}.
#' 
#' @param correct.output \code{Logical} scalar, option to correct the output 
#' for twisted values and remove combinations with NA-values. See details.
#' 
#' @param ... Further arguments, passed to the function.
#' 
#' @return \code{Numeric} matrix of class \code{"EMMAgeo_q"}, minimum and 
#' maximum numbers of end-members as well as corresponding weight 
#' transformation values as rownames.
#' 
#' @author Michael Dietze, Elisabeth Dietze
#' 
#' @seealso \code{EMMA}, \code{test.parameters}, \code{test.robustness}
#' 
#' @references Dietze E, Hartmann K, Diekmann B, IJmker J, Lehmkuhl F, Opitz S,
#' Stauch G, Wuennemann B, Borchers A. 2012. An end-member algorithm for
#' deciphering modern detrital processes from lake sediments of Lake Donggi
#' Cona, NE Tibetan Plateau, China. Sedimentary Geology 243-244: 169-180.
#' @keywords EMMA
#' @examples
#' 
#' ## load example data set
#' data("example_X")
#' 
#' ## create parameter matrix
#' get.q(X = X, l = c(0, 0.05, 0.10))
#' 
#' @export get.q
get.q <- function(X, 
                  l = 0, 
                  q.min = 2,
                  q.max = 10,
                  criteria.min = 0.5,
                  criteria.max = "local_max",
                  correct.output = TRUE,
                  ...)
{
  
  ## check for l vs. lw
  if("lw" %in% names(list(...))) {
    stop('Parameter "lw" is depreciated. Use "l" instead.')
  }
  
  ## define q test vector
  q_test <- seq(from = q.min,
                to = q.max)
  
  ## run parameter test routine
  P <- test.parameters(X = X, 
                       q = q_test, 
                       l = l,
                       ...)
  
  ## generate output variables
  q_min_out <- rep(x = NA, 
                   times = length(l))
  q_max_out <- rep(x = NA, 
                   times = length(l))
  
  ## identify q_min and q_max for all values of l
  for(i in 1:length(l)) {

    ## test q_min threshold
    q_min_test <- P$mRt[,i] > criteria.min
    
    if(sum(q_min_test, na.rm = TRUE) > 0) {
      q_min_out[i] <- q_test[q_min_test == TRUE][1]
    }
    
    ## test q_max threshold
    if(criteria.max != "local_max") {
      
      q_max_test <- P$mRt[,i] < criteria.max
      
      if(sum(q_min_test) > 0) {
        q_max_out[i] <- q_test[q_max_test == TRUE][1]
      }
    } else {
      mRt_diff <- diff(x = P$mRt[,i])
      q_max_out[i] <- seq(from = q.min, to = q.max)[mRt_diff < 0][1]
    }
  }
  
  ## combine q-values
  q <- cbind(q_min_out, q_max_out)
  colnames(q) <- c("q_min", "q_max")
  rownames(q) <- l
  
  ## optionally, remove artifacts
  if(correct.output == TRUE) {
    
    ## check for cases where q_min > q_max
    q_twist <- q_min_out > q_max_out
    q_twist[is.na(q_twist)] <- FALSE
    
    ## print notification of removed cases
    if(sum(q_twist) > 0) {
      print(paste(sum(q_twist), 
                  " cases of q_min > q_max removed: ",
                  paste(seq(from = 1, 
                            to = length(l))[q_twist],
                        collapse = ", "), 
                  ".", 
                  sep = ""))
    }
    
    ## set twisted cases to na
    q_min_out[q_twist] <- NA
    
    ## check for na-cases
    na_row <- as.logical(colSums(apply(q, 1, is.na)))
    
    ## print notification of removed cases
    if(sum(na_row) > 0) {
      print(paste(sum(na_row), 
                  " cases with na-presence removed: ",
                  paste(seq(from = 1, 
                            to = length(l))[na_row],
                        collapse = ", "), 
                  ".", 
                  sep = ""))
    }
  }
  
  ## remove na-cases
  q <- cbind(na.omit(object = q))
  
  ## set class
  class(q)[1] <- "EMMAgeo_q"
  
  ## return result
  return(q)
}
