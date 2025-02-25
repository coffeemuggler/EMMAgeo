#' Convert between phi and micrometers.
#' 
#' The function converts values from the phi-scale to the micrometer-scale and
#' vice versa.
#' 
#' 
#' @param phi \code{Numeric} vector, grain-size class values in phi to be
#' converted.
#' 
#' @param mu \code{Numeric} vector, grain-size class values in micrometres 
#' to be converted.
#' 
#' @return \code{Numeric} vector, converted grain-size class values.
#' 
#' @author Michael Dietze, Elisabeth Dietze
#' @seealso \code{\link{interpolate.classes}}
#' @keywords EMMA
#' @examples
#' 
#' ## generate phi-values
#' phi <- -2:5
#' 
#' ## convert and show phi to mu
#' mu  <- convert.units(phi = phi)
#' mu
#' 
#' ## convert and show mu to phi
#' convert.units(mu = mu)
#' 
#' @export convert.units
convert.units <- function(
  phi,
  mu
){
  
  if(missing(mu) == TRUE){
    
    ## convert phi to mu
    result <- 1000 * 2^-phi
  } else if(missing(phi) == TRUE){
    
    ## convert mu to phi
    result <- -log2(mu / 1000)
  } else {
    
    ## return error message
    stop("No correct variables provided")
  }

  ## return result
  return(result)
}
