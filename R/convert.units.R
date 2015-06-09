#' Function to convert between phi and micrometers.
#' 
#' The function converts values from the phi-scale to the micrometer-scale and
#' vice versa.
#' 
#' 
#' @param phi Numeric vector with grain-size class values in phi to be
#' converted.
#' @param mu Numeric vector with grain-size class values in micrometres to be
#' converted.
#' @return Numeric vector with converted grain-size class values.
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
#' ## convert and show phi to mu
#' convert.units(mu = mu)
#' 
#' @export convert.units
convert.units <-
structure(function(# Function to convert between phi and micrometers.
### The function converts values from the phi-scale to the 
### micrometer-scale and vice versa.
phi,
### Numeric vector with grain-size class values in phi to be converted.
mu
### Numeric vector with grain-size class values in micrometres to be converted.
){
  if(missing(mu) == TRUE){
    ## convert phi to mu
    result <- 1000 * 2^-phi
  } else if(missing(phi) == TRUE){
    ## convert mu to phi
    result <- -log2(mu / 1000)
  } else {stop("No correct variables provided")}

  ## return result
  return(result)
  ### Numeric vector with converted grain-size class values.
  
  ##seealso<<
  ## \code{\link{interpolate.classes}}
  
  ##keyword<<
  ## EMMA
  
}, ex = function(){
  ## generate phi-values
  phi <- -2:5
  
  ## convert and show phi to mu
  mu  <- convert.units(phi = phi)
  mu
  
  ## convert and show phi to mu
  convert.units(mu = mu)
})
