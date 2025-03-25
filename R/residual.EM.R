#' Calculate a residual end-member loading.
#' 
#' This function calculates an optional residual end-member loading. It uses
#' the modelled end-member loadings as input and evaluates the root of 1 minus
#' the sum of all squared loadings. The residual end-member can be used to 
#' analyse the remaining variance, e.g.  if not all (robust) EMs are included 
#' (cf. Dietze et al., 2012). Negative values are set to zero.
#' 
#' 
#' @param Vqn \code{Numeric} matrix, m unscaled robust end-member loadings.
#' 
#' @return \code{Numeric} vector, residual end-member loading.
#' 
#' @author Michael Dietze, Elisabeth Dietze
#' 
#' @seealso \code{EMMA}, \code{robust.EM}
#' 
#' @references Dietze E, Hartmann K, Diekmann B, IJmker J, Lehmkuhl F, Opitz S,
#' Stauch G, Wuennemann B, Borchers A. 2012. An end-member algorithm for
#' deciphering modern detrital processes from lake sediments of Lake Donggi
#' Cona, NE Tibetan Plateau, China. Sedimentary Geology 243-244: 169-180.
#' 
#' @keywords EMMA
#' 
#' @examples
#' 
#' ## load example data
#' data(example_X)
#' data(example_EMrob)
#' 
#' ## define mean robust end-member loadings
#' Vqn <- EMMA(X = X, q = 2, plot = TRUE)$loadings
#' 
#' ## perform residual end-member loading calculation
#' Vqn.res <- residual.EM(Vqn)
#' 
#' ## model EMMA with the residual end-member
#' E_res <- EMMA(X = X, 
#'               q = 3, 
#'               Vqn = rbind(Vqn, Vqn.res), 
#'               plot = TRUE)
#' 
#' @export residual.EM
residual.EM <- function(
  Vqn
){
  
  ## transpose Vqn matrix
  Vqn <- t(Vqn)
  
  ## calculate squared residual end-member loading
  res.sq <- 1 - apply(X = Vqn^2, 
                      MARGIN = 1, 
                      FUN = sum)
  
  ## set negative values to zero
  res.sq[res.sq < 0] <- 0
  
  ## calculate residual end-member loading
  Vqn.res <- sqrt(res.sq)
  
  ## return result
  return(Vqn.res)
}