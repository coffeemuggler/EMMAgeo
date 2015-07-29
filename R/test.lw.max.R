#' Function to find maximum possible lw value.
#' 
#' This function approximates the highest possible value for lw in a nested
#' loop.
#' 
#' 
#' @param X Numeric matrix with m samples (rows) and n variables (columns).
#' @param n Numeric scalar, number of loop runs and values per loop.
#' @return A numeric scalar, maximal possible lw value.
#' @author Michael Dietze, Elisabeth Dietze
#' @seealso \code{\link{EMMA}}, \code{\link{check.data}},
#' \code{\link{test.parameters}}
#' @references Dietze E, Hartmann K, Diekmann B, IJmker J, Lehmkuhl F, Opitz S,
#' Stauch G, Wuennemann B, Borchers A. 2012. An end-member algorithm for
#' deciphering modern detrital processes from lake sediments of Lake Donggi
#' Cona, NE Tibetan Plateau, China. Sedimentary Geology 243-244: 169-180. \cr
#' Klovan JE, Imbrie J. 1971. An Algorithm and FORTRAN-IV Program for
#' Large-Scale Q-Mode Factor Analysis and Calculation of Factor Scores.
#' Mathematical Geology 3: 61-77.
#' @keywords EMMA
#' @examples
#' 
#' ## load example data set
#' data(X.artificial, envir = environment())
#' 
#' ## create weight transformation limits vector
#' lw <- seq(from = 0, to = 0.6, by = 0.02)
#' 
#' ## test the vector
#' test.lw(X = X.artificial, lw = lw)
#' 
#' @export test.lw.max
test.lw.max <- function(
  X,
  n = 10
){

  for(i in 1:n) {
    
    ## define initial lw vector
    if(i == 1) {
      lw.new <- seq(from = 0, 
                    to = 0.5, 
                    length.out = n)
    }
    
    ## test lw vector for validity
    lw <- test.lw(X = X, 
                  lw = lw.new)
    
    ## update lw vector
    lw.new <- seq(from = lw.new[lw$step], 
                  to = lw.new[lw$step + 1], 
                  length.out = 10)
  }
  
  lw.max <- lw$lw.max
  
  return(lw.max)
}