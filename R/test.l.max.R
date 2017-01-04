#' Find maximum possible wight transformation value.
#' 
#' This function approximates the highest possible value for l in a nested
#' loop. It uses \code{test.l} and does not need any further parameters. It 
#' starts with l between zero and 0.5 and iteratively approximates the last 
#' possible vlaues for which the weight-transformed matrix of the input data
#' still allows eigenspace extraction.
#' 
#' 
#' @param X \code{Numeric} matrix, input data set with m samples (rows) 
#' and n variables (columns).
#' 
#' @param n \code{Numeric} scalar, number of loop runs and values per loop.
#' 
#' @param \dots Further arguments passed to the function.
#' 
#' @return \code{Numeric} scalar, maximal possible l value.
#' 
#' @author Michael Dietze, Elisabeth Dietze
#' @seealso \code{\link{EMMA}}, \code{\link{test.l}}
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
#' data(example_X)
#' 
#' ## create weight transformation limits vector
#' l <- seq(from = 0, to = 0.6, by = 0.02)
#' 
#' ## test l.max
#' l.max <- test.l.max(X = X)
#' 
#' 
#' 
#' @export test.l.max
test.l.max <- function(
  X,
  n = 10,
  ...
){
  
  ## check for l vs. lw
  if("lw" %in% names(list(...))) {
    stop('Parameter "lw" is depreciated. Use "l" instead.')
  }

  for(i in 1:n) {
    
    ## define initial l vector
    if(i == 1) {
      l.new <- seq(from = 0, 
                    to = 0.5, 
                    length.out = n)
    }
    
    ## test l vector for validity
    l <- test.l(X = X, 
                 l = l.new)
    
    ## update l vector
    l.new <- seq(from = l.new[l$step], 
                  to = l.new[l$step + 1], 
                  length.out = 10)
  }
  
  l.max <- l$l.max
  
  return(l.max)
}