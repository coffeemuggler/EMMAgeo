#' Function to test maximum valid lw value.
#' 
#' This function performs the weight transformation of the data matrix after
#' Klovan & Imbrie (1971) and performs EMMA() with different weight limits to
#' check if valied results are yielded. It returns the maximum value for which
#' the transformation remains stable.
#' 
#' 
#' @param X Numeric matrix with m samples (rows) and n variables (columns).
#' @param lw Numeric vector specifying the weight transformation limit, i.e.
#' quantile; default is 0.
#' @return A list with objects \item{step}{Numeric scalar with position of last
#' valid value.} \item{lw.max}{Numeric scalar with last valid value of lw.}
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
#' @export test.lw
test.lw <- function(
  X,
  lw
){
  
  ## check/set default value
  if(missing(lw) == TRUE) {lw = 0}
  
  ## loop through all elements of vector lw
  for(i in 1:length(lw)) {

    ## rescale X constant sum
    X  <- X / apply(X, 1, sum)

    ## calculate weight limit quantiles column-wise
    ls <- sapply(X = 1:ncol(X), FUN = function(j) {
      quantile(x = X[,j], probs = c(lw[i], 1 - lw[i]), type = 5)})

    ## perform weight-transformation
    W <- t((t(X) - ls[1,]) / (ls[2,] - ls[1,]))

    ## optional break when transformation is erroneous
    if (is.na(mean(W))) {
      i = i - 1
      break}
    
    ## optional break when Mqs from EMMA cannot be calculated
    if(is.na(mean(EMMA(X = X, q = 2, lw = lw[i])$Mqs))) {
      i = i - 1
      break
    }
  }
 
  ## assign last valid step number and lw-value
  step  <- i
  lw.max  <- lw[i]
  
  ## return results
  return(list(step = step,
              lw.max = lw.max))
}
