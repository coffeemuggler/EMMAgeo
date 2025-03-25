#' Check correctness and consistency of input data
#' 
#' The input data matrix (\code{X}), number of end-members (\code{q}), 
#' weight transformation limits (l) and constant sum scaling parameter 
#' (\code{c}) are checked. This includes checking for absence of missing 
#' values, columns containing only zero-values and for numeric data type of 
#' all variables. A further check tests if \code{l} is below the maximum 
#' possible value, preventing numerical instability prior to factor rotation.
#' 
#' 
#' @param X \code{Numeric} matrix, input data set with m samples (rows) 
#' and n variables (columns).
#' 
#' @param q \code{Numeric} scalar, number of end-members to be modelled.
#' 
#' @param l \code{Numeric} scalar or vector, weight transformation
#' limit, i.e.  quantile.
#' 
#' @param c \code{Numeric} scalar, constant sum scaling parameter, e.g.
#' 1, 100, 1000.
#' 
#' @param \dots Further arguments passed to the function.
#' 
#' @return \code{Character} vector, verbose test results.
#' 
#' @author Michael Dietze, Elisabeth Dietze
#' 
#' @seealso \code{EMMA}
#' 
#' @keywords EMMA
#' 
#' @examples
#' 
#' ## load example data set
#' data(example_X)
#' 
#' ## perform data set check
#' check.data(X = X, 
#'            q = 6, 
#'            l = seq(from = 0, 
#'                    to = 0.2, 
#'                    by = 0.01), 
#'            c = 1)
#' 
#' @export check.data
check.data <- function(
  X, 
  q, 
  l,
  c,
  ...
){
  
  ## create result vector
  result <- NA
  
  ## check for l vs. lw
  if("lw" %in% names(list(...))) {
    stop('Parameter "lw" is depreciated. Use "l" instead.')
  }
  
  ## test data type
  result[length(result) + 1] <- ifelse(is.numeric(X) == FALSE, 
    "    Warning: data matrix X contains non-numeric values.",
    "Data matrix passed test... OK")
  
  result[length(result) + 1] <- ifelse(is.numeric(q) == FALSE, 
    "    Warning: end-member vector q has non-numeric values.",
    "End-member vector passed test... OK")
 
  result[length(result) + 1] <- ifelse(is.numeric(l) == FALSE, 
    paste("    Warning: weight transformation limit vector l", 
          "contains non-numeric values."),
    "Weight transformation limit vector passed test... OK")
  
  result[length(result) + 1] <- ifelse(is.numeric(c) == FALSE,
    "    Warning: constant sum scaling parameter c is non-numeric.",
    "Scaling parameter passed test... OK")
  
  ## check for samples with NA-values
  result[length(result) + 1] <- ifelse(sum(complete.cases(X)) < nrow(X),
    paste("    Warning: The following samples comprise NA-values: ",
          seq(1, nrow(X))[!complete.cases(X)],
          ".",
          sep = ""), "NA-test passed... OK")
  if(sum(complete.cases(X)) < nrow(X)) {X <- X[complete.cases(X),]}
  
  ## test if columns contain only 0 values
  X.0 <- apply(X = X, MARGIN = 2, FUN = sum, na.rm = TRUE)
  X.unmet <- seq(1, ncol(X))[X.0 == 0]
  m.unmet <- paste(X.unmet, collapse = ", ")
  if(length(X.unmet) > 0) {
    paste("    Warning: the following columns contain only zero values: ",
      m.unmet, ".",
      sep = "")} else {
        result[length(result) + 1] <- "Test for zero-only values passed... OK"
      }
   
  ## test range of vector l
  l.max <- test.l(X, l)$l.max
  if(length(l.max) == 0) {
    result[length(result) + 1] <- paste("    Warning: weight transformation",
          "limit is out of range.")
  } else {
  result[length(result) + 1] <- ifelse(max(l) > l.max,
    paste("    Note: weight transformation limit(s) are out",
          "of range. Maximum value is", l.max),
    "Maximum weight transformation limit value passed test... OK")
  }
  
  ## test if all samples sum up to constant value
  X.c <- round(apply(X, 1, sum) - rep(c, nrow(X)), 5)
  X.unmet <- seq(1, nrow(X))[X.c != 0]
  m.unmet <- paste(X.unmet, collapse = ", ")
  result[length(result) + 1] <- ifelse(length(X.unmet) >= 1,
      paste("    Note: the following rows do not sum up to the specified c: ",
            m.unmet, ".", sep = ""),
      "All samples sum up to constant sum... OK")
  
  ## return result
  return(result[2:length(result)])
}