#' Function to create grain-size distributions
#' 
#' The function creates discrete samples of grain-size distributions from
#' continuous particle size measurements, e.g. obtained by SEM-imagery
#' (Francus, 1998), by moving a window of user-defined size with specified
#' step-sizes over the diameter locations.
#' 
#' 
#' @param data Numeric matrix or data frame with two columns: depth
#' (y-coordinates, in metric units) and diameter of the particles (either
#' metric units or phi-scale).
#' @param window.width Numeric scalar, vertical window width. All particles
#' whose y-coordinates fall into this moving window are sampled to build the
#' grain-size distribution for this window.
#' @param step.size Numeric scalar, amount of vertical distance between each
#' sample, i.e. moving window mid point.
#' @param class.limits Numeric vector, specifying the class limits of the
#' grain-size ditribution (either in metric units or phi-scale). If missing,
#' the limits will be set automatically based on the diameter range of the
#' data, with 20 classes, i.e. 21 limits.
#' @return A list with function output \item{sample}{Numeric matrix, sample
#' grain-size distributions.} \item{n}{Numeric vector, number of speciment per
#' sample. }
#' @author Michael Dietze, Elisabeth Dietze
#' @seealso \code{\link{interpolate.classes}}
#' @references Francus P. 1998. An image-analysis technique to measure
#' grain-size variation in thin sections of soft clastic sediments. Sedimentary
#' Geology 121(3-4): 289-298.
#' @keywords EMMA
#' @examples
#' 
#' ## load example data set and create grain-size class vector
#' data(SEM.data, envir = environment())
#' phi <- seq(from = 3, 9, by = 0.5)
#' 
#' ## create samples each 100 micrometers with a 500 micrometers wide window
#' sample <- create.sample(data = SEM.data,
#'                         window.width = 200,
#'                         step.size = 10,
#'                         class.limits = phi)
#' X <- sample$sample
#' n <- sample$n
#' 
#' ## plot grain-size distributions colorised by sample density
#' plot(NA, xlim = range(phi), ylim = range(X), 
#'      main = "Grain-size distributions",
#'      xlab = "Grain-size class", ylab = "Relative amount [%]")
#' for(i in 1:nrow(X)) {lines(phi, X[i,], col = rgb(0, 0, 1, n[i] / max(n)))}
#' 
#' @export create.sample
create.sample <- function(
  data,
  window.width,
  step.size,
  class.limits
){

  ## check correct format of input diameter data
  if(is.data.frame(data) == TRUE) {as.matrix(data)}
  if(ncol(data) != 2) {stop("Data is not a 2-column matrix.")}
  if(is.numeric(data) == FALSE) {stop("Data is non-numeric.")}
  
  ## check/set class limits
  if(missing(class.limits)) {
    class.limits <- seq(from = min(data[,2]), 
                        to = max(data[,2]),
                        length.out = 21)
  }
  if(min(class.limits) > min(data[,2])) {stop(
    "Lowest class limit does not cover minimum diameter")}
  if(max(class.limits) < max(data[,2])) {stop(
    "Highest class limit does not cover maximum diameter")}
  
  ## create auxiliary variables
  y.0 <- min(data[,1]) + window.width / 2
  y.1 <- max(data[,1]) - window.width / 2
  y <- y.0
  samples <- rep(NA, length(class.limits))
  n <- numeric(1)
  
  ## loop through data set
  while(y <= y.1) {
    subsample <- data[(data[,1] >= y - window.width / 2) & 
                        (data[,1] < y + window.width / 2), 2]
    frequency <- hist(subsample, breaks = class.limits, plot = FALSE)$counts
    n[length(n) + 1] <- sum(frequency)
    samples <- rbind(samples, c(0, frequency / max(frequency) * 100))
    y <- y + step.size
  }
  
  ## remove dummy sample and n
  samples <- samples[2:nrow(samples),]
  n <- n[2:length(n)]
  
  ## return sample data set and respective number of specimen per sample
  list(sample = samples,
       n = n)
}