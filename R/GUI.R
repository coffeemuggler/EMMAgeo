#' Start GUI for EMMA
#' 
#' This function starts a browser-based graphic user interface for EMMA. The 
#' GUI has so far been tested on a Linux system, both with the browser of 
#' RStudio and Mozilla Firefox. It permits basic access to import, display 
#' and model a user-provided data set.
#' 
#' To use own data set, this should be a plain ASCII file with samples 
#' organsised as rows and grain-size classes organised as columns. The ASCII 
#' file can be separated by spaces, commas, semi colons or tab stops. The 
#' file may contain a leading column with sample IDs and/or one leading row 
#' with grain-size class breaks. To run EMMA make sure that there are no 
#' classes that contain only zeros throught all samples (i.e., remove them 
#' beforehand, e.g., by \code{X = X[,colSums(X) > 0]}).
#' 
#' @param ... further arguments to pass to \code{\link{runApp}}
#' 
#' @author Michael Dietze
#' @seealso \code{\link{runApp}}
#' @examples 
#' 
#' \dontrun{
#' # Start the GUI
#' GUI()
#' }
#' 
#' @export GUI
GUI <- function(...) {
  app <- shiny::runApp(system.file("shiny/EMMA", 
                                   package = "EMMAgeo"), 
                       launch.browser = TRUE, 
                       ...)
}