#' Start GUI for EMMA
#' 
#' This function starts a browser-based graphic user interface for EMMA. The 
#' GUI has so far been tested on a Linux system, both with the browser of 
#' RStudio and Mozilla Firefox. It permits basic access to import, display 
#' and model a user-provided data set.
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