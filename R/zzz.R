##/////////////////////////////////////////////////////////////////////////////
##//zzz.R
##/////////////////////////////////////////////////////////////////////////////
##
##=============================================================================
## author: Michael Dietze
##=============================================================================

##==============================================================================

.onAttach <- function(libname, pkgname){
  
  ## show startup message
  try(packageStartupMessage(paste0(
    "EMMAgeo v. 0.9.6. Please cite this package as: ",
    "Dietze, E., and Dietze, M.: Grain-size distribution unmixing using ",
    "the R package EMMAgeo, E&G Quaternary Sci. J., 69, 1-18, ",
    "https://doi.org/10.5194/egqsj-69-1-2019, 2019. (Don't forget to ",
    "mention the package and R version.)")),
      silent = TRUE)
}
