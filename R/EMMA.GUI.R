EMMA.GUI <- structure(function(# Function to run EMMA from a GUI.
  
  path = "~/Documents/projects/R-packages/R_EMMAgeo/EMMAgeo/"
  
){

  ## check if path is correct
  if(missing(path) == TRUE) {
    
    ## search on Linux system
    if(Sys.info()["sysname"] == "Linux") {
      potential.paths <- system('find ~/ -type f -name "X.to.find.EMMAgeo.path.MDT.rda"', 
                                intern = TRUE)
      likely.path <- potential.paths[grepl(x = potential.paths,
                                           pattern = "-pc-linux-")]
      likely.path <- likely.path[length(likely.path)]
    } else if(Sys.info()["sysname"] == "Windows") {
      
      ## search on Windows system
      potential.paths <- system('cd\ & dir /S /P "X.to.find.EMMAgeo.path.MDT.rda"', 
                                intern = TRUE)
      likely.path <- potential.paths[grepl(x = potential.paths,
                                           pattern = "-pc-linux-")]
    }
    
    ## truncate path
    likely.path <- strsplit(x = likely.path, 
                            split = "data/X.to.find.EMMAgeo.path.MDT.rda")[[1]]
    
  } else {
    
    
  }
  
  ## otherwise use provided path
  likely.path <- path
  
  ## check if ui.r and server.r are present under given path
  if(file.exists(paste(likely.path, "R/ui.R", sep = "")) == FALSE |
     file.exists(paste(likely.path, "R/server.R", sep = "")) == FALSE) {
    stop("No Shiny app files present under the suggested path!")
  }
  
  
    ## run shiny app
  shiny::runApp(paste(likely.path, "R", sep = ""))
}, ex = function(){
  ## not run
  #EMMA.GUI()
  
})
  