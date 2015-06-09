EMMA.GUI <-
structure(function(# Function to run EMMA from a GUI.
  
#  path = "~/R/x86_64-pc-linux-gnu-library/3.2/"
  path = "~/Documents/projects/R-packages/R_EMMAgeo/GIT/"
  
){
  ## run shiny app
  shiny::runApp(paste(path, "/EMMAgeo/R", sep = ""))
}, ex = function(){
  ## not run
  #EMMA.GUI()
  
})
  