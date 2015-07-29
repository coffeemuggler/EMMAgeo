#' Function to generate GUI server for EMMA.
#' 
#' A graphical user interface (GUI) server backend is started.
#' 
#' @param input A collection of input objects.
#' @param output A collection of input objects.
#' @return A GUI server.
#' @author Michael Dietze, Elisabeth Dietze
#' @keywords EMMA
#' @examples
#' 
#' ## Not run
#' EMMA.GUI()
#' 
#' @export shinyServer

data(X.artificial, envir = environment())
X <- X.artificial

shinyServer(function(
  input, 
  output
  ) {
  
  ## define reactive file loading
  get.data <- reactive({
    
    ## read file
    file.user <- input$file.user
    
    if(is.null(file.user)) {

      return(NULL)
    } else {
     
      return(read.table(file = file.user$datapath, 
                        header = input$header, 
                        sep = input$sep)) 
    }
  })
  
  ## render tabular output
  output$tabular <- renderTable({
    
    ## check if file is loaded and overwrite example data
    if(!is.null(get.data())) {
      X <- get.data()
    }
    
    data.frame(x = X)
  })
  
  ## render GSD plot
  output$graphical <- renderPlot({
    
    ## check if file is loaded and overwrite example data
    if(!is.null(get.data())) {
      X <- get.data()
    }
    
    ## optionally plot curves
    if(input$show == "curves") {
      
      ## create empty plot
      plot(NA, 
           xlim = c(1, ncol(X)), 
           ylim = range(X, na.rm = FALSE))
      
      for(i in 1:nrow(X)) {
        lines(x = 1:ncol(X),
              y = X[i,])
      }
      
      ## optionally plot map
    } else if (input$show == "map") {
      fields::image.plot(t(X))
      
    } else if (input$show == "EMMA") {
      EMMAgeo::EMMA(X = X,
                    q = input$q,
                    lw = input$lw,
                    plot = T)
    }
  })

})
