#' Model all possible end-member scenarios
#' 
#' This function takes a definition of weight transformation 
#' limits and corresponding minimum and maximum numbers of end-members to 
#' model all end-member scenarios in accordance with these parameters. Based 
#' on the output the user can decide on robust end-members.
#' 
#' The plot output is an overlay of several data. The coloured lines in the 
#' background are end-member loadings (number noted in the plot title), 
#' resulting from all possible model scenarios. If \code{col.q == TRUE} they
#' are coloured according to the number of end-members with which the model 
#' was generated. This colour scheme allows to depict end-members that emerge
#' for model realisations with specific number of end-members. The thick 
#' black line is a kernel density estimate curve, generated from the mode 
#' positions of all end-members. The kernel bandwidth is set to 1 percent of 
#' the number of grain-size classes of the input data set, which gave useful
#' results for most of our test data sets. The cumulaitve dot-line-plot is a
#' further visualisation of end-member mode positions. The function is a 
#' modified wrapper function for the function \code{test.robustness()}.
#' 
#' 
#' @param X \code{Numeric} matrix, input data set with m samples (rows) 
#' and n variables (columns).
#' 
#' @param q \code{Numeric} matrix, definitions of minimum and maximum number
#' of end-members (cf. \code{get.q()}), required.
#' 
#' @param l \code{Numeric} vector, weight transformation limit values, 
#' corresponding to the matrix q, required.
#' 
#' @param classunits \code{Numeric} vector, optional class units 
#' (e.g. micrometers or phi-units) of the same length as columns of \code{X}.
#' 
#' @param plot \code{Logical} scalar, option to plot the results (cf. 
#' details for explanations), default is \code{TRUE}.
#' 
#' @param col.q \code{Logical} scalar, option to colour end-member loadings by  
#' the number of end-members which were used to create the model realisation,
#' default is \code{TRUE}.
#' 
#' @param bw \code{Numeric} scalar, optional manual setting of the kde 
#' bandwidth. By default, bw is calculated as 1 percent of the number of 
#' grain-size classes.
#' 
#' @param \dots Further arguments passed to the function.
#' 
#' @return \code{List} object with all modelled end-members, each described by
#' input parameters, mode position, quality measures and value distributions.
#' 
#' @author Michael Dietze, Elisabeth Dietze
#' 
#' @seealso \code{EMMA}, \code{test.l.max}
#' 
#' @references Dietze E, Hartmann K, Diekmann B, IJmker J, Lehmkuhl F, Opitz S,
#' Stauch G, Wuennemann B, Borchers A. 2012. An end-member algorithm for
#' deciphering modern detrital processes from lake sediments of Lake Donggi
#' Cona, NE Tibetan Plateau, China. Sedimentary Geology 243-244: 169-180.
#' 
#' @keywords EMMA
#' 
#' @examples
#' 
#' ## load example data set
#' data(example_X)
#' 
#' ## define input parameters
#' l <- c(0, 0.05, 0.10)
#' q <- cbind(c(2, 2, 3), c(5, 6, 4))
#' 
#' ## infer l-vector
#' em_pot <- model.EM(X = X, q = q, l = l)
#' 
#' @export model.EM
model.EM <- function(
  X,
  q,
  l,
  classunits,
  plot = TRUE,
  col.q = TRUE,
  bw,
  ...
) {
  
  ## check/set class units
  if(missing(classunits) == TRUE) {
    
    classunits <- seq(from = 1, to = ncol(X)) 
  }
  
  ## check for l vs. lw
  if("lw" %in% names(list(...))) {
    stop('Parameter "lw" is depreciated. Use "l" instead.')
  }
  
  ## check7set log option
  if("log" %in% names(list(...))) {
    
    log_plot <- list(...)$log
  } else {
    
    log_plot <- ""
  }
  
  ## check/set l
  if(class(q)[1] == "EMMAgeo_q") {
    
    l <- as.numeric(rownames(q))
  }
  
  ## check input data
  if(length(l) != nrow(q)) {
    stop("l and q are not of identical length!")
  }
  
  ## check/set bw
  if(missing(bw) == TRUE) {
    bw <- ncol(X) / 100
  }
  
  ## create P-matrix
  P <- cbind(q, l)
  
  ## run test.robustness()
  em <- test.robustness(X = X, P = P, ...)
  
  ## assign class units
  em$modes_classunits <- classunits[em$modes]

  ## assign plot colour
  if(col.q == TRUE) {
    plot_col <- em$q - min(em$q) + 1
  } else {
    plot_col <- rep(x = "grey70",
                    times = nrow(em$loadings))
  }
  
  ## optionally, plot output
  if(plot == TRUE) {
    
    ## create plot object
    em_plot <- em
    
    ## set negative values to zero
    em_plot$loadings[em_plot$loadings < 0] <- 0
    em_plot$Vqsn[em_plot$Vqsn < 0] <- 0
    
    ## get old margins
    mar_old <- par()$mar
    
    ## set new margins
    par(mar = c(5.1, 4.1, 4.1, 4.1))
    
    ## create empty, scaled plot
    plot(NA, 
         xlim = range(classunits), 
         ylim = c(min(em_plot$loadings), 
                  max(em_plot$loadings) * 1.1),
         xlab = "Class",
         ylab = "Contribution", 
         log = log_plot,
         main = paste("Loadings (n = ", 
                      nrow(em_plot$loadings), 
                      ")", 
                      sep = ""))

    ## add all potential loadings
    for(i in 1:nrow(em_plot$loadings)) {
      lines(x = classunits,
            y = em_plot$loadings[i,],
            col = adjustcolor(col = plot_col[i], 
                              alpha.f = 0.3))
    }
    
    ## add cumulate mode position plot
    par(new = TRUE)
    
    plot(x = sort(em_plot$modes),
         y = 1:length(em_plot$modes),
         xlim = c(1,ncol(em_plot$loadings)),
         type = "b",
         lwd = 2,
         col = plot_col,
         ann = FALSE,
         axes = FALSE)
    
    ## add mode position KDE plot
    par(new = TRUE)
    kde <- density(x = em_plot$modes,
                   from = 1, 
                   to = ncol(X),
                   bw = bw)
    plot(kde,
         ylim = c(0, max(kde$y) * 1.2),
         lwd = 2,
         ann = FALSE,
         axes = FALSE)
    axis(side = 4)
    mtext(text = "Density", 
          side = 4, 
          line = 2.5)
  
  if(col.q == TRUE) {
    legend(x = "top",
           legend = paste("q = ", seq(from = min(em_plot$q), 
                                      to = max(em_plot$q))),
           text.col = sort(x = unique(x = plot_col)),
           horiz = TRUE,
           box.lty = 0)
  }
    
    ## create empty plot, scaled to initial plot
    par(new = TRUE)
    
    plot(NA, 
         xlim = range(classunits), 
         ylim = c(min(em_plot$loadings), 
                  max(em_plot$loadings) * 1.1),
         log = log_plot,
         ann = FALSE,
         axes = FALSE)
   
    ## restore old margins
    par(mar = mar_old)
  }
  
  
  ## set class
  class(em)[1] <- "EMMAgeo_empot"

  ## return output
  return(em)
}
