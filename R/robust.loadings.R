#' Extract robust end-member loadings
#' 
#' This function takes a list object with potential end-member loadings and 
#' extracts those with modes in specified limits to describe them by mean and
#' standard deviation.
#' 
#' @param em \code{List} of class \code{"EMMAgeo_empot"}, i.e. the outout of 
#' \code{model.EM()} or \code{test.robustness()}, containing potential 
#' end-members, both in unscaled and rescaled version as well as further 
#' parameters.
#' 
#' @param limits \code{Numeric} matrix with two columns, defining the class limits 
#' for the robust end-members to calculate. The first column defines the lower
#' limits, the second column the upper limits. End-members are organised in 
#' rows.
#' 
#' @param classunits \code{Numeric} vector, optional class units 
#' (e.g. micrometers or phi-units) of the same length as the number of 
#' (grain-size) classes per sample.
#' 
#' @param amount \code{Numeric} matrix with two columns, defining the minimum and 
#' maximum amount of the modal class for each end-member.
#' 
#' @param type \code{Character} scalar, type of statistics. One out of 
#' \code{"mean"} and \code{"median"}. Default is \code{"mean"}.
#' 
#' @param qt \code{Numeric} vector of length two, quantiles to describe 
#' end-member loadings Default is \code{c(0.25, 0.75)} (i.e., 
#' the quartile range).
#' 
#' @param plot \code{Logical} scalar, option to enable plot output. Default 
#' is \code{FALSE}.
#' 
#' @param \dots Additional arguments passed to \code{EMMA} and \code{plot}. 
#'  
#' @return \code{List} with statistic descriptions of unscaled and scaled 
#' end-member loadings.
#'  
#' @author Michael Dietze, Elisabeth Dietze
#' @seealso \code{\link{robust.EM}}, \code{\link{robust.scores}}
#' @keywords EMMA
#' @examples
#' 
#' ## load example data set, potential end-members, output of model.EM()
#' data(example_EMpot)
#' 
#' ## define limits for robust end-members
#' limits <- cbind(c(61, 74, 95, 102), 
#'                 c(64, 76, 100, 105))
#' 
#' ## get robust end-member loadings with plot output
#' robust_loadings <- robust.loadings(em = EMpot,
#'                                    limits = limits,
#'                                    plot = TRUE)
#'                     
#' @export robust.loadings
robust.loadings <- function(
  em,
  limits,
  classunits,
  amount,
  type = "mean",
  qt = c(0.25, 0.75),
  plot = FALSE,
  ...
) {
  
  ## check/set amount
  if(missing(amount) == TRUE) {
    amount <- cbind(rep(x = 0, 
                        times = nrow(limits)),
                    rep(x = max(em$Vqsn, na.rm = TRUE), 
                        times = nrow(limits)))
  }
  
  ## check/set class units
  if(missing(classunits) == TRUE) {
    
    classunits <- seq(from = 1, to = ncol(em$X_in)) 
  }
  
  ## read out additional EMMA parameters
  extraArgs <- list(...)
  
  if("col" %in% names(extraArgs)) {
    col_in <- extraArgs$col
  } else{
    col_in <- 1:nrow(limits)
  }
  
  if("log" %in% names(extraArgs)) {
    log_in <- extraArgs$log
  } else{
    log_in <- ""
  }
  
  if("c" %in% names(extraArgs)) {
    c_in <- extraArgs$c
  } else{
    c_in <- 100
  }  
  
  if("rotation" %in% names(extraArgs)) {
    rotation <- extraArgs$rotation
  } else{
    rotation <- "Varimax"
  }  
  
  if("ID" %in% names(extraArgs)) {
    ID <- extraArgs$ID
  } else{
    ID <- seq(from = 1, to = nrow(em$X_in))
  }  
  
  if(ncol(em$X_in) != length(classunits)) {
    
    stop("Units vector is not of same length as data set resolution.")
  }
  
  if(nrow(em$X_in) != length(ID)) {
    
    stop("ID vector is not of same length as samples.")
  }
  
  ## create dummy list structures
  EM.Vqn.list <- vector(mode = "list", 
                        length = nrow(limits))
  
  EM.Vqsn.list <- vector(mode = "list", 
                         length = nrow(limits))
  
  ## select modes that fall into limits for all limit pairs
  for(i in 1:nrow(limits)) {
    
    ## get suitable loadings
    i_ok <- em$modes_classunits >= limits[i,1] & 
      em$modes_classunits <= limits[i,2]
    
    ## identify valid loadings
    Vqn_mode <- em$Vqn[i_ok,]
    Vqsn_mode <- em$Vqsn[i_ok,]
   
    Vqsn_mode_max <- rowMaxs(Vqsn_mode)

    Vqn_mode <- Vqn_mode[Vqsn_mode_max > amount[i, 1] &
                           Vqsn_mode_max < amount[i, 2],]
    
    Vqsn_mode <- Vqsn_mode[Vqsn_mode_max > amount[i, 1] &
                             Vqsn_mode_max < amount[i, 2],]
    
    ## assign valid loadings Vqn
    EM.Vqn.list[[i]] <- Vqn_mode
    
    ## assign valid loadings Vqn
    EM.Vqsn.list[[i]] <- Vqsn_mode
  }
  
  ## identify empty end-member clusters
  i_filled <- do.call(c, lapply(X = EM.Vqn.list, FUN = function(x) {
    
    nrow(x) >= 1
  }))
  
  ## remove empty clusters
  EM.Vqn.list <- EM.Vqn.list[i_filled]
  EM.Vqsn.list <- EM.Vqsn.list[i_filled]
  
  ## announce removed clusters
  if(sum(i_filled) < nrow(limits)) {
    
    limits_empty <- limits[!i_filled,]
    
    warning(paste("No end-members found with modes between", 
                   paste(apply(X = rbind(limits_empty), 
                         MARGIN = 1, 
                         FUN = paste, 
                         collapse = "-"), 
                         collapse = ", ")))
    
    limits <- limits[i_filled,]
  }
  
  ## calculate statistic descriptions for loadings
  EM.Vqn.mean <- lapply(X = EM.Vqn.list, FUN = function(x) {
    apply(X = x, MARGIN = 2, FUN = mean, na.rm = TRUE)
  })
  
  EM.Vqn.sd <- lapply(X = EM.Vqn.list, FUN = function(x) {
    apply(X = x, MARGIN = 2, FUN = sd, na.rm = TRUE)
  })
  
  EM.Vqn.median <- lapply(X = EM.Vqn.list, FUN = function(x) {
    apply(X = x, MARGIN = 2, FUN = median, na.rm = TRUE)
  })
  
  EM.Vqn.qt1 <- lapply(X = EM.Vqn.list, FUN = function(x) {
    apply(X = x, MARGIN = 2, FUN = quantile, probs = qt[1], na.rm = TRUE)
  })
  
  EM.Vqn.qt2 <- lapply(X = EM.Vqn.list, FUN = function(x) {
    apply(X = x, MARGIN = 2, FUN = quantile, probs = qt[2], na.rm = TRUE)
  })
  
  EM.Vqsn.mean <- lapply(X = EM.Vqsn.list, FUN = function(x) {
    apply(X = x, MARGIN = 2, FUN = mean, na.rm = TRUE)
  })
  
  EM.Vqsn.sd <- lapply(X = EM.Vqsn.list, FUN = function(x) {
    apply(X = x, MARGIN = 2, FUN = sd, na.rm = TRUE)
  })
  
  EM.Vqsn.median <- lapply(X = EM.Vqsn.list, FUN = function(x) {
    apply(X = x, MARGIN = 2, FUN = median, na.rm = TRUE)
  })
  
  EM.Vqsn.qt1 <- lapply(X = EM.Vqsn.list, FUN = function(x) {
    apply(X = x, MARGIN = 2, FUN = quantile, probs = qt[1], na.rm = TRUE)
  })
  
  EM.Vqsn.qt2 <- lapply(X = EM.Vqsn.list, FUN = function(x) {
    apply(X = x, MARGIN = 2, FUN = quantile, probs = qt[2], na.rm = TRUE)
  })
  
  ## check if any clusters have remained
  if(nrow(limits) > 0) {
    
    ## prepare plot data for mean or median option
    if(type == "mean") {
      
      Vqn_average <- do.call(rbind, EM.Vqn.mean)
      
      Vqsn_average_plot <- EM.Vqsn.mean
      
      Vqsn_scatter_plot <- vector(mode = "list", 
                                  length = length(Vqsn_average_plot))
      
      for(i in 1:nrow(limits)) {
        Vqsn_scatter_plot[[i]] <- c(EM.Vqsn.mean[[i]] - 1 * EM.Vqsn.sd[[i]],
                                    rev(EM.Vqsn.mean[[i]] + 1 * EM.Vqsn.sd[[i]]))
        
        Vqsn_scatter_plot[[i]][Vqsn_scatter_plot[[i]] < 0] <- 0
      }
      
      classes_plot <- c(classunits, rev(classunits))
      
    } else if(type == "median") {
      
      Vqn_average <- do.call(rbind, EM.Vqn.median)
      
      Vqsn_average_plot <- EM.Vqsn.median
      
      Vqsn_scatter_plot <- vector(mode = "list", 
                                  length = length(Vqsn_average_plot))
      
      for(i in 1:nrow(limits)) {
        Vqsn_scatter_plot[[i]] <- c(EM.Vqsn.qt1[[i]],
                                    rev(EM.Vqsn.qt2[[i]]))
        
        Vqsn_scatter_plot[[i]][Vqsn_scatter_plot[[i]] < 0] <- 0
      }
      
      classes_plot <- c(classunits, rev(classunits))
      
    }
    
    ## optionally create plot output
    if(plot == TRUE) {
      
      ## store old plot parameters
      par_old <- par(no.readonly = TRUE)
      
      plot(x = NA, 
           xlim = range(classunits), 
           ylim = range(c(Vqsn_scatter_plot)),
           main = "Robust loadings",
           xlab = "Class",
           ylab = "Relative amount",
           log = log_in)
      
      for(i in 1:nrow(limits)) {
        polygon(x = classes_plot, 
                y = Vqsn_scatter_plot[[i]], 
                border = NA, 
                col = adjustcolor(col = col_in[i], 
                                  alpha.f = 0.3))
        
        lines(x = classunits, 
              y = Vqsn_average_plot[[i]], 
              col = col_in[i], 
              lwd = 2)
      }
      
      ## restore plot parameters
      par(par_old)
    }
  }


  
  ## create output
  data_out <- list(Vqsn = list(mean = do.call(rbind, EM.Vqsn.mean),
                               sd = do.call(rbind, EM.Vqsn.sd),
                               median = do.call(rbind, EM.Vqsn.median),
                               qt1 = do.call(rbind, EM.Vqsn.qt1),
                               qt2 = do.call(rbind, EM.Vqsn.qt2)),
                   Vqn = list(mean = do.call(rbind, EM.Vqn.mean),
                              sd = do.call(rbind, EM.Vqn.sd),
                              median = do.call(rbind, EM.Vqn.median),
                              qt1 = do.call(rbind, EM.Vqn.qt1),
                              qt2 = do.call(rbind, EM.Vqn.qt2)),
                   X_in = em$X_in,
                   l_in = em$l_in,
                   Vqn_data = EM.Vqn.list,
                   Vqsn_data = EM.Vqsn.list)
  
  ## set class
  class(data_out) <- "EMMAgeo_robload"
  
  ## return output
  return(data_out)
}