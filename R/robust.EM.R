#' Extract robust end-members
#' 
#' This function takes a list object with potential end-member loadings and 
#' extracts those with modes in specified limits to describe them by mean and
#' standard deviation and use these descriptions to propagate the uncertainties
#' to end-member scores.
#' 
#' The function is used to extract potential end-member loadings based on their
#' mode positions and, optionally the height of the mode class, and use them to 
#' infer mean and stanard deviation of all 
#' end-members that match the group criteria defined by \code{limits}. These 
#' information are then used to model the uncertainty of the corresponding
#' end-member scores. The function uses input from two preceeding approaches. 
#' In a compact protocol \code{model.em} delivers these data in a predefined 
#' way. In the extended protocol \code{test.robustness} does this.
#' 
#' @param em \code{List} of class \code{"EMMAgeo_empot"}, i.e. the outout of 
#' \code{model.em()} or \code{test.robustness()}, containing potential 
#' end-members, both in unscaled and rescaled version as well as further 
#' parameters.
#' 
#' @param limits \code{Numeric} matrix with two columns, defining the class 
#' limits for the robust end-members to calculate. The first column defines the 
#' lower limits, the second column the upper limits. End-members are organised 
#' in rows.
#' 
#' @param classunits \code{Numeric} vector, optional class units 
#' (e.g. micrometers or phi-units) of the same length as columns of \code{X}.
#' 
#' @param amount \code{Numeric} matrix with two columns, defining the minimum 
#' and maximum amount of the modal class for each end-member.
#' 
#' @param l \code{Numeric} scalar, weight transformation limit for 
#' modelling the average end-member output.
#' 
#' @param mc_n \code{Numeric} scalar, number of Monte Carlo simulations to 
#' estimate end-member scores uncertainty. The default setting is ten times the 
#' product of number of end-members and number of weight transformation limits. 
#' The latter is inherited from \code{model.em()}. To disable modelling of 
#' scores uncertainty, set \code{mc_n = 0}.
#' 
#' @param type \code{Character} scalar, type of oadings statistics. One out of 
#' \code{"mean"} and \code{"median"}. Default is \code{"mean"}.
#' 
#' @param qt \code{Numeric} vector of length two, quantiles to describe 
#' end-member loadings. Default is \code{c(0.25, 0.75)} (i.e., the quartile 
#' range).
#' 
#' @param cores \code{Numeric} scalar, number of CPU cores to be used for 
#' calculations. Only useful in multicore architectures. Default is \code{1} 
#' (single core).
#' 
#' @param plot \code{Logical} scalar, option for plot output. Default is 
#' \code{FALSE}.
#' 
#' @param \dots Additional arguments passed to \code{EMMA} and \code{plot}. 
#'  
#' @return \code{List} with statistic descriptions of end-member loadings 
#' and scores.
#' @author Michael Dietze, Elisabeth Dietze
#' @seealso \code{\link{robust.loadings}}, \code{\link{robust.scores}}
#' @keywords EMMA
#' @examples
#' 
#' \dontrun{
#' 
#' ## load example data set
#' data(example_X)
#' 
#' ## get weight transformation limit vector
#' l <- get.l(X = X)
#' 
#' ## get minimum and maximum number of end-members
#' q <- get.q(X = X, l = l)
#' 
#' ## get all potential model scenarios
#' EM_pot <- model.EM(X = X, q = q, plot = TRUE)
#' 
#' ## define end-member mode class limits
#' limits <- cbind(c(61, 74, 95, 102), 
#'                 c(64, 76, 100, 105))
#' 
#' ## get robust end-members in the default way, with plot output
#' rem <- robust.EM(em = EM_pot,
#'                  limits = limits,
#'                  plot = TRUE)
#'                     
#' ## get robust end-members by only modelling uncertainty in loadings
#' robust_EM <- robust.EM(em = EM_pot, 
#'                        limits = limits, 
#'                        plot = TRUE)
#' 
#' }
#'                     
#' @export robust.EM
robust.EM <- function(
  em,
  limits,
  classunits,
  amount,
  l,
  mc_n,
  type = "mean",
  qt = c(0.25, 0.75),
  cores = 1,
  plot = FALSE,
  ...
) {
  
  ## check/set class units
  if(missing(classunits) == TRUE) {
    
    classunits <- seq(from = 1, to = ncol(em$X_in)) 
  }
  
  ## read out additional EMMA parameters
  extraArgs <- list(...)
  
  ## check/set amount
  if(missing(amount) == TRUE) {
    amount <- cbind(rep(x = 0, 
                        times = nrow(limits)),
                    rep(x = max(em$Vqsn, na.rm = TRUE), 
                        times = nrow(limits)))
  }
  
  ## check/set rotation
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
  
  ##check/set l
  if(missing(l) == TRUE) {
    print("Parameter l missing! Set to 'mRt' by default")
    
    l <- "mRt"
  } 
  
  if(is.numeric(l) == TRUE) {
    
    l_max <- test.l.max(X = em$X_in)
    
    if(l > l_max) {
      warning("l is greater than l_max! Values set to 0.95 * l_max.")
      l <- 0.95 * l_max
    }
  } else {
    
    if(sum(match(x = l, 
                 table = c("mRt", "mRm", "mRn", "mEt", "mEm", "mEn"))) < 1) {
      warning("Keyword for l not defined! Set to 'mRt' automatically.")
      l <- "mRt"
    }
  }

  ## evaluate robust loadings
  robust_loadings <- robust.loadings(em = em, 
                                     limits = limits, 
                                     classunits = classunits,
                                     amount = amount, 
                                     type = type,
                                     qt = qt, 
                                     plot = FALSE)

  ## assign average loadings
  Vqn_average <- robust_loadings$Vqn$mean

  ## optionally evaluate l_opt
  if(is.character(l) == TRUE) {
    
    l <- try(get.l.opt(X = em$X_in, 
                   l = em$l_in, 
                   quality = l,
                   Vqn = Vqn_average, 
                   rotation = rotation, 
                   plot = FALSE), 
             silent = TRUE)
  }
  
  ## check/set Monte Carlo runs
  if(missing(mc_n) == TRUE) {
    mc_n <- 10 * nrow(Vqn_average) * length(em$l_in)
  }
  
  if(class(l) == "try-error") {
    
    stop("No end-members found that match limit definitions!")
  }
  
  ## evaluate robust scores
  robust_scores <- robust.scores(loadings = robust_loadings, 
                                 l = l, 
                                 mc_n = mc_n, 
                                 cores = cores, 
                                 plot = FALSE)

  ## optionally, plot end-member loadings and scores
  if(plot == TRUE) {

    if(type == "mean") {
      
      Vqsn_average_plot <- robust_loadings$Vqsn$mean
      
      Vqsn_scatter_plot <- vector(mode = "list", 
                                  length = nrow(Vqsn_average_plot))
      
      for(i in 1:nrow(Vqsn_average_plot)) {
        Vqsn_scatter_plot[[i]] <- c(
          robust_loadings$Vqsn$mean[i,] - 1 * robust_loadings$Vqsn$sd[i,],
          rev(robust_loadings$Vqsn$mean[i,] + 1 * robust_loadings$Vqsn$sd[i,]))
        
        Vqsn_scatter_plot[[i]][Vqsn_scatter_plot[[i]] < 0] <- 0
      }
      
      classes_plot <- c(classunits, rev(classunits))
      
    } else {
      
      Vqsn_average_plot <- robust_loadings$Vqsn$median
      
      Vqsn_scatter_plot <- vector(mode = "list", 
                                  length = nrow(Vqsn_average_plot))
      
      for(i in 1:nrow(Vqsn_average_plot)) {
        Vqsn_scatter_plot[[i]] <- c(
          robust_loadings$Vqsn$qt1[i,],
          rev(robust_loadings$Vqsn$qt2[i,]))
        
        Vqsn_scatter_plot[[i]][Vqsn_scatter_plot[[i]] < 0] <- 0
      }
      
      classes_plot <- c(classunits, rev(classunits))
    }
    
    if("EM.ID" %in% names(extraArgs)) {
      EM.ID <- extraArgs$EM.ID
    } else {
      EM.ID <- paste("EM", 1:nrow(Vqsn_average_plot), sep = "")
    }
    
    if("main" %in% names(extraArgs)) {
      main <- extraArgs$main
    } else {
      main <- c("Robust end-member loadings",
                "Robust end-member scores")
    }
    
    if("xlab" %in% names(extraArgs)) {
      xlab <- extraArgs$xlab
    } else {
      xlab <- c("Class",
                "Samples")
    }
    
    if("ylab" %in% names(extraArgs)) {
      ylab <- extraArgs$ylab
    } else {
      ylab <- c("Relative amounts",
                "Relative amounts per end-member")
    }
    
    if("ylim" %in% names(extraArgs)) {
      ylim <- extraArgs$ylim
    } else {
      ylim <- rbind(c(0, max(unlist(Vqsn_scatter_plot), na.rm = TRUE)),
                    c(0, 1))
    }
    
    if("log" %in% names(extraArgs)) {
      log <- extraArgs$log
    } else {
      log <- ""
    }
    
    if("col" %in% names(extraArgs)) {
      col <- extraArgs$col
    } else {
      col <- seq(1, nrow(robust_loadings$Vqsn$mean))
    }
    
    if("cex" %in% names(extraArgs)) {
      cex <- extraArgs$cex
    } else {
      cex <- 1
    }
    
    if("lty" %in% names(extraArgs)) {
      lty <- extraArgs$lty
    } else {
      lty <- 1
    }
    
    if("lwd" %in% names(extraArgs)) {
      lwd <- extraArgs$lwd
    } else {
      lwd <- 1
    }
    
    ## setup plot area
    ## save origial parameters
    par.old <- op <- par(no.readonly = TRUE)
    
    ## adjust margins
    par(mar = c(3.9, 4.5, 3.1, 1.1))
    
    ## define layout
    layout(rbind(c(1, 1, 2, 2), 
                 c(1, 1, 2, 2), 
                 c(3, 3, 4, 4), 
                 c(3, 3, 4, 4),
                 c(3, 3, 4, 4), 
                 c(3, 3, 4, 4),
                 c(5, 5, 5, 5))
    )
    
    ## plot col-wise explained variance
    plot(x = classunits, 
         y = robust_scores$Rn, 
         type = "b", 
         main = paste("Class-wise explained variance (mean = ", 
                      signif(x = mean(robust_scores$Rn * 100), 
                             digits = 2), " %)", sep = ""),
         xlab = xlab[1],
         ylab = expression(R^2),
         log = log)
    
    ## plot row-wise explained variance
    plot(x = 1:nrow(em$X_in), 
         y = robust_scores$Rm, 
         type = "b", 
         main = paste("Sample-wise explained variance (mean = ", 
                      signif(x = mean(robust_scores$Rm * 100), 
                             digits = 2), " %)", sep = ""),
         xlab = xlab[2],
         ylab = expression(R^2))
    
    ## plot end-member loadings
    plot(NA,
         main = main[1],
         xlab = xlab[1],
         ylab = ylab[1],
         xlim = range(classes_plot),
         ylim = as.vector(ylim[1,]),
         log = log,
         col = col[1])
    
    if(nrow(Vqsn_average_plot) > 1) {
      for(i in 1:nrow(Vqsn_average_plot)) {
        
        polygon(x = classes_plot, 
                y = Vqsn_scatter_plot[[i]], 
                border = NA, 
                col = adjustcolor(col = col[i], 
                                  alpha.f = 0.3))
        
        lines(x = classunits, 
              y = Vqsn_average_plot[i,], 
              col = col[i])
      }
    }
    
    ## plot end-member scores
    
    Mqs_scatter_plot <- cbind(robust_scores$mean - robust_scores$sd,
                              robust_scores$mean + robust_scores$sd)
    Mqs_scatter_plot <- ifelse(Mqs_scatter_plot < 0, 0, Mqs_scatter_plot)
    Mqs_scatter_plot <- ifelse(Mqs_scatter_plot > 1, 1, Mqs_scatter_plot)
    
    plot(NA, 
         xlim = c(1, nrow(robust_scores$mean)), 
         ylim = c(0, ncol(robust_scores$mean)), 
         axes = FALSE,
         main = main[2],
         xlab = xlab[2],
         ylab = ylab[2])
    box(which = "plot")
    axis(side = 1, at = axTicks(1), labels = ID)
    ticks_scores <- seq(from = 0, 
                        to = ncol(robust_scores$mean), 
                        by = 0.5)
    labels_scores <- c(rep(c(0, NA), 
                           times = ncol(robust_scores$mean)), NA)
    axis(side = 2, 
         at = ticks_scores, 
         labels = labels_scores)
    
    for(i in 1:ncol(robust_scores$mean)) {
      segments(x0 = 1:length(ID),
               y0 = Mqs_scatter_plot[,i] + (i - 1),
               x1 = 1:length(ID),
               y1 = Mqs_scatter_plot[,i + ncol(robust_scores$mean)] + (i - 1),
               col = col[i])
      points(x = 1:length(ID), 
             y = robust_scores$mean[,i] + (i - 1), 
             cex = 0.5, 
             col = col[i])
    }
    
    ## plot legend-like information
    par(mar = c(1, 1, 1, 1))
    
    plot(NA, 
         xlim = c(0, 1), 
         ylim = c(0, 1), 
         axes = FALSE, 
         ann = FALSE, 
         frame.plot = TRUE)
    
    mtext(line = -2, 
          text = "End-member ID (mode position | explained variance)", 
          cex = 0.9 * cex)
    
    ## bugfix included 2017-01-04 | assign classunits for legend output
    classunits_legend <- classunits[robust_scores$modes]

    legend(x = "bottom", 
           legend = paste(EM.ID, " (", signif(x = classunits_legend, 
                                              digits = 2), " | ", 
                          signif(x = robust_scores$Mqs.var, 
                                 digits = 2), " %)", sep = ""), 
           col = col, 
           lty = lty, 
           lwd = lwd, 
           horiz = TRUE, 
           box.lty = 0)
    
    ## reset plot parameters
    par(par.old)
  }
  
  ## return output
  return(list(loadings = robust_loadings$Vqsn,
              scores = list(mean = robust_scores$mean,
                            sd = robust_scores$sd),
              Vqn = robust_loadings$Vqn,
              Xm = robust_scores$Xm,
              modes = robust_scores$modes,
              Mqs.var = robust_scores$Mqs.var,
              Em = robust_scores$Em,
              En = robust_scores$En,
              RMSEm = robust_scores$RMSEm,
              RMSEn = robust_scores$RMSEn,
              Rm = robust_scores$Rm,
              Rn = robust_scores$Rn,
              mRm = robust_scores$mRm,
              mRn = robust_scores$mRn,
              mRt = robust_scores$mRt,
              ol = robust_scores$ol))
}
