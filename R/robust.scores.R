#' Extract robust end-member scores.
#'
#' This function takes a list object with statistics of end-member loadings and 
#' propagates these uncertainties to end-member scores using Monte Carlo 
#' methods.
#' 
#' @param loadings \code{List} of class \code{"EMMAgeo_robload"}, i.e. the 
#' outout of \code{robust.loadings()}, containing statistic descriptions of 
#' robust end-member loadings.
#' 
#' @param l \code{Numeric} scalar, weight transformation limit to use for 
#' modelling the average end-member output. Can be output of 
#' \code{get.l.opt()}. If omitted, it is set to \code{0}.
#' 
#' @param mc_n \code{Numeric} scalar, number of Monte Carlo simulations to 
#' estimate end-member scores uncertainty. The default setting is ten times the 
#' product of number of end-members and number of weight transformation limits. 
#' The latter is inherited from \code{model.em()}. To disable modelling of 
#' scores uncertainty, set \code{mc_n = 0}.
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
#' @return \code{List} with statistic descriptions of robust end-member 
#' scores.
#' 
#' @author Michael Dietze, Elisabeth Dietze
#' @seealso \code{\link{robust.EM}}, \code{\link{robust.loadings}}
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
#' ## get robust end-member loadings
#' robust_loadings <- robust.loadings(em = EMpot, limits = limits)
#' 
#' ## model end-member scores uncertainties with minimum Monte Carlo runs
#' robust_scores <- robust.scores(loadings = robust_loadings, 
#'                                mc_n = 5, 
#'                                plot = TRUE)
#'                     
#' @export robust.scores
robust.scores <- function(
  loadings,
  l,
  mc_n,
  cores = 1,
  plot = FALSE,
  ...
) {
  
  ## check/set l
  if(missing(l) == TRUE) {
    l <- 0
  }
  
  ## read out additional EMMA parameters
  extraArgs <- list(...)
  
  if("col" %in% names(extraArgs)) {
    col_in <- extraArgs$col
  } else{
    col_in <- 1:nrow(loadings$Vqn$mean)
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
  
  if("classunits" %in% names(extraArgs)) {
    classunits <- extraArgs$classunits
  } else{
    classunits <- seq(from = 1, to = ncol(loadings$X_in))
  }  
  
  if("ID" %in% names(extraArgs)) {
    ID <- extraArgs$ID
  } else{
    ID <- seq(from = 1, to = nrow(loadings$X_in))
  }  
  
  if(ncol(loadings$X_in) != length(classunits)) {
    
    stop("Units vector is not of same length as variables.")
  }
  
  if(nrow(loadings$X_in) != length(ID)) {
    
    stop("ID vector is not of same length as variables.")
  }
  
  ## check/set Monte Carlo runs
  if(missing(mc_n) == TRUE) {
    mc_n <- 10 * nrow(loadings$Vqn$mean) * length(loadings$l_in)
  }
  
  ## set mc parameter l
  mc_l <- rep(x = l, times = mc_n)

  ## set mc parameter Vqn
  mc_Vqn <- vector(mode = "list", length = mc_n)
  
  # for(i in 1:length(mc_Vqn)) {
  #   
  #   Vqn_i <-  matrix(
  #     data = rnorm(n = nrow(loadings$Vqn$mean) * ncol(loadings$Vqn$mean), 
  #                  mean = unlist(loadings$Vqn$mean), 
  #                  sd = unlist(loadings$Vqn$sd)),
  #     nrow = nrow(loadings$Vqn$mean), 
  #     byrow = F)
  #   
  #   mc_Vqn[[i]] <- t(apply(X = Vqn_i, 
  #                          MARGIN = 1, 
  #                          FUN = caTools::runmean, 
  #                          k = 3, 
  #                          endrule = "mean"))
  # }

  for(i in 1:mc_n) {
    
    mc_Vqn[[i]] <- do.call(rbind, lapply(X = loadings$Vqn_data, 
                                         FUN = function(x) {
      x[sample(x = 1:nrow(x), size = 1),]
    }))
  }
  
    
  ## just for investigation
  # plot(NA, xlim = c(1, 116), ylim = c(0, 1), type = "l")
  # for(i in 1:length(mc_Vqn)) lines(mc_Vqn[[i]][1,])
  # for(i in 1:length(mc_Vqn)) lines(mc_Vqn[[i]][2,], col = 2)
  # for(i in 1:length(mc_Vqn)) lines(mc_Vqn[[i]][3,], col = 3)
  # for(i in 1:length(mc_Vqn)) lines(mc_Vqn[[i]][4,], col = 4)
  # 
  
  ## create MC parameter list
  mc_parameters <- vector(mode = "list", length = mc_n)
  
  for(i in 1:length(mc_parameters)) {
    
    mc_parameters[[i]] <- list(Vqn = mc_Vqn[[i]],
                               l = mc_l[i],
                               X = loadings$X_in,
                               c = c_in,
                               rotation = rotation)
  }
  
  ## multicore evaluation
  if(cores > 1) {
    
    ## detect cores
    cores_system <- parallel::detectCores()
    
    ## adjust core number
    if(cores > cores_system) {
      cores <- cores_system
    }
    
    ## initiate cluster
    cl <- parallel::makeCluster(getOption("mc.cores", cores))  
    
    ## do calculations
    mc_scores <- parallel::parLapply (
      cl = cl, 
      X = mc_parameters, 
      fun = function(x) {
        EMMAgeo::EMMA(X = x$X,
                      q = nrow(x$Vqn),
                      l = x$l,
                      c = x$c,
                      Vqn = x$Vqn,
                      rotation = x$rotation)$scores
      })
    
    ## stop cluster
    parallel::stopCluster(cl = cl)
  } else {
    
    ## single-core evaluations
    mc_scores <- lapply(X = mc_parameters, 
                        FUN = function(x) {
                          EMMAgeo::EMMA(X = x$X,
                                        q = nrow(x$Vqn),
                                        l = x$l,
                                        c = x$c,
                                        Vqn = x$Vqn,
                                        rotation = x$rotation)$scores
                        })
  }
  
  ## convert list output to matrix
  mc_scores_unlist <- do.call(cbind, mc_scores)
  
  ## group output by end-members
  mc_scores_by_em <- vector(mode = "list", 
                            length = nrow(loadings$Vqn$mean))
  
  for(i in 1:nrow(loadings$Vqn$mean)) {
    
    index_em <- seq(from = 1, 
                    to = ncol(mc_scores_unlist), 
                    by = nrow(loadings$Vqn$mean))
    index_em <- index_em[-length(index_em)] + (i - 1)
    
    mc_scores_by_em[[i]] <- mc_scores_unlist[,index_em]
  }
  
  
  ## evaluate statistic descriptions for scores
  EM.Mqs.sd <- lapply(X = mc_scores_by_em, FUN = function(x) {
    apply(X = x, MARGIN = 1, FUN = sd, na.rm = TRUE)
  })

  ## evaluate statistic descriptions for scores
  EM.Mqs.mean <- lapply(X = mc_scores_by_em, FUN = function(x) {
    apply(X = x, MARGIN = 1, FUN = mean, na.rm = TRUE)
  })
  
  ## calculate average scores
  EMMA_opt <- EMMAgeo::EMMA(X = loadings$X_in, 
                            q = nrow(loadings$Vqn$mean), 
                            l = l, 
                            c = c_in, 
                            Vqn = loadings$Vqn$mean, 
                            rotation = rotation)
  
  EM.Mqs.average_matrix <- EMMA_opt$scores

  EM.Mqs.average <- vector(mode = "list", length = nrow(loadings$Vqn$mean))
  
  for(i in 1:nrow(loadings$Vqn$mean)) {
    
    EM.Mqs.average[[i]] <- EM.Mqs.average_matrix[,i]
  }
  
  ## prepare plot data for mean or median option
  Mqs_scatter_plot <- vector(mode = "list", 
                             length = nrow(loadings$Vqn$mean))
  
  for(i in 1:nrow(loadings$Vqn$mean)) {
    Mqs_scatter_plot[[i]] <- cbind(EM.Mqs.average[[i]] - 1 * EM.Mqs.sd[[i]],
                                   EM.Mqs.average[[i]] + 1 * EM.Mqs.sd[[i]])
    
    # Mqs_scatter_plot[[i]] <- cbind(EM.Mqs.mean[[i]] - 1 * EM.Mqs.sd[[i]],
    #                                EM.Mqs.mean[[i]] + 1 * EM.Mqs.sd[[i]])
    
    Mqs_scatter_plot[[i]][Mqs_scatter_plot[[i]] < 0] <- 0
    Mqs_scatter_plot[[i]][Mqs_scatter_plot[[i]] > 1] <- 1
  }
  
  classes_plot <- c(classunits, rev(classunits))
  
  ## optionally create plot output
  if(plot == TRUE) {
    
    plot(NA, 
         xlim = c(1, nrow(loadings$X_in)), 
         ylim = c(0, length(EM.Mqs.average)), 
         axes = FALSE,
         main = "Robust scores",
         xlab = "Class",
         ylab = "End-member abundance (0-1)")
    box(which = "plot")
    axis(side = 1)
    ticks_scores <- seq(from = 0, 
                        to = length(EM.Mqs.average), 
                        by = 0.5)
    labels_scores <- c(rep(c(0, NA), 
                           times = length(EM.Mqs.average)), NA)
    axis(side = 2, 
         at = ticks_scores, 
         labels = labels_scores)
    
    for(i in 1:length(Mqs_scatter_plot)) {
      
      segments(x0 = ID,
               y0 = Mqs_scatter_plot[[i]][,1] + (i - 1),
               x1 = ID,
               y1 = Mqs_scatter_plot[[i]][,2] + (i - 1),
               col = col_in[i])
      points(x = ID, 
             y = EM.Mqs.average[[i]] + (i - 1), 
             cex = 0.5, 
             col = col_in[i])
     }
  }
  
  ## return output
  return(list(mean = do.call(cbind, EM.Mqs.average),
              sd = do.call(cbind, EM.Mqs.sd),
              Xm = EMMA_opt$Xm,
              modes = EMMA_opt$modes,
              Mqs.var = EMMA_opt$Mqs.var,
              Em = EMMA_opt$Em,
              En = EMMA_opt$En,
              RMSEm = EMMA_opt$RMSEm,
              RMSEn = EMMA_opt$RMSEn,
              Rm = EMMA_opt$Rm,
              Rn = EMMA_opt$Rn,
              mRm = EMMA_opt$mRm,
              mRn = EMMA_opt$mRn,
              mRt = EMMA_opt$mRt,
              ol = EMMA_opt$ol))
}