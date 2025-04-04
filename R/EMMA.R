#' End-member modelling analysis algorithm.
#' 
#' A multivariate data set (m samples composed of n variables) is decomposed by
#' eigenspace analysis and modelled with a given number of end-members (q).
#' Several steps of scaling, transformation, normalisation, eigenspace
#' decomposition, factor rotation, data modelling and evaluation are performed.
#' 
#' The parameter \code{Vqn} is useful when \code{EMMA} shall be performed with 
#' a set of prior unscaled end-members, e.g. from other data sets that are to 
#' be used as reference or when modelling a data set with mean end-members, as
#' in the output of \code{robust.loadings}.\cr
#' The rotation type \code{Varimax} was used by Dietze et al. (2012). In this 
#' R package, one out of the rotations provided by the package GPArotation 
#' is possible, as well. However, tests showed that the rotation type has no 
#' dramatic consequences for the result.\cr 
#' The function values \code{$loadings} and \code{$scores} are redundant. They 
#' are essentially the same as \code{$Vqsn} and \code{$Mqs}. However, they are 
#' included for user convenience. 
#' 
#' @param X \code{Numeric} matrix, input data set with m samples (rows) 
#' and n variables (columns).
#' 
#' @param q \code{Numeric} scalar, number of end-members to be modelled.
#' 
#' @param l \code{Numeric} scalar or vector, weight transformation
#' limit, i.e.  quantile. Set to zero if omitted.
#' 
#' @param c \code{Numeric} scalar, constant sum scaling parameter, e.g.
#' 1, 100, 1000. Set to 100 if omitted.
#' 
#' @param Vqn \code{Numeric} matrix, optional unscaled user-defined
#' end-member loadings. If provided, these are used instead of model-derived
#' ones. See details.
#' 
#' @param classunits \code{Numeric} vector, optional class units 
#' (e.g. micrometers or phi-units) of the same length as columns of \code{X}.
#' 
#' @param ID \code{Numeric} or character vector, optional sample IDs of the 
#' same length as rows of X.
#' 
#' @param EM.ID \code{Character} vector, end-member names. If present, 
#' these will be set as row-names of the output data set and used in the 
#' legend text.
#' 
#' @param rotation \code{Character} scalar, rotation type, default is 
#' "Varimax". See details.
#' 
#' @param plot \code{Logical} scalar, optional graphical output of the results,
#' default is FALSE. If set to TRUE, end-member loadings and end-member scores
#' are plotted.
#' 
#' @param \dots Additional arguments passed to the plot function. Since the
#' function returns two plots some additional graphical parameters must be
#' specified as vector with the first element for the first plot and the second
#' element for the second plot. 
#' 
#' @return A list with numeric matrix objects. \item{loadings}{Normalised
#' rescaled end-member loadings.} \item{scores}{Rescaled end-member scores.}
#' \item{Vqn}{Normalised end-member loadings.} \item{Vqsn}{Normalised rescaled
#' end-member loadings.} \item{Mqs}{Rescaled end-member scores.}
#' \item{Xm}{Modelled data.} \item{modes}{Mode class of end-member loadings.}
#' \item{Mqs.var}{Explained variance of end-members} \item{Em}{Absolute
#' row-wise model error.} \item{En}{Absolute column-wise model error.} 
#' \item{RMSEm}{row-wise root mean square erroe} \item{RMSEn}{column-wise root 
#' mean square erroe} \item{Rm}{Row-wise (sample-wise) explained variance.} 
#' \item{Rn}{Column-wise (variable-wise) explained variance.} \item{ol}{Number 
#' of overlapping end-members.}
#' 
#' @author Michael Dietze, Elisabeth Dietze
#' 
#' @references Dietze E, Hartmann K, Diekmann B, IJmker J, Lehmkuhl F, Opitz S,
#' Stauch G, Wuennemann B, Borchers A. 2012. An end-member algorithm for
#' deciphering modern detrital processes from lake sediments of Lake Donggi
#' Cona, NE Tibetan Plateau, China. Sedimentary Geology 243-244: 169-180.\cr
#' Klovan JE, Imbrie J. 1971. An Algorithm and FORTRAN-IV Program for
#' Large-Scale Q-Mode Factor Analysis and Calculation of Factor Scores.
#' Mathematical Geology 3: 61-77. Miesch AT. 1976. Q-Mode factor analysis of
#' geochemical and petrologic data matrices with constant row sums. U.S.
#' Geological Survey Professsional Papers 574.
#' 
#' @keywords EMMA
#' 
#' @examples
#' 
#' ## load example data and set phi-vector
#' data(example_X)
#' phi <- seq(from = 1, to = 10, length.out = ncol(X))
#' 
#' ## perform EMMA with 5 end-members
#' EM <- EMMA(X = X, q = 5, l = 0.05, c = 100, plot = TRUE)
#' 
#' ## perform EMMA with 4 end-members and more graphical settings
#' EM <- EMMA(X = X, q = 4, l = 0.05, c = 100, 
#'            plot = TRUE,
#'            EM.ID = c("EM 1", "EM 2", "EM 3", "EM 4"),
#'            classunits = phi,
#'            xlab = c(expression(paste("Class [", phi, "]")), "Sample ID"),
#'            cex = 0.7,
#'            col = rainbow(n = 4))
#' 
#' @export EMMA
EMMA <- function(
 X,
 q, 
 l,
 c,
 Vqn,
 classunits,
 ID,
 EM.ID,
 rotation = "Varimax",
 plot = FALSE,
 ...
) {
  
  ## check for l vs. lw
  if("lw" %in% names(list(...))) {
    
    stop('Parameter "lw" is depreciated. Use "l" instead.')
  }
   
  ## check/set default values
  if(missing(l) == TRUE) {
    
    l <- 0
  }
  
  if(missing(c) == TRUE) {
    
    c <- 100
  }
  
  ## check/set class units vector and test for consistency
  if(missing(classunits) == TRUE) {
    
    classunits <- 1:ncol(X)
  }
  
  if(ncol(X) != length(classunits)) {
    
    stop("Units vector is not of same length as variables.")
  }
  
  ## check/set ID vector and test for consistency
  if(missing(ID) == TRUE) {
    
    ID <- as.character(1:nrow(X))
  }
  
  if(nrow(X) != length(ID)) {
    
    stop("ID vector is not of same length as variables.")
  }
  
  ## check that no zero-only columns exist
  if(any(colSums(X) == 0)) {
    
    stop("X contains columns with only zeros.")
  }
  
  ## End-member modelling
  ## rescale X to constant sum c
  X  <- X / rowSums(X) * c

  ## calculate weight limit quantiles column-wise
  ls <- apply(X <- X, 
              MARGIN = 2, 
              FUN = quantile, 
              probs = c(l, 1 - l), 
              type = 5)
  
  ## perform weight-transformation
  W <- t((t(X) - ls[1,]) / (ls[2,] - ls[1,]))

  ## create similarity matrix as outer product
  A <- t(W) %*% W

  ## perform eigenspace decomposition
  EIG <- try(eigen(A), silent = TRUE)
  
  if(class(EIG)[1] == "try-error") {
    
    stop("Cannot compute eigen space! Consider decreasing l.")
  }

  ## assign raw eigenvectors V
  V <- EIG$vectors[,order(seq(ncol(A), 1, -1))]
  Vf <- V[,order(seq(ncol(A), 1, -1))]

  ## rotate eigenvectors and assign factor loadings
  Vr <- do.call(what = rotation, 
                args = list(Vf[,1:q]))
  Vq <- Vr$loadings[,seq(q, 1, -1)]

  ## rescale factor loadings and transpose matrix
  Vqr <- t(Vq) / colSums(Vq) * c

  ## optionally calculate normalised factor loadings
  if(missing(Vqn) == TRUE) {
    
    ## calculate row-wise min and max values
    Vqr.min <- matrixStats::rowMins(x = Vqr)
    Vqr.max <- matrixStats::rowMaxs(x = Vqr)
    
    ## calculate Vqn
    Vqn <- (Vqr - Vqr.min) / (Vqr.max - Vqr.min)
  } else {
    
    ## or use respective values from input data  
    Vqn <- Vqn
  }

  ## get eigenvector scores as non-negative least squares estimate
  Mq <- t(apply(W, MARGIN = 1, FUN = function(row) {
    
    limSolve::nnls(t(Vqn), as.vector(t(row)))$X
  }))
  
  ## modelled values matrix
  Wm <- Mq %*% Vqn

  ## calculate scaling vector after Miesch (1976)
  ls <- apply(X <- X, 
              MARGIN = 2, 
              FUN = quantile, 
              probs = c(l, 1 - l), 
              type = 5)
  
  s <- c - sum(ls[1,]) / rowSums(Vqn * (ls[2,] - ls[1,]))

  ## rescale end-member loadings after Miesch (1976)
  Vqs <- t(apply(X = cbind(t(t(s)), Vqn), 
                 MARGIN = 1, 
                 FUN = function(row) {
                   
                   row[1] * row[2:length(row)] * (ls[2, ] - ls[1, ]) + ls[1, ]
                   }))

  ## normalise end-member loadings
  Vqsn <- Vqs / rowSums(Vqs) * c

  ## rescale end-member scores
  Mqs <- t(t(Mq) / s) / rowSums(t(t(Mq) / s))
  
  ## Model evaluation
  ## get number of overlapping end-member loadings
  ol <- sum(apply(X = Vqsn, 
                  MARGIN = 1, 
                  FUN = function(row) {
                    if (row[row == max(row)] < max(row[row == max(row)])) { 
                      1 
                    } else { 
                      0 
                      }
                    }))
  
  ## calculate modelled data set
  Xm <- Mqs %*% Vqs

  ## normalise modelled data set
  Xm  <- Xm / rowSums(Xm) * c
  
  ## evaluate absolute error and explained variances
  Mqs.var <- colVars(x = Mqs)
  Mqs.var <- Mqs.var / sum(Mqs.var) * 100 # explained scores variance
  Em <- rowMeans(x = abs(X - Xm))         # absolute row-wise model error
  En <- rowMeans(x = abs(X - Xm))         # absolute column-wise model error
  Rm <- diag(cor(t(X), t(Xm))^2)          # row-wise explained variance
  Rn <- diag(cor(X, Xm)^2)                # column-wise explained variance
  RMSEm <- sqrt(mean(x = Em^2, na.rm = TRUE))
  RMSEn <- sqrt(mean(x = En^2, na.rm = TRUE))
  mRm <- mean(Rm, na.rm = TRUE)
  mRn <- mean(Rn, na.rm = TRUE)
  mRt <- mean(c(mRm, mRn), na.rm = TRUE)

  ## Sort output by mode position for consistent output
  ## create auxiliary x-unit and index variable, calculate index
  x.unit <- 1:ncol(X)
  ind <- numeric(q)
  for (i in 1:q) {
    
    ind[i] <- x.unit[Vqsn[i,] == max(Vqsn[i,])]
  }
  ind.sort <- seq(1, nrow(Vqsn))[order(ind)]
  
  Vqsn <- Vqsn[ind.sort,]
  Vqn <- Vqn[ind.sort,]
  Mqs <- Mqs[,ind.sort]
  Mqs.var <- Mqs.var[ind.sort]
  
  if(missing(EM.ID) == TRUE) {
    EM.ID <- paste("EM", 1:q, sep = "")
  }
  
  ## determine mode class for all end-member loadings
  modes <- numeric(q)
  for(i in 1:q) {modes[i] <- classunits[Vqsn[i,1:ncol(
    Vqsn)] == max(Vqsn[i,1:ncol(Vqsn)])]}
 
  ## optionally, plot end-member loadings and scores
  if(plot == TRUE) {
    
    ## read additional arguments list and check/set default values
    extraArgs <- list(...)
    
    if("main" %in% names(extraArgs)) {
      main <- extraArgs$main
    } else {
      main <- c("End-member loadings",
                "End-member scores")
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
                "Relative amounts")
    }
    
    if("ylim" %in% names(extraArgs)) {
      ylim <- extraArgs$ylim
    } else {
      ylim <- rbind(c(0, max(Vqsn, na.rm = TRUE)),
                    c(-0.04, 1.04))
    }
    
    if("log" %in% names(extraArgs)) {
      log <- extraArgs$log
    } else {
      log <- ""
    }

    if("col" %in% names(extraArgs)) {
      col <- extraArgs$col
    } else {
      col <- seq(1, q)
    }
    
    if("cex" %in% names(extraArgs)) {
      cex <- extraArgs$cex
    } else {
      cex <- 1
    }
    
    if("lty" %in% names(extraArgs)) {
      lty <- extraArgs$lty
    } else {
      lty <- rep(1, q)
    }

    if("lwd" %in% names(extraArgs)) {
      lwd <- extraArgs$lwd
    } else {
      lwd <- rep(1, q)
    }

    ## setup plot area
    ## save origial parameters
    par.old <- op <- par(no.readonly = TRUE)
    
    ## adjust margins
    par(mar = c(3.9, 4.5, 3.1, 1.1),
        cex.axis = cex, 
        cex.lab = cex, 
        cex.main = 1.2 * cex, 
        cex.sub = cex)
    
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
         y = Rn, 
         type = "b", 
         main = paste("Class-wise explained variance (mean = ", 
                      signif(x = mean(Rn * 100), digits = 2), " %)", sep = ""),
         xlab = xlab[1],
         ylab = expression(R^2),
         log = log,
         cex = cex)
    
    ## plot row-wise explained variance
    plot(x = 1:nrow(X), 
         y = Rm, 
         type = "b", 
         main = paste("Sample-wise explained variance (mean = ", 
                      signif(x = mean(Rm * 100), digits = 2), " %)", sep = ""),
         xlab = xlab[2],
         ylab = expression(R^2),
         cex = cex)
    
    ## plot end-member loadings
    plot(classunits, Vqsn[1,], type = "l", 
         main = main[1],
         xlab = xlab[1],
         ylab = ylab[1],
         ylim = as.vector(ylim[1,]),
         log = log,
         col = col[1],
         lty = lty[1],
         lwd = lwd[1],
         cex = cex)
    
    if(nrow(Vqsn) >= 2) {
      for(i in 2:nrow(Vqsn)) {
        lines(x = classunits, 
              y = Vqsn[i,], 
              col = col[i],
              lty = lty[i],
              lwd = lwd[i],
              cex = cex)
      }
    }
    
    ## plot end-member scores
    barplot(t(Mqs), 
            main = main[2],
            xlab = xlab[2], 
            ylab = ylab[2],
            ylim = as.vector(ylim[2,]),
            col = col,
            border = NA,
            space = 0, 
            names.arg = ID,
            horiz = FALSE)
    box(which = "plot")
    
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
          cex = 0.8 * cex, 
          font = 2)
    
    legend(x = "bottom", 
           legend = paste(EM.ID, " (", round(x = modes, digits = 2), " | ", 
                          signif(x = Mqs.var, digits = 2), " %)", sep = ""), 
           col = col, 
           lty = lty, 
           lwd = lwd, 
           horiz = TRUE, 
           box.lty = 0, 
           cex = cex)
    
    ## reset plot parameters
    par(par.old)
  }

  ## optionally, assign EM.IDs to Vqsn and Vqn matrices
  if(missing(EM.ID) == FALSE) {
    rownames(Vqn) <- EM.ID
    rownames(Vqsn) <- EM.ID
  }
  
  ## Output list with numeric matrix objects.
  list(loadings = Vqsn,    ## Normalised rescaled end-member loadings.
       scores   = Mqs,     ## Rescaled end-member scores.
       Vqn      = Vqn,     ## Normalised end-member loadings.
       Vqsn     = Vqsn,    ## Normalised rescaled end-member loadings.
       Mqs      = Mqs,     ## Rescaled end-member scores.
       Xm       = Xm,      ## Modelled data.
       modes    = modes,   ## Mode class of end-member loadings.
       Mqs.var  = Mqs.var, ## Explained variance of end-members
       Em       = Em,      ## Absolute row-wise model error.
       En       = En,      ## Absolute column-wise model error.
       RMSEm    = RMSEm,   ## Sample-wise RMSE.
       RMSEn    = RMSEn,   ## Class-wise RMSE.
       Rm       = Rm,      ## Row-wise (sample-wise) explained variance.
       Rn       = Rn,      ## Column-wise (variable-wise) explained variance.
       mRm      = mRm,      ## Average class-wise explained variance.
       mRn      = mRn,      ## Average sample-wise average explained variance.
       mRt      = mRt,      ## Total average explained variance.
       ol       = ol)      ## Number of overlapping end-members.
}
