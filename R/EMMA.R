EMMA <-
structure(function # Function for end-member modelling analysis.
### A multivariate data set (m samples composed of n variables) is decomposed
### by eigenspace analysis and modelled with a given number of end-members 
### (q). Several steps of scaling, transformation, normalisation, eigen 
### space decomposition, factor rotation, data modelling and evaluation are
### performed.
(X,
### Numeric matrix with m samples (rows) and n variables (columns).
q, 
### Numeric scalar with number of end-members to be modelled.
lw,
### Numeric scalar with the weight tranformation limit, i.e. quantiles, 
### cf. Klovan & Imbrie (1971); default is 0.
c,
### Numeric scalar specifying the constant sum scaling parameter, e.g. 1, 
### 100, 1000; default is 100.
Vqn,
### Numeric matrix specifying optional unscaled user-defined end-member
### loadings. If provided, these are used instead of model-derived ones.
EM.ID,
### Character vector with end-member names. If present, these will be set
### as row-names of the output data set and used in the legend text if
### a legend is enabled.
classunits,
### Numeric vector, optional class units (e.g. micrometers or phi-units) of 
### the same length as columns of X.
ID,
### Numeric or character vector, optional sample IDs of the same
### length as rows of X.
rotation = "Varimax",
### Character scalar, rotation type, default is "Varimax" (cf. Dietze et 
### al., 2012). One out of the rotations provided in GPArotation is 
### possible (cf. \code{\link{rotations}}).
plot = FALSE,
### Logical scalar, optional graphical output of the results, default is 
### FALSE. If set to TRUE, end-member loadings and end-member scores are 
### plotted.
legend,
### Character scalar, legend position (cf. \code{\link{legend}}). If
### omitted, no legend will be plotted, default is no legend.
...,
### Additional arguments passed to the plot function. Since the function 
### returns two plots some additional graphical parameters must be
### specified as vector with the first element for the first plot
### and the second element for the second plot. Legend contents and 
### colours apply to both plots. If colours are specified, \code{colour} 
### should be used instead of \code{col}. See example section for 
### further advice.
pm = FALSE
### Logical scalar to enable pm.
) {
   
  ## check/set default values
  if(missing(lw) == TRUE) {lw <- 0}
  if(missing(c) == TRUE) {c <- 100}
  if(missing(EM.ID) == TRUE) {EM.ID <- paste("EM", 1:q, sep = "")}
  
  ## check/set class units vector and test for consistency
  if(missing(classunits) == TRUE) classunits <- 1:ncol(X)
  if(ncol(X) != length(classunits)) stop(
    "Units vector is not of same length as variables.")
  
  ## check/set ID vector and test for consistency
  if(missing(ID) == TRUE) ID <- 1:nrow(X)
  if(nrow(X) != length(ID)) stop(
    "ID vector is not of same length as variables.")
  
  ## End-member modelling
  ## rescale X to constant sum c
  X  <- X / apply(X, 1, sum) * c

  ## calculate weight limit quantiles column-wise
  ls <- sapply(X = 1:ncol(X), FUN = function(i) {
    quantile(x = X[,i], probs = c(lw, 1 - lw), type = 5)})

  ## perform weight-transformation
  W <- t((t(X) - ls[1,]) / (ls[2,] - ls[1,]))

  ## create similarity matrix as outer product
  A <- t(W) %*% W

  ## perform eigenspace decomposition
  EIG <- eigen(A)

  ## assign raw eigenvectors V
  V <- EIG$vectors[,order(seq(ncol(A), 1, -1))]
  Vf <- V[,order(seq(ncol(A), 1, -1))]

  ## rotate eigenvectors and assign factor loadings
  Vr <- do.call(rotation, list(Vf[,1:q]))
  Vq <- Vr$loadings[,seq(q, 1, -1)]

  ## rescale factor loadings and transpose matrix
  Vqr <- t(Vq) / apply(Vq, 2, sum) * c

  ## optionally calculate normalised factor loadings
  if(missing(Vqn) == TRUE) {
    Vqn <- ((Vqr - apply(Vqr, 1, min)) / (
      apply(Vqr, 1, max) - apply(Vqr, 1, min)))
  } else {
  ## or use respective values from input data  
    Vqn <- Vqn}

  ## get eigenvector scores as non-negative least squares estimate
  Mq <- matrix(nrow = nrow(X), ncol = q)
  for (i in 1:nrow(X)) {Mq[i,] = limSolve::nnls(t(Vqn), as.vector(t(W[i,])))$X}

  ## modelled values matrix
  Wm <- Mq %*% Vqn

  ## calculate scaling vector after Miesch (1976)
  ls <- sapply(X = 1:ncol(X), FUN = function(x) {
    quantile(x = X[,x], probs = c(lw, 1 - lw), type = 5)})
  s <- c - sum(ls[1,]) / apply(Vqn * (ls[2,] - ls[1,]), 1, sum)

  ## rescale end-member loadings after Miesch (1976)
  Vqs <- Vqn
  for(i in 1:q) {Vqs[i,] <- s[i] * Vqn[i,] * ((ls[2,] - ls[1,]) + ls[1,])}

  ## normalise end-member loadings
  Vqsn <- Vqs / apply(Vqs, 1, sum) * c

  ## rescale end-member scores
  Mqs <- t(t(Mq) / s) / apply(t(t(Mq) / s), 1, sum)
  
  ## Model evaluation
  ## get number of overlapping end-member loadings
  ol <- numeric(q)
  for (i in 1:q) {
    ol[i] <- if(Vqsn[i, (Vqsn[i,] == max(Vqsn[i,]))] < 
      max(Vqsn[i, (Vqsn[i,] == max(Vqsn[i,]))])) {1} else {0}
  }
  ol <- sum(ol)

  ## evaluate absolute error and explained variances
  Xm <- Mqs %*% Vqs                       # modelled data values
  Mqs.var <- diag(var(Mqs)) / 
    sum(diag(var(Mqs))) * 100             # explained scores variance
  Em <- as.vector(apply(X - Xm, 1, mean)) # absolute row-wise model error
  En <- as.vector(apply(X - Xm, 2, mean)) # absolute column-wise model error
  Rm <- diag(cor(t(X), t(Xm))^2)          # row-wise explained variance
  Rn <- diag(cor(X, Xm)^2)                # column-wise explained variance

  ## Sort output by mode position for consistent output
  ## create auxiliary x-unit and index variable, calculate index
  x.unit <- 1:ncol(X)
  ind <- numeric(q)
  for (i in 1:q) {ind[i] <- x.unit[Vqsn[i,] == max(Vqsn[i,])]}
  ind.sort <- seq(1, nrow(Vqsn))[order(ind)]
  
  Vqsn <- Vqsn[ind.sort,]
  Vqn <- Vqn[ind.sort,]
  Mqs <- Mqs[,ind.sort]
  Mqs.var <- Mqs.var[ind.sort]
  EM.ID <- EM.ID[ind.sort]
  
  ## determine mode class for all end-member loadings
  modes <- numeric(q)
  for(i in 1:q) {modes[i] <- classunits[Vqsn[i,1:ncol(
    Vqsn)] == max(Vqsn[i,1:ncol(Vqsn)])]}
 
  ## optionally, plot end-member loadings and scores
  if(plot == TRUE) {
    
    ## read additional arguments list and check/set default values
    extraArgs <- list(...)
    main <- if("main" %in% names(extraArgs)) {extraArgs$main} else
    {c("End-member loadings",
      "End-member scores")}
    xlab <- if("xlab" %in% names(extraArgs)) {extraArgs$xlab} else
    {c("Classes",
       "Samples")}
    ylab <- if("ylab" %in% names(extraArgs)) {extraArgs$ylab} else
    {c("Amount, relative",
       "Amount, relative")}
    ylim <- if("ylim" %in% names(extraArgs)) {extraArgs$ylim} else
    {rbind(c(0, max(Vqsn, na.rm = TRUE)),
           c(0, 1))}
    log <- if("log" %in% names(extraArgs)) {extraArgs$log} else
    {""}
    colour <- if("colour" %in% names(extraArgs)) {extraArgs$colour} else
    {seq(1, q)}
    if("legend" %in% names(extraArgs)) {extraArgs$legend} else
    {legend.text <- rep(NA, q)
      for(i in 1:q) legend.text[i] <- paste(EM.ID[i], " (", 
        round(modes[i], 2), ")", sep = "")}
    legend.cex <- if("cex" %in% names(extraArgs)) {extraArgs$cex} else
    {1}
    legend.lty <- if("lty" %in% names(extraArgs)) {extraArgs$lty} else
    {1}

    ## setup plot area
    par(mfcol = c(1, 2),
        oma=c(0, 1, 0, 0))
    
    ## plot end-member loadings
    plot(classunits, Vqsn[1,], type = "l", 
         main = main[1],
         xlab = xlab[1],
         ylab = ylab[1],
         ylim = as.vector(ylim[1,]),
         log = log,
         col = colour[1])
    if(nrow(Vqsn) >= 2) {for(i in 2:nrow(Vqsn)) {lines(
      classunits, Vqsn[i,], col = colour[i])}}
    
    if(missing(legend) == FALSE) {legend.position <- legend
                                  legend(x = legend.position,
             legend = legend.text,
             col = colour,
             cex = legend.cex,
             lty = legend.lty)}
    
    ## plot end-member scores
    barplot(t(Mqs), names.arg = ID,
            main = main[2],
            xlab = xlab[2], 
            ylab = ylab[2],
            ylim = as.vector(ylim[2,]),
            col = colour,
            horiz = FALSE)
  }
  ## reset plot area format
  par(mfcol = c(1, 1),
      oma=c(0, 0, 0, 0))
  
  ## optionally add pm
  if(pm == TRUE) {pm <- check.data(matrix(runif(4), ncol = 2),
                                   5, 0.01, 100, invisible = FALSE)}
  
  ## readjust plot margins
  par(oma = c(0, 0, 0, 0))

  ## optionally, assign EM.IDs to Vqsn and Vqn matrices
  if(missing(EM.ID) == FALSE) {
    rownames(Vqn) <- EM.ID
    rownames(Vqsn) <- EM.ID
  }
  
  ##value<< A list with numeric matrix objects.
  list(loadings = Vqsn,    ##<< Normalised rescaled end-member loadings.
       scores   = Mqs,     ##<< Rescaled end-member scores.
       Vqn      = Vqn,     ##<< Normalised end-member loadings.
       Vqsn     = Vqsn,    ##<< Normalised rescaled end-member loadings.
       Mqs      = Mqs,     ##<< Rescaled end-member scores.
       Xm       = Xm,      ##<< Modelled data.
       modes    = modes,   ##<< Mode class of end-member loadings.
       Mqs.var  = Mqs.var, ##<< Explained variance of end-members
       Em       = Em,      ##<< Absolute row-wise model error.
       En       = En,      ##<< Absolute column-wise model error.
       Rm       = Rm,      ##<< Row-wise (sample-wise) explained variance.
       Rn       = Rn,      ##<< Column-wise (variable-wise) explained variance.
       ol       = ol)      ##<< Number of overlapping end-members.
       
  ##end<<
  
  ##details<<
  ## The function values \code{$loadings} and \code{$scores} are redundant.
  ## They are essentially the same as \code{$Vqsn} and \code{$Mqs}. However,
  ## they are included for user convenience.
  
  ##references<<
  ## Dietze E, Hartmann K, Diekmann B, IJmker J, Lehmkuhl F, Opitz S, 
  ## Stauch G, Wuennemann B, Borchers A. 2012. An end-member algorithm for 
  ## deciphering modern detrital processes from lake sediments of Lake Donggi 
  ## Cona, NE Tibetan Plateau, China. Sedimentary Geology 243-244: 169-180.\cr
  ## Klovan JE, Imbrie J. 1971. An Algorithm and FORTRAN-IV Program for 
  ## Large-Scale Q-Mode Factor Analysis and Calculation of Factor Scores. 
  ## Mathematical Geology 3: 61-77.
  ## Miesch AT. 1976. Q-Mode factor analysis of geochemical and petrologic  
  ## data matrices with constant row sums. U.S. Geological Survey 
  ## Professsional Papers 574.
  
  ##seealso<<
  ## \code{\link{test.parameters}}, \code{\link{rotations}},
  ## \code{\link{eigen}}, \code{\link{nnls}}
  
  ##keyword<<
  ## EMMA
}, ex = function(){
  ## load example data and set phi-vector
  data(X.artificial, envir = environment())
  phi <- seq(from = 1, to = 10, length.out = ncol(X.artificial))
  
  ## perform EMMA with 5 end-members
  EM <- EMMA(X = X.artificial, q = 5, lw = 0.05, c = 100, plot = TRUE)
  
  ## perform EMMA with 4 end-members and more graphical settings
  EM <- EMMA(X = X.artificial, q = 4, lw = 0.05, c = 100, 
             plot = TRUE,
             EM.ID = c("EM 1", "EM 2", "EM 3", "EM 4"),
             classunits = phi,
             xlab = c(expression(paste("Class [", phi, "]")), "Sample ID"),
             legend = "topleft", 
             cex = 0.7,
             colour = colors()[c(441, 496, 499, 506)])
})
