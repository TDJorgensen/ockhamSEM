### Terrence D. Jorgensen
### Last updated: 19 October 2021
### try out lavaanList functionality


## -----------------
## Class and Methods
## -----------------

##' Class for Fitting a Model to Random Data Patterns
##'
##' This class extends \code{\linkS4class{lavaanList}} to specialize in fitting
##' a hypothesized model to random data patterns, thus accounting for model
##' complexity not merely as \emph{df} or number of estimated parameters, but
##' as fitting propensity (Bonifay & Cai, 2017; Preacher, 2003, 2006).
##'
##' @name FitProp-class
##' @aliases FitProp-class  show,FitProp-method  summary,FitProp-method
##   plot,FitProp-method  anova,FitProp-method
##' @docType class
##'
##' @slot lavaanList_slots All slots are available from
##'   \code{\linkS4class{lavaanList}} via \code{...}
##' @slot lavListCall call to \code{\link[lavaan]{lavaanList}} used to fit the
##'   model to the list of random data patterns in \code{@@momentList} or
##'   \code{@@ranDataList}, stored as a \code{list} of arguments
##' @slot convergence \code{list} of \code{logical} vectors indicating whether,
##'   for each random data patterns, (1) the model converged on a solution, (2)
##'   \emph{SE}s could be calculated, (3) the (residual) covariance matrix of
##'   latent variables (\eqn{\Psi}) is non-positive-definite, and (4) the
##'   residual covariance matrix of observed variables (\eqn{\Theta}) is
##'   non-positive-definite.
##' @slot custom \code{logical} indicating whether a custom data-generating
##'   function was provided, rather than using a procedure implemented in the
##'   \code{\pkg{ockhamSEM}} package.
##' @slot momentList \code{list} of summary statistics for random data patterns
##' @slot ranDataList \code{list} of raw data, typically generated from summary
##'   statistics in \code{@@momentList} (e.g., \code{if (!@@custom)}).
##' @slot obsFitMeas \code{\link[lavaan]{fitMeasures}()} output from real data
##' @slot fitMeasDist \code{\link[lavaan]{fitMeasures}()} output from each
##'   random data pattern
##' @slot lavobjects \code{list} of \code{\linkS4class{lavaan}} objects (always
##'   the \code{$target} model, optionally a \code{$baseline.model}).
##'
##' @param object an object of class \code{FitProp}
##' @param ... arguments passed to/from other methods (not used)
##' @param fit.measures see \code{\link[lavaan]{fitMeasures}()}
##' @param lower.tail \code{logical} indicating whether lower values imply
##'   better fit (default: \code{TRUE}, for "badness of fit" indices) or whether
##'   higher values (\code{FALSE}) imply better ("goodness of") fit.
##' @param NML,UIF \code{logical} indicating whether to calculate comparative
##'   fit indices that adjust for a model's fitting propensity
##' @param probs,type see \code{\link[stats]{quantile}}
##' @param omit.reps \code{character} vector of conditions for removing
##'   replications.  \code{"no.conv"} removes replications for which the model
##'   did not converge.  \code{"no.se"} removes replications for which
##'   \emph{SE}s could not be calculated.  \code{"no.hpd"} removes replications
##'   with improper solutions (a nonpositive definite matrix).  \code{integer}s
##'   can also be included in the vector, indicating particular replications
##'   to remove.
##'
##'
##' @return
##'   \item{show}{\code{signature(object = "FitProp")}: Prints a message
##'     about data-generation options and how many models converged.}
##'   \item{summary}{\code{signature(object = "FitProp", ...,
##'     fit.measures = "srmr", lower.tail = TRUE, NML = FALSE, UIF = FALSE,
##'     probs = seq(0, 1, .1), type = 7, omit.reps = c("no.conv","no.se"))}:
##'     For any models in \code{...}, quantiles from the empirical CDFs of
##'     each \code{fit.measures} are provided, along with summary statistics
##'     for model comparison (e.g., UIF).}
##'
##' @author
##' Terrence D. Jorgensen (University of Amsterdam; \email{TJorgensen314@@gmail.com})
##'
##' @references
##' Bonifay, W. E., & Cai, L. (2017). On the complexity of item response theory
##' models. \emph{Multivariate Behavioral Research, 52}(4), 465--484.
##' \doi{10.1080/00273171.2017.1309262}
##'
##' Preacher, K. J. (2003).
##' \emph{The role of model complexity in the evaluation of structural equation models}
##'  (Unpublished doctoral dissertation). The Ohio State University:
##' \url{http://rave.ohiolink.edu/etdc/view?acc_num=osu1054130634}
##'
##' Preacher, K. J. (2006). Quantifying parsimony in structural equation modeling.
##' \emph{Multivariate Behavioral Research, 41}(3), 227--259.
##' \doi{10.1207/s15327906mbr4103_1}
##'
##' @examples
##' ## See ?fitprop help page for examples
##'
setClass("FitProp", contains = "lavaanList",
         slots = c(lavListCall = "list",    # store call to lavaanList()
                   convergence = "list",    # also check SEs and Heywood cases
                   custom      = "logical", # logical: custom data-generator?
                   momentList  = "list",    # store random moments (e.g., correlations)
                   ranDataList = "list",    # store random raw data (e.g., from mvrnorm)
                   obsFitMeas  = "vector",  # store observed fitMeasures()
                   fitMeasDist = "list",    # store propensity distribution for all fitMeasures()
                   lavobjects  = "list"))   # original lavaan object(s)


##' @rdname FitProp-class
##' @aliases show,FitProp-method
##' @export
setMethod("show", "FitProp", function(object) {
  if (!object@custom) {
    dataGen <- paste("Options for R generation\n",
                     "  method =          ", object@call$rmethod, "\n",
                     "  Only positive R?  ", object@call$onlypos, "\n",
                     collapse = "")
  } else dataGen <- "Custom function supplied to generate random patterns.\n"

  nData <- object@meta$ndat
  nConverged <- sum(sapply(object@convergence, "[[", i = "converged"))

  SE <- sapply(object@convergence, "[[", "SE")
  SE[is.na(SE)] <- FALSE

  Heywood.ov <- sapply(object@convergence, "[[", "Heywood.ov")
  Heywood.ov[is.na(Heywood.ov)] <- FALSE

  Heywood.lv <- sapply(object@convergence, "[[", "Heywood.lv")
  Heywood.lv[is.na(Heywood.lv)] <- FALSE

  cat('FitProp object based on ', nData, ' random data patterns. \n', dataGen,
      'See class?FitProp help page for available methods. \n\n',
      'Convergence information:\n', 'The model converged on ',
      nConverged, ' random data patterns. \n\n', sep = "")

  if (!all(SE)) cat('Standard errors could not be computed for data set(s)',
                    paste(which(!SE), collapse = ", "), '\nTry fitting the',
                    'model to individual data patterns to diagnose problems.\n\n')

  if (any(Heywood.ov | Heywood.lv))
    cat('Heywood cases detected for data set(s)',
        paste(which(Heywood.ov | Heywood.lv), collapse = ", "),
        '\nThese are not necessarily a cause for concern because they can',
        'occur due to sampling error, and the model is (by definition)',
        'incorrectly specified for random data patterns. \n\n')

  object
})

##' @importFrom stats quantile na.omit
##' @importFrom matrixStats logSumExp
summarizeFitProp <- function(object, ..., fit.measures = "srmr", lower.tail = TRUE,
                             NML = FALSE, UIF = FALSE,
                             probs = seq(0, 1, .1), type = 7,
                             omit.reps = c("no.conv","no.se")) {
  useReps <- chooseReps(object, omit.reps)

  ## copy/check lower.tail= if necessary
  if (length(lower.tail) == 1L) lower.tail <- rep(lower.tail, length(fit.measures))
  if (length(lower.tail) != length(fit.measures))
    stop('lower.tail should be as long as the specified number of fit.measures=,',
         ' or one logical value to be recycled across fit measures.')
  names(lower.tail) <- fit.measures

  ## save reps-by-fit.measures matrix
  fitDist <- do.call(data.frame, sapply(fit.measures, function(fm) {
    sort(sapply(object@fitMeasDist[useReps], function(x) { x[[fm]] }),
         decreasing = !lower.tail[fm])
  }, simplify = FALSE))

  ## save fit.measures-by-quantiles matrix
  quantiles <- do.call(data.frame, sapply(fitDist, quantile, simplify = FALSE,
                                          probs=probs, type=type, na.rm = TRUE))
  class(quantiles) <- c("lavaan.data.frame","data.frame") # for print() method

  out <- list(obs = object@obsFitMeas[fit.measures],
              propensity = fitDist, quantiles = quantiles)
  class(out$obs) <- c("lavaan.vector","numeric")

  ## Summarize object's fit relative to empirical propensity distribution
  if (!("logl" %in% fit.measures)) NML <- FALSE
  if (NML) {
    logl <- object@obsFitMeas[["logl"]]
    logl.dist <- na.omit(fitDist$logl[fitDist$logl < 0])
    n.ll <- length(logl.dist)
    out$NML <- -2*(logl - (logSumExp(logl.dist) - log(n.ll)))
    class(out$NML) <- c("lavaan.vector","numeric")
  }

  if (UIF) {
    out$UIF <- sapply(fit.measures, function(fm) {
      obs <- object@obsFitMeas[[fm]]
      fpd <- fitDist[[fm]]
      uif <- if (lower.tail[fm]) {
        mean(obs < fpd, na.rm = TRUE)
      } else mean(obs > fpd, na.rm = TRUE)
      uif
    })
    class(out$UIF) <- c("lavaan.vector","numeric")
  }

  class(out) <- c("FitPropSummary","list")
  out
}
##' @rdname FitProp-class
##' @aliases summary,FitProp-method
##' @export
setMethod("summary", "FitProp", summarizeFitProp)

print.FitPropSummary <- function(x, ...) {
  printThese <- c("obs","quantiles")
  if (!is.null(x$NML)) printThese <- c(printThese, "NML")
  if (!is.null(x$UIF)) printThese <- c(printThese, "UIF")
  print(x[printThese])
  invisible(x)
}



##' Compare Fitting Propensities
##'
##' Compare the empirical distributions of models' fitting propensity, or their
##' fit to the data adjusted for fitting propensity
##'
##' @param ... objects of class \code{\linkS4class{FitProp}}
##' @param sameReps \code{logical} indicating whether the whether to use only
##'   replications in which all selected models converged
##' @param omit.reps \code{character} vector of conditions for removing
##'   replications.  \code{"no.conv"} removes replications for which the model
##'   did not converge.  \code{"no.se"} removes replications for which
##'   \emph{SE}s could not be calculated.  \code{"no.hpd"} removes replications
##'   with improper solutions (a nonpositive definite matrix).  \code{integer}s
##'   can also be included in the vector, indicating particular replications
##'   to remove.
##' @param fit.measures see \code{\link[lavaan]{fitMeasures}()}
##' @param lower.tail \code{logical} indicating whether lower values imply
##'   better fit (default: \code{TRUE}, for "badness of fit" indices) or whether
##'   higher values (\code{FALSE}) imply better ("goodness of") fit.
##' @param probs,type see \code{\link[stats]{quantile}}
##' @param NML,UIF \code{logical} indicating whether to calculate comparative
##'   fit indices that adjust for a model's fitting propensity
##' @param conf confidence level for interval estimates of standardized
##'   differences in central tendency. \code{FALSE} indicates to exclude them.
##' @param nd \code{integer}. Number of decimal places rounded in printed output
##'
##' @return
##' For any models in \code{...}, quantiles are provided from the empirical
##' distribution(s) of their \code{fit.measures} to random data, along with
##' pairwise differences between distributions. Summary statistics for model
##' comparison (e.g., UIF) are also available.
##'
##' @author
##' Terrence D. Jorgensen (University of Amsterdam; \email{TJorgensen314@@gmail.com})
##'
##' @references
##' Botha, J. D., Shapiro, A., \& Steiger, J .H. (1988). Uniform indices-of-fit
##' for factor analysis models. \emph{Multivariate Behavioral Research, 23}(4),
##'  443--450. \doi{10.1207/s15327906mbr2304_2}
##'
##' Rissanen, J. (2001). Strong optimality of the normalized ML models as
##' universal codes and information in data.
##' \emph{IEEE Transactions on Information Theory, 47}(5), 1712--1717.
##' \doi{10.1109/18.930912}
##'
##' @examples
##' ## See ?fitprop help page for examples
##'
##' @export
compareFitProp <- function(..., sameReps = TRUE, omit.reps = c("no.conv","no.se"),
                           fit.measures = "srmr", lower.tail = TRUE,
                           probs = seq(0, 1, .1), type = 7,
                           NML = FALSE, UIF = FALSE, conf = FALSE, nd = 3L) {
  if (conf) stopifnot(conf < 1 & conf > 0)

  ## separate models from lists of models
  dots <- list(...)
  idx.list <- sapply(dots, is.list)
  modLists <- dots[ idx.list]
  mods     <- dots[!idx.list]
  ## capture names of any arguments passed via dots
  allnames  <- sapply(substitute(list(...))[-1], deparse)
  listnames <- allnames[ idx.list]
  modnames  <- allnames[!idx.list]

  ## make sure models are named
  if (length(mods) && is.null(names(mods))) {
    names(mods) <- modnames
  } else for (nn in seq_along(mods)) {
    if (names(mods)[nn] == "") names(mods)[nn] <- modnames[nn]
  }
  ## make sure lists are named
  if (length(modLists) && is.null(names(modLists))) {
    names(modLists) <- listnames
  } else for (nn in seq_along(modLists)) {
    if (names(modLists)[nn] == "") names(modLists)[nn] <- listnames[nn]
  }
  ## within each list, make sure models are named
  for (i in seq_along(modLists)) {
    if (length(modLists[[i]]) && is.null(names(modLists[[i]]))) {
      names(modLists[[i]]) <- seq_along(modLists[[i]])
    } else for (nn in seq_along(modLists[[i]])) {
      if (names(modLists[[i]])[nn] == "") names(modLists[[i]])[nn] <- nn
    }
  }

  ## collapse into a single list of models
  if (length(modLists)) mods <- c(mods, unlist(modLists))

  ## check class
  not.FitProp <- !sapply(mods, inherits, what = "FitProp")
  if (any(not.FitProp)) stop("The following are not FitProp objects:\n",
                             paste0(names(which(not.FitProp)), collapse = ", "))
  ## check convergence
  nonConv <- !sapply(mods, function(fit) {
    any(sapply(fit@convergence, "[", i = "converged"))
  })
  if (all(nonConv)) {
    stop('No models converged on any random data patterns')
  } else if (any(nonConv)) {
    message('The following models did not converge on any random data patterns',
            ', so they are ignored:\n',
            paste(names(nonConv)[nonConv], collapse = ",\t"))
    mods <- mods[which(!nonConv)]
  }

  ## Save indices to subset
  convMat <- do.call(cbind, sapply(mods, chooseReps,
                                   omit.reps = omit.reps, simplify = FALSE))
  if (sameReps) {
         conv <- apply(convMat, 1, all)
  } else conv <- apply(convMat, 1, any)

  ## copy/check lower.tail= if necessary
  if (length(lower.tail) == 1L) lower.tail <- rep(lower.tail, length(fit.measures))
  if (length(lower.tail) != length(fit.measures))
    stop('lower.tail should be as long as the specified number of fit.measures=,',
         ' or one logical value to be recycled across fit measures.')
  names(lower.tail) <- fit.measures


  ## UNIVARIATE INFO COMES FROM summary() METHOD
  sumList <- lapply(mods, summarizeFitProp, fit.measures = fit.measures,
                    lower.tail = lower.tail, NML = NML, UIF = UIF,
                    probs = probs, type = type, omit.reps = omit.reps)
  quantiles <- sapply(fit.measures, function(fm) {
    qq <- do.call(data.frame, sapply(sumList, function(m) m$quantiles[,fm],
                                simplify = FALSE))
    class(qq) <- c("lavaan.data.frame","data.frame")
    rownames(qq) <- rownames(sumList[[1]]$quantiles)
    qq
  }, simplify = FALSE)


  ## get degrees of freedom for each model to label output matrix
  DFs <- sapply(mods, function(x) x@obsFitMeas[["df"]] )
  DFs <- paste0(names(DFs), " (df = ", DFs, ")")
  template <- matrix(NA, length(mods), length(mods),
                     dimnames = list(DFs, names(mods)))
  diag(template) <- TRUE # can't compare a model with itself
  class(template) <- c("lavaan.matrix.symmetric","matrix") # for print() method

  ## Loop through sorted models in sequence of most to least restricted model
  out.pc <- list()
  for (fm in fit.measures) {
    ## copy template to store each comparison (using this fit.measure)
    KS <- Cd <- Hg <- Cliff <- template
    for (R in 2:nrow(template)) {
      for (C in (R - 1):1) {
        pairwise <- x.vs.y(mods[[R]], mods[[C]], sameReps=sameReps, useReps=conv,
                           fit.measures = fm, lower.tail = lower.tail[fm])
        if (conf) {
          ## K-S test
          pv <- round(pairwise[[fm]]$KS$p.value, nd)
          pv <- ifelse(pv > 0, yes = paste0("(p = ", pv, ")"),
                       no = paste0("(p < .", rep(0, nd - 1L), 1, ")"))
          KS[R, C] <- KS[C, R] <- paste(pairwise[[fm]]$KS$statistic, pv)

          ## Effect Sizes
          Cd[R, C] <- paste0(round(pairwise[[fm]]$Cohen.d$estimate, nd), " [",
                             round(pairwise[[fm]]$Cohen.d$conf.int[1], nd), ", ",
                             round(pairwise[[fm]]$Cohen.d$conf.int[2], nd), "]")
          Hg[R, C] <- paste0(round(pairwise[[fm]]$Hedges.g$estimate, nd), " [",
                             round(pairwise[[fm]]$Hedges.g$conf.int[1], nd), ", ",
                             round(pairwise[[fm]]$Hedges.g$conf.int[2], nd), "]")
          Cliff[R, C] <- paste0(round(pairwise[[fm]]$Cliff.d$estimate, nd), " [",
                                round(pairwise[[fm]]$Cliff.d$conf.int[1], nd), ", ",
                                round(pairwise[[fm]]$Cliff.d$conf.int[2], nd), "]")
          ## set above-diagonal value to opposite sign
          Cd[C, R] <- paste0(round(-1*pairwise[[fm]]$Cohen.d$estimate, nd), " [",
                             round(-1*pairwise[[fm]]$Cohen.d$conf.int[2], nd), ", ",
                             round(-1*pairwise[[fm]]$Cohen.d$conf.int[1], nd), "]")
          Hg[C, R] <- paste0(round(-1*pairwise[[fm]]$Hedges.g$estimate, nd), " [",
                             round(-1*pairwise[[fm]]$Hedges.g$conf.int[2], nd), ", ",
                             round(-1*pairwise[[fm]]$Hedges.g$conf.int[1], nd), "]")
          Cliff[C, R] <- paste0(round(-1*pairwise[[fm]]$Cliff.d$estimate, nd), " [",
                                round(-1*pairwise[[fm]]$Cliff.d$conf.int[2], nd), ", ",
                                round(-1*pairwise[[fm]]$Cliff.d$conf.int[1], nd), "]")

        } else {
          KS[R, C]    <- KS[C, R] <- pairwise[[fm]]$KS$statistic
          Cd[R, C]    <- pairwise[[fm]]$Cohen.d$estimate
          Hg[R, C]    <- pairwise[[fm]]$Hedges.g$estimate
          Cliff[R, C] <- pairwise[[fm]]$Cliff.d$estimate
          ## set above-diagonal value to opposite sign
          Cd[C, R]    <- -1*Cd[R, C]
          Hg[C, R]    <- -1*Hg[R, C]
          Cliff[C, R] <- -1*Cliff[R, C]
        }
        ## end nested for loops
      }
    }
    ## save all for this fit.measure, before erasing in next iteration
    out.pc[[fm]] <- list(KS = KS, Cohen.d = Cd, Hedges.g = Hg, Cliff.d = Cliff)
  }

  out <- list(quantiles = quantiles, pairwise = out.pc)

  ## NML (using "logLik") & UIF (using any index)
  if (!("logl" %in% fit.measures)) NML <- FALSE
  if (NML) {
    out$NML <- sapply(sumList, function(m) m$NML)
    class(out$NML) <- c("lavaan.vector","numeric")
  }
  if (UIF) {
    out$UIF <- as.data.frame(do.call(rbind, sapply(sumList, function(m) m$UIF,
                                                   simplify = FALSE)))
    class(out$UIF) <- c("lavaan.data.frame","data.frame")
  }


  class(out) <- "compareFitProp"
  out
}

print.compareFitProp <- function(x, ..., nd = 3) {
  cat("
    These are quantiles from each model's empirical distribution
    of fitting propensity.  For the same quantile, models with
    better fit have better fitting propensity.\n\n")
  ## print table of quantiles
  print(x$quantiles, ...)

  ## loop over fit.measures for matrices of pairwise comparisons
  cat("
    The KS statistic is always positive, but positive (standarized)
    differences in central tendency (Cohen's d, Hedges' g, Cliff's delta)
    indicate the model in Column C has greater propensity to fit random
    data patterns than the model in Row R.

    Missing values (NA) indicate models did not jointly converge on any
    random data patterns.\n\n")
  for (fm in seq_along(x$pairwise)) {
    cat("========== ", toupper(names(x$pairwise)[fm]), ":\n\n")
    matList <- x$pairwise[[fm]]
    matList <- sapply(matList, function(mat) {
      if (is.numeric(mat)) mat <- round(mat, nd)
      class(mat) <- "matrix" # "lavaan.matrix" won't work
      mat[upper.tri(mat, diag = TRUE)] <- ""
      mat
    }, simplify = FALSE)
    print(matList, quote = FALSE, ...)
    cat("\n")
  }

  if (!all(is.null(x$NML), is.null(x$UIF))) {
    cat("
    These indices are designed to compare the fit of models, adjusting
    for model complexity (i.e., fitting propensity).\n\n")
  }

  ## print table of NML/UIF(s)
  if (!is.null(x$NML)) {
    cat("========== (-2*log) NML: lower is preferred \n")
    print(x$NML, ..., nd = nd)
    cat("\n\n")
  }
  if (!is.null(x$UIF)) {
    cat("========== Uniform Indices of Fit: higher is preferred \n")
    print(x$UIF, ..., nd = nd)
    cat("\n")
  }

  invisible(x)
}



## --------------------
## Constructor Function
## --------------------

##' Fit Model to Random Data Patterns
##'
##' Estimate a model's fitting propensity to evaluate model complexity
##' beyond merely degrees of freedom
##'
##' @details
##' Inspired by work by Preacher (2003, 2006) and Bonifay & Cai (2017),
##' this function performs three steps for analyses to assess the fit propensity of competing
##' structural equation models:
##' 1. Randomly generate correlation (or covariance matrices);
##' 2. Fit models to each correlation matrix; and
##' 3. Save a indices that could be used for evaluating model fit in subsequent summaries.
##' Conceptually, models that exhibit better fit to such randomly generated data
##' may have better fit propensity, and are therefore potentially less parsimonious.
##'
##' Models must be fitted with the \code{\pkg{lavaan}} package (e.g., created
##' from \code{\link[lavaan]{cfa}}, \code{\link[lavaan]{sem}}, or
##' \code{\link[lavaan]{lavaan}} functions).
##' Currently, only models using ML estimation are supported.
##' The underlying \code{\link[lavaan]{lavOptions}} from the fitted
##' \code{\linkS4class{lavaan}} models will be re-used when fitting the
##' same model to each random data pattern, which will be saved in the
##' returned object.
##' A \code{summary} methods is available to summarize results for a single
##' model, and \code{\link{compareFitProp}} is available for model comparison.
##'
##' Generation of random correlation matrices is provided using several approaches. The \code{"mcmc"}
##' algorithm implements a Markov Chain Monte Carlo approach and was ported from Fortran code
##' in Preacher (2003). For details on the algorithm's actual implementation, see Preacher (2003),
##' Falk and Muthukrishna (in press), or the source code for the \code{mcmc()} function. If this algorithm
##' is chosen, \code{mcmc.args} accepts a list that can modify some default settings. In particular,
##' \code{iter} sets the total number of iterations to run (default = 5000000). If parallel processing
##' is enabled, this number will be divided among the number of chains. \code{miniter} sets a
##' minimum number of iterations per chain to avoid many processors leading to too few iterations per
##' chain (default = 10000). \code{jmpsize} overrides the step size for each update to the candidate
##' correlation matrix. Smaller step sizes typically lead to more acceptance and may be necessary for
##' larger correlation matrices (default jump size depends on the number of variables). Though, in
##' general the MCMC algorithm becomes more difficult to work well with many variables.
##'
##' The \code{"onion"} method is one approach that relies on work of Joe (2006) and
##' Lewandowski, Kurowick, and Joe (2009); matrices are generated recursively, one variable at a
##' time. The onion method is computationally more efficient than the MCMC algorithm. Under the
##' hood, the \code{\link[clusterGeneration]{genPositiveDefMat}} function in the
##' \code{\pkg{clusterGeneration}} package is used, with default
##' arguments of \code{covMethod="onion"}, \code{eta=1}, and \code{rangeVar=c(1,1)}. These arguments ensure that the
##' Onion method is used, generation is uniform over the space of positive definite matrices
##' (but see note on positive manifold below), and with unit variances.
##'
##' An additional option \code{"clustergen"} is provided for direct interface with the \code{\link[clusterGeneration]{genPositiveDefMat}}
##' function in the \code{\pkg{clusterGeneration}} package. A named list can be passed to \code{clustergen.args} to
##' override any defaults used by \code{\link[clusterGeneration]{genPositiveDefMat}}, and the user is referred to documentation
##' for that function. This allows, for example, generation using C-Vines, covariance matrices
##' (i.e., variables that do not all have unit variances), and several other covariance/correlation
##' matrix generation techniques.
##'
##' The \code{onlypos=} argument controls whether correlation matrices are
##' restricted to have only positive correlations.
##' The original MCMC algorithm by Preacher (2003, 2006) generated correlation
##' matrices with positive manifold only (i.e., only positive correlations).
##' The algorithm is easily changed to allow also negative correlations.
##' The \code{"onion"} method and any functions from \code{\pkg{clusterGeneration}}
##' by default generate matrices with both positive and negative correlations.
##' To obtain matrices with positive manifold only, an ad-hoc correction is
##' implemented for these latter approaches where the matrix \eqn{R} is transformed:
##'
##' \deqn{R = (R+1)/2}.
##'
##' To our knowledge, there is no guarantee that this will result in uniform
##' sampling from the space of all correlation matrices with positive manifold,
##' yet fit propensity results for some examples are very similar to those
##' of the MCMC algorithm.
##'
##'
##' @importFrom lavaan lavInspect lavNames parTable lavaanList
##' @importFrom methods as
##'
##' @param object an object of class \code{\linkS4class{lavaan}}
##' @param use.FitProp Optional object of class \code{\linkS4class{FitProp}}.
##'   Assuming \code{object} is fitted to the same variables, providing this
##'   can save computing time by fitting the model in \code{object} to the same
##'   generated data patterns in \code{use.FitProp}, in which case no
##'   further arguments should be specified (they will be recycled).
##' @param reps \code{integer} indicating how many random data patterns to generate
##' @param seed \code{integer} passed to \code{\link{set.seed}}
##' @param ... Additional arguments passed to \code{\link[lavaan]{lavaanList}}
##' @param rmethod Optional \code{character} indicating which among a set of
##'   methods to use, which are implemented in the \pkg{ockhamSEM} package
##'   (see \strong{Details}). Required when not specifying a custom
##'   \code{dataFunction=} via \dots
##' @param onlypos \code{logical} indicating whether to restrict random data
##'   patterns to be positive manifold (see \strong{Details}).  Ignored when a
##'   custom \code{dataFunction=} is provided via \dots
##' @param mcmc.args Optional named \code{list} of arguments that controls
##'   options when \code{rmethod = "mcmc"} (see \strong{Details}).
##' @param clustergen.args Optional named \code{list} of arguments that controls
##'   options when \code{rmethod = "onion"} or \code{"clustergen"}
##'   (see \strong{Details}).
##' @param baseline.model Optional object of class \code{\linkS4class{lavaan}}
##'   to be passed to \code{\link[lavaan]{fitMeasures}}
##'
##' @return
##'   An object of class \code{\linkS4class{FitProp}}, which inherits from
##'   \code{\linkS4class{lavaan}}.
##'
##' @seealso
##' \code{\linkS4class{FitProp}}, \code{\link{compareFitProp}}
##'
##' @author
##' Terrence D. Jorgensen (University of Amsterdam; \email{TJorgensen314@@gmail.com})
##'
##' @references
##' Bonifay, W. E., & Cai, L. (2017). On the complexity of item response theory
##' models. \emph{Multivariate Behavioral Research, 52}(4), 465--484.
##' \doi{10.1080/00273171.2017.1309262}
##'
##' Falk, C. F., & Muthukrishna, M. (in press). Parsimony in model selection:
##' Tools for assessing fit propensity. \emph{Psychological Methods}.
##' \doi{10.1037/met0000422}
##'
##' Lewandowski, D., Kurowicka, D., & Joe, H. (2009). Generating random
##' correlation matrices based on vines and extended onion method.
##' \emph{Journal of Multivariate Analysis, 100}(9), 1989--2001.
##' \doi{10.1016/j.jmva.2009.04.008}
##'
##' Joe, H. (2006). Generating random correlation matrices based on partial
##' correlations. \emph{Journal of Multivariate Analysis, 97}(10), 2177--2189.
##' \doi{10.1016/j.jmva.2005.05.010}
##'
##' Preacher, K. J. (2003).
##' \emph{The role of model complexity in the evaluation of structural equation models}
##'  (Unpublished doctoral dissertation). The Ohio State University:
##' \url{http://rave.ohiolink.edu/etdc/view?acc_num=osu1054130634}
##'
##' Preacher, K. J. (2006). Quantifying parsimony in structural equation modeling.
##' \emph{Multivariate Behavioral Research, 41}(3), 227--259.
##' \doi{10.1207/s15327906mbr4103_1}
##'
##' @examples
##'
##' # Set up a covariance matrix to fit models to
##' p <- 3 # number of variables
##' temp_mat <- diag(p) # identity matrix
##' colnames(temp_mat) <- rownames(temp_mat) <- paste0("V", seq(1, p))
##'
##' # Define and fit path models using lavaan package
##' mod1 <- 'V3 ~ V1 + V2
##'   V1 ~~ 0*V2'
##' mod2 <- 'V3 ~ V1
##'   V2 ~ V3'
##' mod3 <- 'V3 ~ V2
##'   V2 ~ V1'
##'
##' fit1 <- sem(mod1, sample.cov=temp_mat, sample.nobs=500)
##' fit2 <- sem(mod2, sample.cov=temp_mat, sample.nobs=500)
##' fit3 <- sem(mod3, sample.cov=temp_mat, sample.nobs=500)
##'
##' ## obtain fitting propensity for each model
##' fp1 <- fitprop(fit1, reps = 1000, onlypos = TRUE)
##' fp2 <- fitprop(fit2, reps = 1000, onlypos = TRUE)
##' fp3 <- fitprop(fit3, reps = 1000, onlypos = TRUE)
##'
##' ## summarize 1 model's fitting propensity
##' summary(fp1, fit.measures = c("srmr","logl","cfi"),
##'         lower.tail = c(TRUE, FALSE, FALSE), NML = TRUE, UIF = TRUE)
##'
##' ## compare fitting propensity of a set of models
##' compareFitProp(fp1, fp2, fp3, fit.measures = c("srmr","logl","cfi"),
##'                lower.tail = c(TRUE, FALSE, FALSE), conf = .95,
##'                ## compare fit, adjusted for fitting propensity:
##'                NML = TRUE, UIF = TRUE)
##'
##'
##' ## Factor Models
##'
##' HS.model <- ' visual  =~ x1 + x2 + x3
##'               textual =~ x4 + x5 + x6
##'               speed   =~ x7 + x8 + x9 '
##' mod.hi <- ' g =~ visual + textual + speed '
##' mod.bi <- paste0('g =~ x', 1:9)
##' mod.par <- c(paste0('g =~ 1*x', 1:9), paste0('x', 1:9, ' ~~ evar*x', 1:9))
##'
##' fit.hi <- cfa(c(HS.model, mod.hi), data = HolzingerSwineford1939,
##'               std.lv = TRUE, estimator = "uls")
##' fit.bi <- cfa(c(HS.model, mod.bi), data = HolzingerSwineford1939,
##'               orthogonal = TRUE, std.lv = TRUE, estimator = "uls")
##' fit.uni <- cfa(mod.bi, data = HolzingerSwineford1939,
##'                std.lv = TRUE, estimator = "uls")
##' fit.par <- cfa(mod.par, data = HolzingerSwineford1939, estimator = "uls")
##'
fp.hi  <- fitprop(fit.hi, reps = 100, onlypos = TRUE,
                  parallel = "multicore", ncpus = 10,
                  baseline.model = fit.par)
fp.bi  <- fitprop(fit.bi, reps = 100, onlypos = TRUE,
                  parallel = "multicore", ncpus = 10,
                  baseline.model = fit.par)
fp.uni <- fitprop(fit.uni, reps = 100, onlypos = TRUE,
                  parallel = "multicore", ncpus = 10,
                  baseline.model = fit.par)
##' summary(fp.hi)
##' ## can name the models passed to ...
##' compareFitProp(Higher.Order = fp.hi, Bifactor = fp.bi, Unidimensinoal = fp.uni,
##'                fit.measures = c("srmr","logl","cfi"),
##'                lower.tail = c(TRUE, FALSE, FALSE), conf = .95,
##'                NML = TRUE, UIF = TRUE)

##' @importFrom stats setNames
##' @export
fitprop <- function(object,
                    use.FitProp = NULL, # option to provide existing
                    reps = 1000,
                    seed = 1234,
                    ...,
                    rmethod = c("onion","mcmc","clustergen"),
                    onlypos = FALSE,
                    mcmc.args = list(),
                    clustergen.args = list(),
                    baseline.model = NULL) {

  if (!inherits(object, "lavaan")) {
    stop("Input `object=' must be a fitted lavaan model")
  }

  if (!is.null(baseline.model))   if (!inherits(baseline.model, "lavaan")) {
    stop("Input `baseline.model=' must be a fitted lavaan model")
  }
  ## Check for multigroup or multilevel models
  #FIXME: This should only be necessary when using built-in genmat() options.
  #       Users could also supply their own dataFunction=
  if (lavInspect(object,"ngroups") > 1 || lavInspect(object,"nlevels") > 1) {
    stop("Only single-group, single-level models are currently supported")
  }

  ## Use existing FitProp data?
  #FIXME: error for any nonconverged reps (lavInspect can't return data)
  if (!is.null(use.FitProp)) {
    ## recycle using new model
    lavListCall <- use.FitProp@lavListCall
    lavListCall$model <- parTable(object)
    lavListCall$dataFunction <- NULL
    lavListCall$dataFunction.args <- NULL
    lavListCall$dataList <- use.FitProp@ranDataList#[!sapply(use.FitProp@ranDataList, is.null)]
    baseline.model <- use.FitProp@lavobjects$baseline.model

    ## otherwise, go through steps below
    dots <- list(FUN = use.FitProp@call$FUN)
    isCustom <- use.FitProp@custom

  } else {
    isCustom <- FALSE
    rmethod <- as.character(rmethod)[1L]

    ## Extract variable names
    vnames <- lavNames(object)
    ## number of variables
    d <- length(vnames)

    dots <- list(...)

    ## seed for generation used in parallel?
    if (!is.null(dots$parallel)) if (dots$parallel[1] %in% c("multicore","snow")) {
      if (is.null(dots$iseed)) dots$iseed <- as.integer(seed)
    } else set.seed(seed)

    ## fit model to random patterns using lavaanList()
    lavListCall <- object@call
    lavListCall[[1]] <- lavaan::lavaanList
    lavListCall[c("data","sample.cov","sample.mean","sample.th","sample.nobs")] <- NULL
    lavListCall$ndat <- reps
    lavListCall <- c(lavListCall, dots)
    ## did the user provide a custom data-generating function?
    if (is.null(lavListCall$dataFunction)) {
      ## check arguments for default method
      if (rmethod == "onion") {
        control <- onion.args.default(d, clustergen.args)
      } else if (rmethod == "mcmc") {
        control <- mcmc.args.default(1, d, mcmc.args)
      } else if (rmethod == "clustergen") {
        control <- clustergen.args.default(d, clustergen.args)
      }
      ## save arguments for data-generator
      genMatArgs <- list(rmethod = rmethod, control = control, onlypos = onlypos,
                         # nmat = 1L, # already the default
                         sample.nobs = lavInspect(object, "nobs"),
                         empirical = TRUE, varnames = vnames)
      ## save canned function and arguments in the call
      lavListCall$dataFunction <- ockhamSEM::genmat #FIXME? ockhamSEM::genmat
      lavListCall$dataFunction.args <- genMatArgs
    }
  }


  ## Function to get requested fitMeasures() output for FitProp object, along
  ## with info about convergence etc.
  ## NOTE: Need "lavaan::" to allow for parallel computations.
  .getOutput. <- function(obj) {
    converged <- lavaan::lavInspect(obj, "converged")
    if (converged) {
      se <- lavaan::parTable(obj)$se
      se.test <- all(!is.na(se)) & all(se >= 0) & any(se != 0)
      if (lavaan::lavInspect(obj, "ngroups") == 1L && lavaan::lavInspect(obj, "nlevels") == 1L) {
        Heywood.lv <- det(lavaan::lavInspect(obj, "cov.lv")) <= 0
        Heywood.ov <- det(lavaan::lavInspect(obj, "theta")) <= 0
      } else {
        ## multiple groups/levels, whenever the package allows it
        Heywood.lv <- !all(sapply(lavaan::lavInspect(obj, "cov.lv"), det) > 0)
        Heywood.ov <- !all(sapply(lavaan::lavInspect(obj, "theta"), det) > 0)
      }
    } else {
      se.test <- Heywood.lv <- Heywood.ov <- NA
    }
    list(ranCov = lavaan::lavInspect(obj, "sampstat"),
         ranDat = lavaan::lavInspect(obj, "data"),
         fitMeas = lavaan::fitMeasures(obj, baseline.model = baseline.model),
         converged = converged, SE = se.test,
         Heywood.lv = Heywood.lv, Heywood.ov = Heywood.ov)
  }

  ## did the user provide a custom output-extracting function?
  lavListCall$FUN <- if (is.null(dots$FUN)) .getOutput. else function(obj) {
    temp1 <- .getOutput.(obj)
    temp2 <- dots$FUN(obj)
    if (!is.list(temp2)) temp2 <- list(userFUN1 = temp2)
    if (is.null(names(temp2))) names(temp2) <- paste0("userFUN", 1:length(temp2))
    duplicatedNames <- which(sapply(names(temp2), function(x) {
      x %in% c("ranCov","ranDat","fitMeas",
               "converged","SE","Heywood.lv","Heywood.ov")
    }))
    for (i in duplicatedNames) names(temp2)[i] <- paste0("userFUN", i)
    c(temp1, temp2)
  }
  fit <- eval(as.call(lavListCall))

  ## assign class and add new slots
  fit <- as(fit, "FitProp")
  fit@call <- match.call() # overwrites lavaanList's default @call
  fit@call$reps <- reps
  fit@call$seed <- seed
  fit@call$rmethod <- rmethod
  fit@call$onlypos <- onlypos
  fit@lavListCall <- lavListCall # save lavaanList() call
  fit@custom <- isCustom
  ## Store requested output per random data pattern
  fit@momentList  <- lapply(fit@funList, "[[", i = "ranCov")
  fit@ranDataList <- lapply(fit@funList, "[[", i = "ranDat")
  obsFitMeas <- fitMeasures(object, baseline.model = baseline.model)
  fit@obsFitMeas  <- setNames(as.numeric(obsFitMeas), names(obsFitMeas))
  fit@fitMeasDist <- lapply(fit@funList, "[[", i = "fitMeas")
  convList <- lapply(fit@funList, "[", i = c("converged","SE",
                                             "Heywood.lv","Heywood.ov"))
  nonConv <- which(sapply(convList, is.null))
  if (length(nonConv)) for (i in nonConv) {
    convList[[i]] <- list(converged = FALSE, SE = NA, Heywood.lv = NA, Heywood.ov = NA)
  }
  fit@convergence <- lapply(convList, function(x) do.call(c, x))

  ## to catch problems early?
  conv <- which(sapply(fit@convergence, "[", i = "converged"))
  if (!length(conv)) warning('The model did not converge for any random data sets.')

  ## keep any remaining @funList slots (when users supply custom FUN=)
  funNames <- names(fit@funList[[1]])
  keepIndex <- which(!sapply(funNames, function(x) {
    x %in% c("ranCov","ranDat","fitMeas",
             "converged","SE","Heywood.lv","Heywood.ov")
  }))
  if (length(keepIndex)) {
    fit@funList <- lapply(fit@funList, "[", i = keepIndex)
    if (length(keepIndex) > 1L) {
      keepNames <- funNames[keepIndex]
      noNames <- which(keepNames == "")
      for (i in seq_along(noNames)) keepNames[ noNames[i] ] <- paste0("userFUN", i)
      fit@funList <- lapply(fit@funList, "names<-", value = keepNames)
    }
  } else fit@funList <- list()

  ## save original lavaan object(s)
  fit@lavobjects <- list(target = object, baseline.model = baseline.model)

  fit
}



## ----------------
## Hidden Functions
## ----------------

## Function to compare propensity distributions for 1 pair of models.
## Returns a list (one per fit measure), although compareFitProp() above only
## calls it for one fit.measure at a time.
##' @importFrom stats ks.test
##' @importFrom effsize cliff.delta cohen.d
x.vs.y <- function(x, y, sameReps, useReps, fit.measures, lower.tail = TRUE) {
  ## check class
  stopifnot(inherits(x, "FitProp") || inherits(y, "FitProp"))

  ## match lower.tail to fit.measures
  if (length(lower.tail) == 1L) lower.tail <- rep(lower.tail, length(fit.measures))
  if (length(lower.tail) != length(fit.measures))
    stop('lower.tail should be as long as the specified number of fit.measures=,',
         ' or one logical value to be repeated across fit measures.')
  names(lower.tail) <- fit.measures

  ## return nested list of pairwise comparisons (within fit.measures=)
  sapply(fit.measures, function(fm) {
    xf <- sapply(x@fitMeasDist[useReps], function(x) x[[fm]] )
    yf <- sapply(y@fitMeasDist[useReps], function(x) x[[fm]] )

    ##        when there are ties
    list(KS = suppressWarnings(ks.test(xf, yf)[c("statistic","p.value")]),
         Cohen.d = cohen.d(xf, yf, paired = sameReps,
                           na.rm = TRUE)[c("estimate","conf.int")],
         Hedges.g = cohen.d(xf, yf, hedges.correction = TRUE,
                            na.rm = TRUE)[c("estimate","conf.int")],
         Cliff.d = cliff.delta(xf, yf)[c("estimate","conf.int")])
  }, simplify = FALSE)
}


## function to weed out replications that could not:
## - converge ("no.conv"),
## - compute SEs ("no.se")
## - yield a proper solution ("no.npd")
## Also accepts integers to weed out particular replications, although that
## might not come up in this context.
chooseReps <- function(object, omit.reps = c("no.conv","no.se")) {
  useReps <- rep(TRUE, length(object@call$reps))
  if ("no.conv" %in% omit.reps) useReps <- sapply(object@convergence, "[[", i = "converged")
  if ("no.se" %in% omit.reps) useReps <- useReps & sapply(object@convergence, "[[", i = "SE")
  if ("no.npd" %in% omit.reps) {
    Heywood.lv <- sapply(object@convergence, "[[", i = "Heywood.lv")
    Heywood.ov <- sapply(object@convergence, "[[", i = "Heywood.ov")
    useReps <- useReps & !(Heywood.lv | Heywood.ov)
  }
  ## custom removal by replication number (leftover from lavaan.mi code)
  rm.imps <- omit.reps[ which(omit.reps %in% 1:length(useReps)) ]
  if (length(rm.imps)) useReps[as.numeric(rm.imps)] <- FALSE

  useReps
}





