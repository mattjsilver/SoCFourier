#^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
# Methylation array-wide analysis of the association between date of conception
# and DNA methylation, adjusted for specified covariates.
# 
# The relationship between DNA methylation (outcome) and date of conception
# (predictor) is modelled using Fourier regression with two Fourier terms: 
# sin(DoC) and cos(DoC) where DoC is date of conception measured in radians 
# (0=1st Jan; 2*pi=31st Dec). Adjustment covariates are specified by the user 
# in design matrix required as input.
#
# Fourier model with 2 Fourier terms assumes a sinusoidal relationship between 
# methylation and DoC with a single methylation maximum and minimum within a
# single year. The amplitude (distance between maximum and minimum) 
# and phase (time of year of methylation maximum/minimum) are determined by
# the estimated coefficients for each Fourier term.
#
# The significance of any 'seasonal' (date of conception) association is deter-
# mined by likelihood ratio test (LRT) comparing the full model with both 
# Fourier terms against a nested, covariates-only model
#
# See Silver et al. Environmentally sensitive hotspots in the methylome of the 
# early human embryo for further details:
# https://www.biorxiv.org/content/10.1101/777508v2
#^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^


#^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
#' run Fourier regression models and compute LRT pvals
#'
#' @param meth df methylation values (cpgs x samples)
#' @param design df model design matrix (cpgs x covs); no intercept required
#' @param n_threads number of cpu threads to use for parallel computation
#'
#' @return df with LRT chi-squared stat, pval and FDR for all tested cpgs
#'
generate_fourier_lrt_pvals = function(meth=NULL,design=NULL,n_threads=10){

  require(lmtest)
  require(lubridate)
  require(doParallel) # required for multi-processor execution of array-wide LRTs
  require(tidyverse)

  cat('\n** fitting model...\n')
  cat('design matrix covariates:\n',colnames(design),'\n')

  # covert date of conception to radians (doc.theta) in [0,2*pi]
  day.of.yr = yday(design$doc)
  design$doc.theta = day.of.yr/366*2*pi
  design$sin.doc.theta = sin(design$doc.theta)
  design$cos.doc.theta = cos(design$doc.theta)
  design$doc = NULL
  design$doc.theta = NULL

  # design matrix for full model (1 pair of Fourier terms + covariates)
  design.fourier = design

  # design matrix for baseline model (covariates only)
  design.base = design %>% 
    dplyr::select(
      -c(sin.doc.theta,cos.doc.theta)
    )
  
  cat('performing likelihood ratio tests to identify
      models with significant seasonality...\n')

  lrt.pval = function(h0,h1){
    test = lrtest(h0,h1)
    return(list(chisq=test$Chisq[2],pval=test$`Pr(>Chisq)`[2]))
  }

  cl <- makeCluster(n_threads)
  registerDoParallel(cl)

  lrt.results =

    foreach (cg = row.names(meth), .combine='c', .packages = 'lmtest') %dopar%

    {

      # compare best Fourier model against baseline
      fit.base = lm(meth[cg,]~.,data=design.base)
      fit.fourier = lm(meth[cg,]~.,data=design.fourier)
      lrt.res = lrt.pval(fit.base,fit.fourier)
      c(cg,lrt.res$chisq,lrt.res$pval)

    }
  stopCluster(cl)

  cat('lrt completed\n\n')

  lrt.results = as.data.frame(matrix(lrt.results,ncol = 3,byrow = T),stringsAsFactors = F)
  colnames(lrt.results) = c('cg','lrt.chisq','lrt.pval')
  row.names(lrt.results) = lrt.results$cg
  lrt.results$lrt.chisq = as.numeric(lrt.results$lrt.chisq)
  lrt.results$lrt.pval = as.numeric(lrt.results$lrt.pval)
  lrt.results$lrt.fdr = p.adjust(lrt.results$lrt.pval,method = 'fdr')

  return(lrt.results)
}


#^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
#' Generate fitted values from Fourier regression model over [0,2*pi].
#' Values are mean-centered
#'
#' @param model.coefs Fourier model coefficients produced by coef(lm.fit) 
#'   where lm.fit is object of class 'lm' produced by lm()
#' @param design df model design matrix (cpgs x covs) including Fourier terms;
#'   no intercept required
#'
#' @return df of fitted values (cpgs x fitted values)
#'
get_fourier_fitted_vals <- function(model.coefs = NULL, design = NULL) {

  # generate fitted M-values with intercept
  doc.theta.vals = seq(0,2*pi,0.02)
  design.fitted.vals = data.frame(
    intercept = rep(1,length(doc.theta.vals)),
    sin.doc.theta = sin(doc.theta.vals),
    cos.doc.theta = cos(doc.theta.vals)
  )

  coefs = t(model.coefs)
  fitted.vals =
    coefs[,c("`(Intercept)`","`sin(doc.theta)`","`cos(doc.theta)`")] %*%
    t(design.fitted.vals[,c('intercept','sin.doc.theta','cos.doc.theta')])

  # convert to methylation Betas and center
  fitted.vals = m2beta(fitted.vals)
  fitted.vals=as.data.frame(fitted.vals-rowMeans(fitted.vals))
  colnames(fitted.vals) = doc.theta.vals

  return(fitted.vals)

}


#^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
#' Get summary stats for fitted Fourier regression model
#'   These comprise date of methylation maximum, methylation minimum and 
#'   amplitude (difference between methylation values at maximum and minimum)
#'
#' @param fitted.vals fitted.vals object from get_fourier_fitted_vals()
#' @param lrt.results lrt.results object from generate_fourier_lrt_pvals()
#' @param soc.cpgs chr vector of cg.ids for SoC-CpGs
#' @param ctrl.cpgs chr vector of cg.ids for matched control cpgs
#' @param rand.cpgs chr vector of cg.ids for random control cpgs
#'
#' @return summary stats
#'
generate_model_fit_stats <- function(
  fitted.vals=NULL, lrt.results=NULL, ctrl.cpgs=NULL, rand.cpgs=NULL, soc.cpgs=NULL)
{
  fitted.vals.max = apply(fitted.vals,1,function(x) max(as.numeric(x[1:(length(x)-2)])))
  fitted.vals.max.id = apply(fitted.vals,1,function(x) which.max(as.numeric(x[1:(length(x)-2)])))

  fitted.vals.min = apply(fitted.vals,1,function(x) min(as.numeric(x[1:(length(x)-2)])))
  fitted.vals.min.id = apply(fitted.vals,1,function(x) which.min(as.numeric(x[1:(length(x)-2)])))

  model.stats = data.frame(
    cg = row.names(fitted.vals),
    max.doc.theta =
      colnames(fitted.vals)[fitted.vals.max.id],
    min.doc.theta =
      colnames(fitted.vals)[fitted.vals.min.id],
    amplitude = abs(fitted.vals.max - fitted.vals.min),
    stringsAsFactors = F
  )

  # add LRT pvals
  model.stats$lrt.pval = lrt.results[row.names(model.stats),'lrt.pval']
  model.stats$lrt.fdr = lrt.results[row.names(model.stats),'lrt.fdr']

  model.stats$locus = NA
  model.stats[model.stats$cg %in% soc.cpgs,'locus'] = 'soc'
  model.stats[model.stats$cg %in% ctrl.cpgs,'locus'] = 'ctrl'
  model.stats[model.stats$cg %in% rand.cpgs,'locus'] = 'rand'

  return(model.stats)
}


#^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
#' Transform methylation M-values to Beta values
#'
#' @param M df or matrix of M-values
#'
#' @return df or matrix of methylation Beta values
m2beta = function (M)
{
  beta <- 2^M/(2^M + 1)
  return(beta)
}