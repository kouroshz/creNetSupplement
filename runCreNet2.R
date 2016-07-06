#!/usr/bin/env Rscript
##
## Runner Script for the new creNet2
###############################################################################

suppressPackageStartupMessages(library("optparse"))

## Parsing the commandline input
option_list = list(
  make_option(c("-m", "--method"), default = "cv",
              help="method to be applied: Either nested CV on training data (cv) or
              model selection and test on test data (test): Default Value: cv"),
  make_option(c("-n", "--model"),default = 'group.lasso',
              help="one of group.lasso or lasso"),
  make_option(c("-e", "--entries"),
              help="path to the KB entries File"),
  make_option(c("-r", "--relations"),
              help="path to KB relations File"),
  make_option(c("-d", "--traindata"),
              help="path to training data"),
  make_option(c("-t", "--testdata"), default = NULL,
              help="path to testing data"),
  make_option(c("-a", "--nalphas"), default = 1,
              help="number of alphas to try: Default Value: 1 corresponding to alpha = 0 for pure group lasso"),
  make_option(c("-l", "--nlambdas"), default = 100,
              help="number of lambdas to try: Default Value: 100"),
  make_option(c("-w", "--weight"), default = 'cre',
              help="weight type. 1) cre, 2) log.cre, 3) sqrt, 4) both, 5) none"),
  make_option(c("-f", "--filter"), default = TRUE,
              help="apply cre filter. Default value: TRUE"),
  make_option(c("-i", "--iternum"), default = 10,
              help="Number of iteration to reduce CV vriability and get confidence intervals. Default: 10"),
  make_option(c("-c", "--crecutoff"), default = 0.01,
              help="Cutoff for cre prior pvalues: Default Value: 0.01"),
  make_option(c("-p", "--pvaluecutoff"), default = 0.05,
              help="cutoff for DEGs pvalue: Default Value: 0.05"),
  make_option(c("-s", "--standardize"), default = 'all',
              help="standardization method: 1) self, 2) all, 3) train. Default Value: all"),
  make_option(c("-v", "--verbose"), default = TRUE,
              help="verbose: Default Value: TRUE"),
  make_option(c("-o", "--output"), default = NULL,
              help="path to the output File"))

## Parsing the input argument
opt <- parse_args(OptionParser(option_list=option_list))

method          <- opt$method
model           <- opt$model
ents.file       <- opt$entries
rels.file       <- opt$relations
data.train.file <- opt$traindata
data.test.file  <- opt$testdata
nalphas         <- as.integer(opt$nalphas)
nlam            <- as.integer(opt$nlambdas)
type.weight     <- opt$weight
filter          <- opt$filter
num.iter        <- opt$iternum
cre.sig         <- opt$crecutoff
de.sig          <- opt$pvaluecutoff
verbose         <- as.logical(opt$verbose)
standardize     <- opt$standardize
output.file     <- opt$output

## Required libraires
require(creNet2)

if(method == 'cv'){
  run.ms        <- FALSE
  run.nested.cv <- TRUE
  x.test <- NULL
  y.test <- NULL
}else if(method == 'test'){
  run.ms        <- TRUE
  run.nested.cv <- FALSE
}

if(verbose)
  cat('\n Prepare KB and Data \n')

L <- processDataAndKB(ents.file, rels.file, data.train.file, data.test.file=data.test.file, verbose = FALSE)
ents <- L$ents
rels <- L$rels
x.train <- L$x.train
y.train <- L$y.train
x.test <- L$x.test
y.test <- L$y.test


if(model == 'lasso'){
  if(run.ms){
    if(verbose)
      cat(paste('\n Running Model', model, 'and method', method, '\n'))
    fit <- cvGlmnet(x.train, y.train, x.test, y.test, num.iter = num.iter, nfold = 5, alpha = 1)
    pred.probs <- fit$test.probs
    
    print(reportResults(pred.probs,y.test, 0.5, 'equal prior'))
    nonzero.genes  <- fit$nonzero.genes
    cat('\n significant hypothesis\n')
    ents.mRNA = ents[ents$type == 'mRNA', ]
    print(ents.mRNA[which(ents.mRNA$id %in% nonzero.genes),])
    
  }else{
    if(verbose)
      cat(paste('\n Running Model', model, 'and method', method, '\n'))
    
    cat('\n running nested cre cv')
    Obj <- nested.cvGlmnet(x.train, y.train, num.iter = num.iter, nfold = 4, verbose = verbose)
    pred.train <- Obj$pred
    outer.indecies <- Obj$outer.indecies
    nonzero.genes  <- Obj$nonzero.genes
    ## equal prior
    print(nested.accuracy(pred.train, y.train, outer.indecies, cf = 0.5))
    print(nested.reportResults(pred.train,y.train, outer.indecies, 0.5, 'equal prior'))
    
    cat('\n significant hypothesis\n')
    ents.mRNA = ents[ents$type == 'mRNA', ]
    print(ents.mRNA[which(ents.mRNA$id %in% nonzero.genes),])
  }
}else{
  alphas <- seq(0, 1, length.out = nalphas)
  
  if(run.ms){
    if(verbose)
      cat(paste('\n Running Model', model, 'and method', method, '\n'))
    
    cat('\n running models selection and prediction \n')
    
    if(type.weight %in% c('none', 'sqrt')){
      filter <- FALSE
    }
    if(verbose)
      cat('\n running CRE \n')
    CF <- creFilter(ents, rels, x.train, y.train, x.test, y.test, cre.sig = cre.sig, de.sig = de.sig,
                    type.weight = type.weight, filter = filter, verbose = FALSE)
    
    
    slice.train <- CF$slice.train
    slice.test  <- CF$slice.test
    
    groups <- CF$groups
    uid.groups <- CF$uid.groups
    sig.hyps <- CF$sig.hyps
    weights = CF$weights
    
    cre.priors.sig <- CF$cre.priors.sig
    
    ## Training data
    data.train <- list(x=slice.train,y=y.train)
    
    
    ## Select alpha and lambda by 5-fold cross-validation
    fit.cv <- cvSGL(data.train, index=groups, weights=weights, standardize = standardize, type="logit", 
                    num.iter = num.iter, alphas=alphas, nlam=nlam, measure="auc", ncores=1, verbose=verbose, nfold=5)
    
    ## Predict responses for testing data with best (alpha,lambda) values
    ## Testing data are self-standardized (use standardize="train" to use the mean & variances of training data)
    pred.test <- predict(fit.cv,newX=slice.test,standardize="self")
    pred.test = do.call(cbind, pred.test)
    ## equal prior
    reportResults(pred.test,y.test, 0.5, 'equal prior')
    
    ## best ROC threshold dist adjustment
    cf.dist <- unlist(lapply(fit.cv$fit, function(x) x$opt.thresh.dist))
    reportResults(pred.test,y.test, cf.dist , 'best ROC threshold dist adjustment')
    
    ## best ROC threshold F1 adjustment
    cf.f1 <- unlist(lapply(fit.cv$fit, function(x) x$opt.thresh.f1))
    reportResults(pred.test,y.test, cf.f1, 'best ROC threshold F1 adjustment')
    
    ## best ROC threshold ba adjustment
    cf.ba <- unlist(lapply(fit.cv$fit, function(x) x$opt.thresh.ba))
    reportResults(pred.test,y.test, cf.ba, 'best ROC threshold ba adjustment')
    
    
    best.lambda = fit.cv$best.lambda
    best.alpha = fit.cv$best.alpha
    best.beta <- lapply(fit.cv$fit, function(x) x$beta)
    nonzero.genes <- unique(unlist(lapply(best.beta, function(x) which(x != 0))))
    
    cat('best alpha\n')
    print(best.alpha)
    cat("\n")
    
    cat('best lambda\n')
    print(best.lambda)
    cat("\n")
    
    #nonzero.genes = which(best.beta[1:ncol(CF$slice.train)] != 0)
    
    cat('non-zero groups\n')
    print(ents[which(ents$uid %in% unique(uid.groups[unique(nonzero.genes)])),])
  }else{
    if(verbose)
      cat(paste('\n Running Model', model, 'and method', method, '\n'))
    
    cat('\n running nested cre cv')
    Obj <- nested.cvSGL(ents, rels, x.train, y.train, type = "logit", alphas=alphas, nlam=nlam, num.iter = num.iter,
                        type.weight = type.weight, standardize=standardize, measure = "auc", nfold=4, ncores=1)
    
    pred.train     <- Obj$pred
    cf.dist        <- Obj$opt.thresh.dist
    cf.f1          <- Obj$opt.thresh.f1
    cf.ba          <- Obj$opt.thresh.ba
    outer.indecies <- Obj$outer.indecies
    sort(unique(Obj$sig.hyps))
    
    ## equal prior
    nested.accuracy(pred.train, y.train, outer.indecies, cf = 0.5)
    nested.reportResults(pred.train,y.train, outer.indecies, 0.5, 'equal prior')
    
    ## dist prior
    nested.accuracy(pred.train, y.train, outer.indecies, cf = cf.dist)
    nested.reportResults(pred.train,y.train, outer.indecies, cf.dist, 'dist prior')
    
    ## f1 prior
    nested.accuracy(pred.train, y.train, outer.indecies, cf = cf.f1)
    nested.reportResults(pred.train,y.train, outer.indecies, cf.f1, 'f1 prior')
    
    ## ba prior
    nested.accuracy(pred.train, y.train, outer.indecies, cf = cf.ba)
    nested.reportResults(pred.train,y.train, outer.indecies, cf.ba, 'ba prior')
    
    cat('\n significant hypothesis\n')
    print(sort(unique(Obj$sig.hyps)))
  }
}
