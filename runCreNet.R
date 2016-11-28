#!/usr/bin/env Rscript
##
## Runner Script for creNet
###############################################################################

suppressPackageStartupMessages(library("optparse"))

## Parsing the commandline input
option_list = list(
  make_option(c("-m", "--method"), default = "cv",
              help="method to be applied: one of cv, test or novel. Default Value: cv"),
  make_option(c("-n", "--model"),default = 'creNet',
              help="one of creNet or lasso. Default Value: creNet"),
  make_option(c("-k", "--risk"),default = 'auc',
              help="one of auc or ll. Default Value: auc"),
  make_option(c("-u", "--uids"),default = NULL,
              help="List of regulators to consider (comma separated string of uids). Default Value: NULL"),
  make_option(c("-e", "--entries"),
              help="path to the KB entries File"),
  make_option(c("-r", "--relations"),
              help="path to KB relations File"),
  make_option(c("-d", "--traindata"),
              help="path to training data"),
  make_option(c("-t", "--testdata"), default = NULL,
              help="path to testing data (required for method test and novel)"),
  make_option(c("-a", "--nalphas"), default = 1,
              help="number of alphas to try: Default Value: 1 corresponding to alpha = 0 for pure group lasso"),
  make_option(c("-l", "--nlambdas"), default = 100,
              help="number of lambdas to try: Default Value: 100"),
  make_option(c("-w", "--weight"), default = 'cre',
              help="weight type. 1) cre, 2) log.cre, 3) sqrt, 4) both, 5) none"),
  make_option(c("-f", "--filter"), default = TRUE,
              help="apply cre filter. Default value: TRUE"),
  make_option(c("-R", "--RNAseq"), default = FALSE,
              help="Data is fpkm gene values. Default Value: FALSE"),
  make_option(c("-L", "--LOOCV"), default = FALSE,
              help="Perform loocv. Default Value: FALSE"),
  make_option(c("-i", "--iternum"), default = 10,
              help="Number of iteration to reduce CV variability and get confidence intervals. Default: 10"),
  make_option(c("-c", "--crecutoff"), default = 0.01,
              help="Cutoff for cre prior pvalues: Default Value: 0.01"),
  make_option(c("-p", "--pvaluecutoff"), default = 0.05,
              help="cutoff for DEGs pvalue: Default Value: 0.05"),
  make_option(c("-s", "--standardize"), default = 'all',
              help="standardization method: 1) self, 2) all, 3) train. Default Value: all"),
  make_option(c("-v", "--verbose"), default = TRUE,
              help="verbose: Default Value: TRUE"),
  make_option(c("-o", "--outfile"), default = NULL,
              help="path to the output prediction results file (optional)"),
  make_option(c("-z", "--groupsfile"), default = NULL,
              help="path to the output nonzero groups files (optional)"))

## Parsing the input argument
opt <- parse_args(OptionParser(option_list=option_list))

method          <- opt$method
model           <- opt$model
measure         <- opt$risk
uids            <- opt$uids
ents.file       <- opt$entries
rels.file       <- opt$relations
data.train.file <- opt$traindata
data.test.file  <- opt$testdata
nalphas         <- as.integer(opt$nalphas)
nlam            <- as.integer(opt$nlambdas)
type.weight     <- opt$weight
filter          <- as.logical(opt$filter)
RNAseq          <- as.logical(opt$RNAseq)
LOOCV           <- as.logical(opt$LOOCV)
num.iter        <- opt$iternum
cre.sig         <- opt$crecutoff
de.sig          <- opt$pvaluecutoff
verbose         <- as.logical(opt$verbose)
standardize     <- opt$standardize
output.file     <- opt$outfile
groups.file     <- opt$groupsfile

## Required libraires
require(creNet)

if(method == 'cv'){
  run.ms        <- FALSE
  run.nested.cv <- TRUE
  x.test <- NULL
  y.test <- NULL
}else if(method %in% c('test', 'novel')){
  run.ms        <- TRUE
  run.nested.cv <- FALSE
}else{
  cat('\n please provide a valid method -m \n')
  quit()
}

if(!is.null(uids)){
  uids <- strsplit(uids, split = ",")[[1]]
}


if(verbose)
  cat('\n Prepare KB and Data \n')

if(model == 'lasso'){
  isLasso = TRUE
}else{
  isLasso = FALSE
}
L <- processDataAndKB(ents.file, rels.file, data.train.file, 
                      data.test.file=data.test.file, verbose = FALSE, 
                      uids = uids, isLasso = isLasso)
ents <- L$ents
rels <- L$rels
x.train <- L$x.train
y.train <- L$y.train
x.test  <- L$x.test
y.test  <- L$y.test

if(LOOCV){
  nfold = length(y.train) - 1
}else{
  nfold = 4
}

if(RNAseq){
  x.train <- x.train + 0.25
  x.train <- log(x.train)
  if(!is.null(x.test)){
    x.test <- x.test + 0.25
    x.test <- log(x.test)
  }
}
if(model == 'lasso'){
  
  if(run.ms){
    if(verbose)
      cat(paste('\n Running Model', model, 'and method', method, '\n'))
    fit <- cvGlmnet(x.train, y.train, x.test, y.test, num.iter = num.iter, 
                    nfold = (nfold+1), alpha = 1)
    pred.probs <- fit$test.probs
    
    if(method == 'test'){
      L <- reportResults(pred.probs,y.test, 0.5, 'equal prior', verbose = verbose)
      L <- data.frame(test = rownames(L), L, stringsAsFactors = F)
    }else if(method == 'novel'){
      L <- reportNovelResults(pred.probs, 0.5, 'equal prior', verbose = verbose)
      L <- data.frame(y = L, probs = pred.probs, stringsAsFactors = F)
    }
   
     
    if(!is.null(output.file)){
      write.table(L, output.file, sep = '\t', col.names = T, row.names = F, quote = F)
    }
    
    nonzero.genes  <- fit$nonzero.genes
    nonzero.coeffs <- fit$nonzero.coeffs
    
    if(!is.null(uids)){
      ents.mRNA = ents[ents$type == 'mRNA', ]
      
      L <- data.frame(genes = ents.mRNA[match(nonzero.genes, ents.mRNA$id),], 
                      coeffs = nonzero.coeffs, stringsAsFactors = F)
      
      LL <- aggregate(L$coeffs, by = list(L$genes.uid),mean)
      L <- cbind(L[match(LL[,1], L$genes.uid), c(1,2,3,4)], LL[,2])
      colnames(L)[5] <- 'coeffs'
    }else{
      library(org.Hs.eg.db)
      library(dplyr)
      x <- org.Hs.egALIAS2EG
      mapped_genes <- mappedkeys(x)
      xx = as.list(x[mapped_genes])
      MapNam <-  data.frame(symbols = names(unlist(xx)), 
                            entrez = unlist(xx), stringsAsFactors = F)
      xx <- merge(nonzero.genes, MapNam, by.x = 1, by.y = 2)
      
      L <- data.frame(entrez = nonzero.genes,  
                      coeffs = nonzero.coeffs, stringsAsFactors = F)
      L <- merge(xx, L, by.x = 1, by.y = 1)
    }
    
    if(!is.null(groups.file)){
      write.table(L, groups.file, sep = '\t', col.names = T, row.names = F, quote = F)
    }
    
    if(verbose){
      cat('\n significant hypothesis\n')
      print(L)
    }
    
  }else{
    if(verbose)
      cat(paste('\n Running Model', model, 'and method', method, '\n'))
    
    cat('\n running nested cre cv')
    Obj <- nested.cvGlmnet(x.train, y.train, num.iter = num.iter, nfold = nfold, verbose = verbose)
    pred.train <- Obj$pred
    outer.indecies <- Obj$outer.indecies
    
    ## equal prior
    L <- nested.reportResults(pred.train,y.train, outer.indecies, 0.5, 'equal prior', verbose = verbose)
    L <- data.frame(test = rownames(L), L, stringsAsFactors = F)
    
    #if(verbose){
    #  print(nested.accuracy(pred.train, y.train, outer.indecies, cf = 0.5))
    #}
    
    if(!is.null(output.file)){
      write.table(L, output.file, sep = '\t', col.names = T, row.names = F, quote = F)
    }

    nonzero.genes  <- Obj$nonzero.genes
    nonzero.coeffs <- Obj$nonzero.coeffs
    
    if(!is.null(uids)){
      
      ents.mRNA = ents[ents$type == 'mRNA', ]
      
      L <- data.frame(genes = ents.mRNA[match(nonzero.genes, ents.mRNA$id),], 
                      coeffs = nonzero.coeffs[match(nonzero.genes, names(nonzero.coeffs))], stringsAsFactors = F)
      
      LL <- aggregate(L$coeffs, by = list(L$genes.uid),mean)
      L <- cbind(L[match(LL[,1], L$genes.uid), c(1,2,3,4)], LL[,2])
      colnames(L)[5] <- 'coeffs'
    }else{
      library(org.Hs.eg.db)
      library(dplyr)
      x <- org.Hs.egALIAS2EG
      mapped_genes <- mappedkeys(x)
      xx = as.list(x[mapped_genes])
      MapNam <-  data.frame(symbols = names(unlist(xx)), 
                            entrez = unlist(xx), stringsAsFactors = F)
      xx <- merge(unique(nonzero.genes), MapNam, by.x = 1, by.y = 2)
      
      L <- data.frame(entrez = nonzero.genes,  
                      coeffs = nonzero.coeffs, stringsAsFactors = F)
      L <- merge(xx, L, by.x = 1, by.y = 1)
      L = L %>% group_by(x,symbols) %>% summarise(mean_coeff = mean(coeffs))
      colnames(L) = c('entrez', 'name', 'coeff')
      L = L[!duplicated(cbind(L$entrez, L$coeff)),]
    }
        
    if(!is.null(groups.file)){
      write.table(L, groups.file, sep = '\t', col.names = T, row.names = F, quote = F)
    }
    
    if(verbose){
      cat('\n significant hypothesis\n')
      print(L)
    }
    
  }
}else{
  alphas <- seq(0, 1, length.out = nalphas)
  out.tab <- data.frame(matrix(-1, nrow = 7, ncol = 7), stringsAsFactors = F)
  colnames(out.tab) <- c('test', 'ep', 'ep', 'dist', 'dist', 'ba', 'ba')
  out.tab[1,1] <- 'stats'
  out.tab[1, 2:7] = rep(c('mean', 'sd'), 3)
  
  if(type.weight %in% c('none', 'sqrt')){
    filter <- FALSE
  }
  
  if(!is.null(uids)){
    type.weight <- "none"
    filter <- FALSE
  }
  
  if(run.ms){
    if(verbose)
      cat(paste('\n Running Model', model, 'and method', method, '\n'))
    
    cat('\n running models selection and prediction \n')
    
    if(verbose & is.null(uids))
      cat('\n running CRE \n')
    CF <- creFilter(ents, rels, x.train, y.train, x.test, y.test, cre.sig = cre.sig, 
                    de.sig = de.sig, RNAseq = RNAseq, 
                    type.weight = type.weight, filter = filter, verbose = FALSE)
    
    
    slice.train <- CF$slice.train
    slice.test  <- CF$slice.test
    
    groups <- CF$groups
    uid.groups <- CF$uid.groups
    sig.hyps <- CF$sig.hyps
    weights <- CF$weights
    child.uid <- CF$child.uid
    child.sgn <- CF$child.sgn
    cre.priors.sig <- CF$cre.priors.sig
    
    ## Training data
    data.train <- list(x=slice.train,y=y.train)
    
    
    ## Select alpha and lambda by 5-fold cross-validation
    fit.cv <- cvSGL(data.train, index=groups, weights=weights, 
                    standardize = standardize, type="logit", 
                    num.iter = num.iter, alphas=alphas, nlam=nlam, 
                    measure=measure, ncores=1, verbose=verbose, nfold=(nfold + 1))
    
    ## Predict responses for testing data with best (alpha,lambda) values
    ## Testing data are self-standardized (use standardize="train" to use the mean & variances of training data)
    pred.test <- predict(fit.cv,newX=slice.test,standardize="self")
    pred.test = do.call(cbind, pred.test)
    
    if(method == 'test'){
      ## equal prior
      DF <- reportResults(pred.test,y.test, 0.5, 'equal prior', verbose = verbose)
      
      out.tab[2:7, 2:3] <- DF[,1:2]
      out.tab[2:7, 1] <- rownames(DF)
      
      ## best ROC threshold dist adjustment
      cf.dist <- unlist(lapply(fit.cv$fit, function(x) x$opt.thresh.dist))
      DF <- reportResults(pred.test,y.test, cf.dist , 'best ROC threshold dist adjustment', 
                          verbose = verbose)
      
      out.tab[2:7, 4:5] <- DF[,1:2]
      
      ## best ROC threshold F1 adjustment
      #cf.f1 <- unlist(lapply(fit.cv$fit, function(x) x$opt.thresh.f1))
      #reportResults(pred.test,y.test, cf.f1, 'best ROC threshold F1 adjustment')
      
      ## best ROC threshold ba adjustment
      cf.ba <- unlist(lapply(fit.cv$fit, function(x) x$opt.thresh.ba))
      DF <- reportResults(pred.test,y.test, cf.ba, 'best ROC threshold ba adjustment', verbose = verbose)
      out.tab[2:7, 6:7] <- DF[,1:2]
      
      if(!is.null(output.file)){
        write.table(out.tab, output.file, sep = '\t', quote = F, col.names = T, row.names = F)
      }
      
    }else{
      ## equal prior
      #labs <- reportNovelResults(pred.test, 0.5, 'equal prior', verbose = TRUE)
      #cat(labs)
      ## best ROC threshold ba adjustment
      cf.ba <- unlist(lapply(fit.cv$fit, function(x) x$opt.thresh.ba))
      labs <- reportNovelResults(pred.test, cf.ba, 'best ROC threshold ba adjustment', verbose = verbose)
      if(verbose){
        cat('\n predicted labels: \n')
        cat(labs)
        cat('\n')
      }
      
      L <- data.frame(y = labs, probs = pred.test, stringsAsFactors = F)
      if(!is.null(output.file)){
        write.table(L, output.file, sep = '\t', quote = F, col.names = T, row.names = F)
      }
      
    }
    
    best.lambda <- fit.cv$best.lambda
    best.alpha <- fit.cv$best.alpha
    best.beta <- lapply(fit.cv$fit, function(x) x$beta)
    
    nonzero.genes <- unique(unlist(lapply(best.beta, function(x) which(x != 0))))
    nonzero.groups <- ents[which(ents$uid %in% unique(uid.groups[unique(nonzero.genes)])),]
    coeffs <- lapply(best.beta, function(x) x[nonzero.genes])
    coeffs.mean <- apply(do.call(cbind, coeffs), 1, mean)
    coeffs.sd <- apply(do.call(cbind, coeffs), 1, sd)
    ind.sig.coeffs <- which(abs(coeffs.mean) >= quantile(abs(coeffs.mean), prob = 0.75))
    sig.coeffs = coeffs.mean[ind.sig.coeffs]
    nonzero.genes <- nonzero.genes[ind.sig.coeffs]
    nonzero.genes <- data.frame(uid = unlist(child.uid)[nonzero.genes], beta = sig.coeffs, stringsAsFactors = F)
    nonzero.genes <- merge(ents,nonzero.genes, by.x = 1, by.y = 1) 
    
    if(verbose){
      cat('best alpha\n')
      print(best.alpha)
      cat("\n")
      
      cat('best lambda\n')
      print(best.lambda)
      cat("\n")
      
      #nonzero.genes = which(best.beta[1:ncol(CF$slice.train)] != 0)
      
      if(method == 'test'){
        print(out.tab)
      }
      cat('non-zero groups\n')
      print(nonzero.groups)
      
      cat('non-zero genes\n')
      print(nonzero.genes)
    }
    nonzero.groups = data.frame(nonzero.groups, beta = NA, stringsAsFactors = F)
    group.tab = data.frame(rbind(nonzero.groups, nonzero.genes), stringsAsFactors = F)
    if(!is.null(groups.file))
      write.table(group.tab, groups.file, sep = '\t', quote = F, col.names = T, row.names = F)
  }else{
    if(verbose)
      cat(paste('\n Running Model', model, 'and method', method, '\n'))
    
    cat('\n running nested cre cv')
    Obj <- nested.cvSGL(ents, rels, x.train, y.train, type = "logit", alphas=alphas, LOOCV=LOOCV, 
                        nlam=nlam, num.iter = num.iter,RNAseq=RNAseq,
                        type.weight = type.weight, filter = filter,
                        standardize=standardize, measure = "auc", 
                        verbose = verbose,
                        nfold=nfold, ncores=1)
    
    pred.train     <- Obj$pred
    cf.dist        <- Obj$opt.thresh.dist
    cf.f1          <- Obj$opt.thresh.f1
    cf.ba          <- Obj$opt.thresh.ba
    outer.indecies <- Obj$outer.indecies
    ##sort(unique(Obj$sig.hyps))
    
    ## equal prior
    ##nested.accuracy(pred.train, y.train, outer.indecies, cf = 0.5)
    DF <- nested.reportResults(pred.train,y.train, outer.indecies, 0.5, 'equal prior', verbose = verbose)
    
    out.tab[2:7, 2:3] <- DF[,1:2]
    out.tab[2:7, 1] <- rownames(DF)
    
    ## dist prior
    ##nested.accuracy(pred.train, y.train, outer.indecies, cf = cf.dist)
    DF <- nested.reportResults(pred.train,y.train, outer.indecies, cf.dist, 'dist prior', verbose = verbose)
    
    out.tab[2:7, 4:5] <- DF[,1:2]
    
    ## f1 prior
    #nested.accuracy(pred.train, y.train, outer.indecies, cf = cf.f1)
    #nested.reportResults(pred.train,y.train, outer.indecies, cf.f1, 'f1 prior')
    
    ## ba prior
    ## nested.accuracy(pred.train, y.train, outer.indecies, cf = cf.ba)
    DF <- nested.reportResults(pred.train,y.train, outer.indecies, cf.ba, 'ba prior', verbose = verbose)
    
    out.tab[2:7, 6:7] <- DF[,1:2]
    
    if(!is.null(output.file))
      write.table(out.tab, output.file, sep = '\t', quote = F, col.names = T, row.names = F)
    
    nonzero.genes  <- Obj$nonzero.genes
    nonzero.coeffs <- Obj$nonzero.coeffs
    
    ents.mRNA = ents[ents$type == 'mRNA', ]
    
    L <- data.frame(genes = ents.mRNA[match(nonzero.genes, ents.mRNA$uid),], 
                    coeffs = nonzero.coeffs, stringsAsFactors = F)
    LL <- aggregate(L$coeffs, by = list(L$genes.uid),mean)
    L <- cbind(L[match(LL[,1], L$genes.uid), c(1,2,3,4)], LL[,2])
    colnames(L)[5] <- 'coeffs'
    L <- L[which(abs(L$coeffs) >= quantile(abs(L$coeffs), prob = 0.75)),]
    nonzero.groups <- cbind(ents[match(Obj$sig.hyps, ents$name),], NA)
    colnames(nonzero.groups) <- colnames(L)
    L <- rbind(nonzero.groups, L)
    
    if(!is.null(groups.file)){
      write.table(L, groups.file, sep = '\t', col.names = T, row.names = F, quote = F)
    }
    
    if(verbose){
      cat('\n significant hypothesis\n')
      print(L)
    }
  }
}
