#!/usr/bin/env Rscript
##
## Bootstrap analysis for robustness of selected regulators
###############################################################################

suppressPackageStartupMessages(library("optparse"))

## Parsing the commandline input
option_list = list(
  make_option(c("-e", "--entries"),
              help="path to the KB entries File"),
  make_option(c("-r", "--relations"),
              help="path to KB relations File"),
  make_option(c("-d", "--traindata"),
              help="path to training data"),
  make_option(c("-a", "--nalphas"), default = 1,
              help="number of alphas to try: Default Value: 1 corresponding to alpha = 0 for pure group lasso"),
  make_option(c("-l", "--nlambdas"), default = 100,
              help="number of lambdas to try: Default Value: 100"),
  make_option(c("-w", "--weight"), default = 'cre',
              help="weight type. 1) cre, 2) log.cre, 3) sqrt, 4) both, 5) none"),
  make_option(c("-f", "--filter"), default = TRUE,
              help="apply cre filter. Default value: TRUE"),
  make_option(c("-i", "--iternum"), default = 100,
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
              help="path to the output prediction results file (optional)"))

## Parsing the input argument
opt <- parse_args(OptionParser(option_list=option_list))

ents.file       <- opt$entries
rels.file       <- opt$relations
data.train.file <- opt$traindata
nalphas         <- as.integer(opt$nalphas)
nlam            <- as.integer(opt$nlambdas)
type.weight     <- opt$weight
filter          <- opt$filter
num.iter        <- opt$iternum
cre.sig         <- opt$crecutoff
de.sig          <- opt$pvaluecutoff
verbose         <- as.logical(opt$verbose)
standardize     <- opt$standardize
output.file     <- opt$outfile

## Required libraires
require(creNet)

L <- processDataAndKB(ents.file, rels.file, data.train.file, data.test.file=NULL, 
                      verbose = FALSE, uids = NULL)
ents <- L$ents
rels <- L$rels
x.train <- L$x.train
y.train <- L$y.train
x.test <- L$x.test
y.test <- L$y.test

if(type.weight %in% c('none', 'sqrt')){
  filter <- FALSE
}

n <- length(y.train)
alphas <- seq(0, 1, length.out = nalphas)

G <- {}
for(i in 1:num.iter){
  ## sample the training data.
  ind <- unique(sample(1:n, size = n, replace = T))
  x <- x.train[ind,]
  y <- y.train[ind]
  
  cat(paste('\n boot strap sample:', i, '\n'))
  
  CF <- creFilter(ents, rels, x, y, x.test, y.test, cre.sig = cre.sig, de.sig = de.sig,
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
  data.train <- list(x=slice.train,y=y)
  
  
  ## Select alpha and lambda by 5-fold cross-validation
  fit.cv <- cvSGL(data.train, index=groups, weights=weights, standardize = standardize, 
                  type="logit", num.iter = 1, alphas=alphas, nlam=nlam, 
                  measure=measure, ncores=1, verbose=FALSE, nfold=5)
  
  best.lambda <- fit.cv$best.lambda
  best.alpha <- fit.cv$best.alpha
  best.beta <- lapply(fit.cv$fit, function(x) x$beta)
  
  nonzero.genes <- unique(unlist(lapply(best.beta, function(x) which(x != 0))))
  nonzero.groups <- ents[which(ents$uid %in% unique(uid.groups[unique(nonzero.genes)])),]

  
  cat('\n non-zero groups \n')
  print(nonzero.groups)
  
  G <- c(G,nonzero.groups$name)
}


rle(G)$lengths
rle(G)$values

nonzero.groups = data.frame(nonzero.groups = rle(G)$values, 
                            frequencies = rle(G)$lengths, stringsAsFactors = F)
if(!is.null(output.file))
  write.table(nonzero.groups, groups.file, sep = '\t', quote = F, col.names = T, 
              row.names = F)

