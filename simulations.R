### The following functions perform some statistics and post-processing:

########## Stat, pos-processing Start
## This function returns the network directly connected to a given node.
##
## Arguments:
##
## node     -    id of the node
##
## ents     -    entries of the network
##
## rels     -    relationships
##
## levels   -    boolean (include other level for proteins)
nodeNet = function(node, ents, rels, levels = F){
  node.ents = {}
  node.rels = {}
  
  ## Find the targets of the node
  targ = rels[which(rels[,2] == node),]
  node.ents = rbind(node.ents, ents[which(ents[,1] %in% targ[,3]), ])
  node.rels = rbind(node.rels, targ)
  
  ## Find the sources of the node
  src = rels[which(rels[,3] == node),]
  node.ents = rbind(node.ents, ents[which(ents[,1] %in% src[,2]), ])
  node.rels = rbind(node.rels, src)
  
  ## add the node
  node.ents = rbind(node.ents, ents[which(ents[,1] == node), ])
  
  ## consider other level
  if(levels){
    if(ents[which(ents[,1] == node), 4] == 'protein'){
      if(substr(node,nchar(node), nchar(node)) == '1'){
        node2 = paste(substr(node,1, (nchar(node) - 1)), '0', sep = '')
      }else{
        node2 = paste(substr(node,1, (nchar(node) - 1)), '1', sep = '')
      }
    }
    ## Find the targets of the node
    targ = rels[which(rels[,2] == node2),]
    node.ents = rbind(node.ents, ents[which(ents[,1] %in% targ[,3]), ])
    node.rels = rbind(node.rels, targ)
    
    ## Find the sources of the node
    src = rels[which(rels[,3] == node2),]
    node.ents = rbind(node.ents, ents[which(ents[,1] %in% src[,2]), ])
    node.rels = rbind(node.rels, src)
    
    ## add the node
    node.ents = rbind(node.ents, ents[which(ents[,1] == node2), ])
    
  }
  
  L = list(ents = unique(node.ents), rels = unique(node.rels))
  
}

## This function returns the network directly connected to a list of given node.
##
## Arguments:
##
## nodes    -    list of id of the nodes
##
## ents     -    entries of the network
##
## rels     -    relationships
##
## levels   -    boolean (include other level for proteins)
nodeNetList = function(nodes, ents, rels, levels = F){
  node.ents = {}
  node.rels = {}
  
  for(node in nodes){
    ## Find the targets of the node
    targ = rels[which(rels[,2] == node),]
    node.ents = rbind(node.ents, ents[which(ents[,1] %in% targ[,3]), ])
    node.rels = rbind(node.rels, targ)
    
    ## Find the sources of the node
    src = rels[which(rels[,3] == node),]
    node.ents = rbind(node.ents, ents[which(ents[,1] %in% src[,2]), ])
    node.rels = rbind(node.rels, src)
    
    ## add the node
    node.ents = rbind(node.ents, ents[which(ents[,1] == node), ])
    
    ## consider other level
    if(levels){
      if(ents[which(ents[,1] == node), 4] == 'protein'){
        if(substr(node,nchar(node), nchar(node)) == '1'){
          node2 = paste(substr(node,1, (nchar(node) - 1)), '0', sep = '')
        }else{
          node2 = paste(substr(node,1, (nchar(node) - 1)), '1', sep = '')
        }
      }
      ## Find the targets of the node
      targ = rels[which(rels[,2] == node2),]
      node.ents = rbind(node.ents, ents[which(ents[,1] %in% targ[,3]), ])
      node.rels = rbind(node.rels, targ)
      
      ## Find the sources of the node
      src = rels[which(rels[,3] == node2),]
      node.ents = rbind(node.ents, ents[which(ents[,1] %in% src[,2]), ])
      node.rels = rbind(node.rels, src)
      
      ## add the node
      node.ents = rbind(node.ents, ents[which(ents[,1] == node2), ])
      
    }
    
  }
  
  L = list(ents = unique(node.ents), rels = unique(node.rels))
  
}

## This function prints the stat of a node
##
## Arguments:
##
## reg      -    direction of regulation
##
## node     -    id of the node
##
## ents     -    entries of the network
##
## rels     -    relationships
##
## evidence -    gene values
node.stat = function(reg, node, ents, rels, evidence){
  
  targ = rels[which(rels[,2] == node),3]
  dirs = ifelse(rels[which(rels[,2] == node), 4] == 'increase', 1, -1)
  
  targ.vals = rep(0, length(targ))
  for(t in targ[which(targ %in% evidence[,1])]){
    targ.vals[which(targ == t)] = as.numeric(as.character(evidence[which(evidence[,1] == t),2]))
  }
  
  prediction = reg * as.numeric(dirs)
  r = prediction - targ.vals
  c = length(which(r == 0))
  z = length(which(targ.vals == 0))
  i = length(targ.vals) - (z + c)
  
  
  L = list(c = c, i = i, z = z)
  L
}

## This function prints the stat of a list of nodes
##
## Arguments:
##
## nodes      -    list of nodes
##
## node.vals  -    node values (direction of regulation)
##
## ents       -    entries of the network
##
## rels       -    relationships
##
## evidence   -    gene values

node.stat.list = function(nodes, node.vals, ents, rels, evidence){
  
  L = nodeNetList(nodes, ents, rels, levels = F)
  on.trans = merge(L$ents[which(L$ents[,4] == 'mRNA'),1, drop = F], evidence, by.x = 1, by.y = 1)
  no.pred.num = dim(evidence)[1] - dim(on.trans)[1]
  
  correct.num = 0 #explained by at least one node
  incorrect.num = 0 #wrong prediction by all of the nodes
  for(t in 1:dim(on.trans)[1]){
    tr = on.trans[t,1]
    tr.val = on.trans[t,2]
    tr.par = L$rels[which(L$rels[,3] == tr),2]
    tr.dir = ifelse(L$rels[which(L$rels[,3] == tr),4]  == 'increase', 1, -1)
    par.vals = node.vals[which(nodes %in% tr.par)]
    pred.vals = par.vals * tr.dir
    if(any(pred.vals == tr.val)){
      correct.num = correct.num + 1
    }else{
      incorrect.num = incorrect.num + 1
    }
  }
  
  L = list(correct = correct.num, incorrect = incorrect.num, zero = no.pred.num)
}


########## Stat, pos-processing End

########## Functions for generating evidence data and simulating data:

###############################################################
## Given a (one-level) network, this function simulates the Differentially expressed 
## gene values values using the
## given parameters. 
##
## Arguments:
##
## ents    -    one-level network ents (generated from a given KB using genOneLevel)
##
## rels    -    one-level network rels (generated from a given KB using genOneLevel)
##
## nodes   -    list of uid of selected hypothesis
##
## values  -    values of selected hypothesis
##
## p.c     -   probability of the edge being wrong for nonzero edges
##
## p.m     -   prior probability of the ture state of the genes being -1
##
## p.z     -   prior probability of the ture state of the genes being 0
##
## p.p     -   prior probability of the ture state of the genes being 1
##
## h.z     -   prior probability of the hypothesis being off
##
## alpha   -   false positive rate
##
## beta    -   false negative rate
simDEGs <- function(ents, rels, nodes, values, p.c, p.m, p.z, p.p, alpha, beta){
  ## Simulate values.
  mrna.ind = which(ents[,'type'] == 'mRNA')
  num.mrna = length(mrna.ind)
  sim.data = data.frame(id = ents[mrna.ind, 'uid'], val = rep(0,num.mrna), stringsAsFactors = F)
  colnames(sim.data) = c('id', 'val')
  
  ## simulate the value of each mRNA
  for(i in 1:dim(sim.data)[1]){
    i.ind = which(rels[,'trguid'] == sim.data[i,1])
    direction = ifelse(rels[i.ind,'type'] == 'increase', 1, -1)
    
    src.ind = which(nodes %in% rels[i.ind,'srcuid'])
    pred.val = sum(values[src.ind] * direction)
    
    if(pred.val >= 1){
      ph = c(p.c * p.m, p.c * p.z, 1 - p.c + p.c * p.p, 0)
    }else if(pred.val <= -1){
      ph = c(1 - p.c + p.c * p.m, p.c * p.z, p.c * p.p, 0)
    }else{ ## Conflict case with equal 1s and -1s: Flip a coin.
      z = runif(1)
      if(z > 0.5){
        ph = c(p.c * p.m, p.c * p.z, 1 - p.c + p.c * p.p, 0)
      }else{
        ph = c(1 - p.c + p.c * p.m, p.c * p.z, p.c * p.p, 0)
      }
    }
    
    pzm1 = c(1-2*beta,alpha,beta,(1/3))
    pz0  = c(beta,1-2*alpha,beta,(1/3))
    pzp1 = c(beta,alpha,1-2*beta,(1/3))
    
    sim.prob = c(sum(pzm1*ph), sum(pz0*ph), sum(pzp1*ph))
    
    
    sim.data[i,2] = sample(c(-1,0,1), 1, prob=sim.prob)
  }
  
  ## change uid to entrez id
  sim.data[,1] = ents$id[match(sim.data[,1], ents$uid)]
  
  sim.data
}


## Simulating RNA-seq counts using a negative binomila distribution
## Negative binomial captures the vairation technical and biological
## replicates. This function simulates the data as follows. A base-line
## expression mean mu0 is provided for simulating the counts of the genes
## with differential expression status 0. Next given a fold change fc
## the mean of the upregulated genes is set to fc * mu0 and the mean
## of the down regulated genes is set to (1/fc) * mu0 for the case where
## fc > 1. For the case fc < 1 the situation is reveresed. The variance
## of the negative binomial is mu + (mu^2)/r where r is the dispersion
## parameter. The r value is set to mu/3 by default. 
##
## INPUT
##
## ents    -    network entities
##
## rels    -    network relations
##
## n       -    number of positive examples
## 
## m       -    number of negative examples
##
## degs    -    list of differentially expressed genes (output of simDEGs)
##
## fc      -    fold change between the two conditions
##
## mu0     -    base-level mean (e.g., 100)
##
## r       -    dispersion parameter.
simRNAseq <- function(ents, rels, n, m, degs, fc = 2, mu0 = 3, r = NA){
  ents.mRNA <- ents[which(ents[,'type'] == 'mRNA'),]
  numers <- unlist(lapply(ents.mRNA$id, pystr_isnumeric))
  ents.mRNA <- ents.mRNA[numers,]
  x <- matrix(0, nrow = (n+m), ncol = nrow(ents.mRNA))
  if(is.na(r)){
    r <- mu0/3
  }

  if(fc > 1){
    mu1 <- fc * mu0
    mu2 <- (1 / fc) * mu0
  }else{
    mu1 <- (1 / fc) * mu0
    mu2 <- fc * mu0
  }
  
  ind.up   <- match(degs$id[degs$val == 1], ents.mRNA$id)
  ind.down <- match(degs$id[degs$val == -1], ents.mRNA$id)
  nn1 <- length(ind.up)
  nn2 <- length(ind.down)
  nn  <- nrow(ents.mRNA) - (nn1 + nn2)
  if(length(ind.up) > 0){
    x[1:n, ind.up] = matrix(rnbinom(n = n*nn1, size = r, mu = mu1), nrow = n)
    x[(n+1):(n+m), ind.up] = matrix(rnbinom(n = m*nn1, size = r, mu = mu2), nrow = m)
  }
  if(length(ind.down) > 0){
    x[1:n, ind.down] = matrix(rnbinom(n = n*nn2, size = r, mu = mu2), nrow = n)
    x[(n+1):(n+m), ind.down] = matrix(rnbinom(n = m*nn2, size = r, mu = mu1), nrow = m)
  }
  
  x[,-c(ind.up, ind.down)] = matrix(rnbinom(n = (n+m)*nn, size = r, mu = mu0), nrow = (n+m))
  
  dat <- data.frame(y = c(rep(1, n), rep(0, m)), x, stringsAsFactors = F)
  colnames(dat)[1] <- 'y'
  colnames(dat)[2:ncol(dat)] <- ents.mRNA$id
  
  return(dat)
}

## Simulating RNA-seq counts using a negative binomila distribution
## Negative binomial captures the vairation technical and biological
## replicates. This function simulates the data as follows. A base-line
## expression mean mu0 is provided for simulating the counts of the genes
## with differential expression status 0. Next given a fold change fc
## the mean of the upregulated genes is set to fc * mu0 and the mean
## of the down regulated genes is set to (1/fc) * mu0 for the case where
## fc > 1. For the case fc < 1 the situation is reveresed. The variance
## of the negative binomial is mu + (mu^2)/r where r is the dispersion
## parameter. The r value is set to mu/3 by default. 
##
## Note: in this version fc will be selected at random between fc.min
## and fc.max per pation, per gene.
##
## INPUT
##
## ents    -    network entities
##
## rels    -    network relations
##
## n       -    number of positive examples
## 
## m       -    number of negative examples
##
## degs    -    list of differentially expressed genes (output of simDEGs)
##
## fc.min  -    min fold change between the two conditions
##
## fc.max  -    max fold change between the two conditions
##
## mu0     -    base-level mean (e.g., 100)
##
## r       -    dispersion parameter.
simRNAseq2 <- function(ents, rels, n, m, degs, fc.min = 1.1, fc.max = 2, mu0 = 3, r = NA){
  ents.mRNA <- ents[which(ents[,'type'] == 'mRNA'),]
  numers <- unlist(lapply(ents.mRNA$id, pystr_isnumeric))
  ents.mRNA <- ents.mRNA[numers,]
  x <- matrix(0, nrow = (n+m), ncol = nrow(ents.mRNA))
  if(is.na(r)){
    r <- mu0/3
  }

  ind.up   <- match(degs$id[degs$val == 1], ents.mRNA$id)
  ind.down <- match(degs$id[degs$val == -1], ents.mRNA$id)
  
  nn1 <- length(ind.up)
  nn2 <- length(ind.down)
  nn  <- nrow(ents.mRNA)

  if(length(ind.up) > 0){
    for(jj in ind.up){
      for(ii in 1:n){
        fc = runif(1, fc.min, fc.max)
        mu1 <- fc * mu0
        x[ii, jj] = rnbinom(n = 1, size = r, mu = mu1)
      }
      for(ii in (n+1):(n+m)){
        fc = runif(1, fc.min, fc.max)
        mu2 <- (1/fc) * mu0
        x[ii, jj] = rnbinom(n = 1, size = r, mu = mu2)
      }
    }
  }
  if(length(ind.down) > 0){
    for(jj in ind.down){
      for(ii in 1:n){
        fc = runif(1, fc.min, fc.max)
        mu1 <- (1/fc) * mu0
        x[ii, jj] = rnbinom(n = 1, size = r, mu = mu1)
      }
      for(ii in (n+1):(n+m)){
        fc = runif(1, fc.min, fc.max)
        mu2 <- fc * mu0
        x[ii, jj] = rnbinom(n = 1, size = r, mu = mu2)
      }
    }
  }
  
  x[,-c(ind.up, ind.down)] = matrix(rnbinom(n = (n+m)* (nn - nn1 - nn2), size = r, mu = mu0), nrow = (n+m))
  
  dat <- data.frame(y = c(rep(1, n), rep(0, m)), x, stringsAsFactors = F)
  colnames(dat)[1] <- 'y'
  colnames(dat)[2:ncol(dat)] <- ents.mRNA$id
  
  return(dat)
}

#### Randomizing the network: The network of the selected regulators (L$rels)
#### will be randomized as fllows. A fraction (n%) of the network is kept fixed and
#### the rest of the network is randomized with the relations in the bigger network
randomizeKB <- function(ents, rels, rel.uid.fix = NA)
{
  if(!is.na(rel.uid.fix[1])){
    fix.rels  <- rels[rels$uid %in% rel.uid.fix, ]
    rand.rels <- rels[!(rels$uid %in% rel.uid.fix), ]
    rand.rels$trguid = rand.rels$trguid[sample(sample(nrow(rand.rels)))]
    rand.rels$type   = rand.rels$type[sample(sample(nrow(rand.rels)))]
    rels <- rbind(fix.rels, rand.rels)
    rownames(rels) = 1:nrow(rels)
  }else{
    rels$trguid = rels$trguid[sample(sample(nrow(rels)))]
    rels$type   = rels$type[sample(sample(nrow(rels)))]
    rownames(rels) = 1:nrow(rels)
  }
  
  L = list(ents = ents, rels = rels)  
  L
}
