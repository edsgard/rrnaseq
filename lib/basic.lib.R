
bigcor.par <- function(MAT, nblocks = 10, verbose = TRUE, ncores=20, yMAT = NA, ...){
#https://gist.github.com/bobthecat/5024079
    
  library(ff, quietly = TRUE)
  library(doMC)
  
  if(ncores=="all"){
      ncores = multicore:::detectCores()
      registerDoMC(cores = ncores)
  } else{
      registerDoMC(cores = ncores)
  }
  

  #add dummy cols such that modulus nblocks == 0
  NCOL.predummy = ncol(MAT)
  if (NCOL.predummy %% nblocks != 0){
      n.dummy = nblocks - (NCOL.predummy %% nblocks)
      dummy.mat = matrix(NA, nrow = nrow(MAT), ncol = n.dummy)
      MAT = cbind(MAT, dummy.mat)
  }
  
  NCOL = ncol(MAT)

  if(length(yMAT) == 1){
      if(is.na(yMAT)){
          yMAT = MAT
          singlemat = TRUE
      }
  }else{
      singlemat = FALSE
      y.ncol.predummy = ncol(yMAT)
      if(y.ncol.predummy != NCOL.predummy){
          stop('number of columns of MAT is not equal to that of yMAT')
      }else{
          if(NCOL.predummy %% nblocks != 0){
              yMAT = cbind(yMAT, dummy.mat)
          }
      }
  }
  
  ## preallocate square matrix of dimension
  ## ncol(MAT) in 'ff' single format
  corMAT = ff(vmode = "single", dim = c(NCOL, NCOL))
 
  ## split column numbers into 'nblocks' groups.
  #NB: requires that NCOL %% nblocks == 0
  SPLIT = split(1:NCOL, rep(1:nblocks, each = NCOL / nblocks))
 
	## create all unique combinations of blocks
	COMBS = expand.grid(1:length(SPLIT), 1:length(SPLIT))
	COMBS = t(apply(COMBS, 1, sort))
	COMBS = unique(COMBS)
 
	## iterate through each block combination, calculate correlation matrix
	## between blocks and store them in the preallocated matrix on both
	## symmetric sides of the diagonal
	results = foreach(i = 1:nrow(COMBS)) %dopar% {
		COMB = COMBS[i, ]
		G1 = SPLIT[[COMB[1]]]
		G2 = SPLIT[[COMB[2]]]
		if (verbose) cat("Block", COMB[1], "with Block", COMB[2], "\n")
		flush.console()
		COR = cor(MAT[, G1], yMAT[, G2], ...)
		corMAT[G1, G2] = COR
                if(singlemat){ #cor-mat is then symmetric
                    corMAT[G2, G1] = t(COR)
                }else{                
                    COR = cor(MAT[, G2], yMAT[, G1], ...)
                    corMAT[G2, G1] = COR
                }
		COR = NULL
	}
 
	gc()

  #remove dummy cols
  corMAT = corMAT[1:NCOL.predummy, 1:NCOL.predummy]
  MAT = MAT[, 1:NCOL.predummy]
  
  #set row- and colnames
  colnames(corMAT) = colnames(MAT)
  rownames(corMAT) = colnames(MAT)
  
  return(corMAT)
}

rm.const.vec <- function(data.mat, col.rm = TRUE, row.rm = TRUE){

    #rm cols with constant variance
    if(col.rm){
        const.cols = which(apply(data.mat, 2, function(x){var(x) == 0}))
        if(length(const.cols) != 0){
            warning('There were constant columns. These were removed.')
            data.mat = data.mat[, setdiff(1:ncol(data.mat), const.cols)]
        }
    }

    #rm rows (genes) with constant variance
    if(row.rm){
    
        const.rows = which(apply(data.mat, 1, function(x){var(x) == 0}))
        if(length(const.rows) != 0){
            warning('There were constant rows. These were removed.')
            data.mat = data.mat[setdiff(1:nrow(data.mat), const.rows), ]
        }
    }

    return(data.mat)
}

spear.cor.col.dist <- function(x, ...){
#For use by pvclust
    
    cor.res = cor(x, use="pairwise.complete.obs", method = 'spearman')
    
    #get dist
    dist.mat = as.dist(((cor.res * -1) + 1) / 2)

    attr(dist.mat, "method") = 'spearman'
    
    return(dist.mat)
}

spear.cor.dist <- function(data.mat, ...){
#For use by heatmap.2
#This function computes and returns the distance matrix computed by
#     using the specified distance measure (here: correlation method) to compute the distances
#     between the rows of a data matrix.
    
    cor.res = cor(t(data.mat), use="pairwise.complete.obs", method = 'spearman')
    
    #get dist
    dist.mat = as.dist(((cor.res * -1) + 1) / 2)
    
    return(dist.mat)
}

pears.cor.dist <- function(data.mat, ...){
#This function computes and returns the distance matrix computed by
#     using the specified distance measure (here: correlation method) to compute the distances
#     between the rows of a data matrix.
    
    cor.res = cor(t(data.mat), use="pairwise.complete.obs", method = 'pearson')
    
    #get dist
    dist.mat = as.dist(((cor.res * -1) + 1) / 2)
    
    return(dist.mat)
}

plot.pairs.pca <- function(pca.basis, e.var, meta.mat = NA, factor.color = 'stage', point.text.cex = 0.7, points.cex = 1){

    #set up color mappings
    samples = rownames(pca.basis)
    if(!is.logical(meta.mat)){
        cell.color.col = paste(factor.color, 'color', sep = '.')
        cell.color = meta.mat[samples, cell.color.col]
        cell2color.map = unique(meta.mat[samples, c(factor.color, cell.color.col)])
        cell2color.map.col = cell2color.map[, cell.color.col]
        cell2color.map.labels = cell2color.map[, factor.color]
    }else{
        cell.color = par('col')
    }
        
    #set up panels
    point.text.panel <- function(x, y){text(x, y, labels = samples, col = cell.color, cex = point.text.cex)}
    point.panel <- function(x, y){points(x, y, col = cell.color, pch = 16, cex = points.cex)}
    diag.panel <- function(x, y, labels, ...){
        #density histogram
        pu <- par("usr")
        d <- density(x,...)
        par("usr" = c(pu[1:2], 0, max(d$y)*1.5))
        lines(d)
    }
    text.panel <- function(x, y, labels, ...){

        if(labels == 'PC1'){
            txt.y = 0.8
            if(!is.logical(meta.mat)){
                txt.y = 0.3
                legend('topleft', legend = cell2color.map.labels, col = cell2color.map.col, lty = 1, cex = 0.7)
            }
            text(0.5, txt.y, sprintf('%s: %.1f%s', labels, e.var[grep(paste('^', labels, '$', sep = ''), colnames(pca.basis))] * 100, '%'));
        }
        else{
            text(0.5, 0.8, sprintf('%s: %.1f%s', labels, e.var[grep(paste('^', labels, '$', sep = ''), colnames(pca.basis))] * 100, '%'));
        }
    }
                                                 
    #plot
    pairs(pca.basis, upper.panel = point.text.panel, lower.panel = point.panel, diag.panel = diag.panel, text.panel = text.panel)
    
}

plot.hclust.my <- function(col.clust, meta.mat, strat.factor, cex){
     
    d = dendrapply(as.dendrogram(col.clust), labelCol, cex = cex, meta.mat = meta.mat, strat.factor = strat.factor)
    plot(d)

    #legend
    strat.factor.col = paste(strat.factor, 'color', sep = '.')
    strat.factor2color.map = unique(meta.mat[, c(strat.factor, strat.factor.col)])
    legend('topright', legend = strat.factor2color.map[, strat.factor], col = strat.factor2color.map[, strat.factor.col], lty = 1)
    
}

labelCol <- function(x, cex = 0.3, meta.mat, strat.factor = 'stage') {

    strat.factor.col = paste(strat.factor, 'color', sep = '.')
    if (is.leaf(x)) {
                
        ## fetch label
        label = attr(x, "label")
            
        ## set label color
        factor.color = as.character(meta.mat[label, strat.factor.col])
        attr(x, "nodePar") <- list(lab.col = factor.color, lab.cex = cex, cex = cex)
    }
    return(x)
}
