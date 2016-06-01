

winsorize <- function(x, n.win = 1){
###Winsorize extreme values
###Ex: Two most extreme values (from each side):
###t(apply(ed, 1, winsorize, 2))

    n.vals = length(x)
    fraction = n.win / n.vals
    if(length(fraction) != 1 || fraction < 0 || fraction > 0.5){
        stop("bad value for 'fraction'")
    }

    win.sorted.ind = order(x)
    x[win.sorted.ind[1:n.win]] = x[win.sorted.ind[n.win + 1]]
    x[win.sorted.ind[(n.vals - n.win + 1):n.vals]] = x[win.sorted.ind[n.vals - n.win]]

    ##fraction-based
    if(0){
        lim = quantile(x, probs = c(fraction, 1 - fraction))

        ##set extreme values
        x[ x < lim[1] ] = lim[1]
        x[ x > lim[2] ] = lim[2]
    }
    
    return(x)
}


pca.twodim.plot <- function(pca, pc.x, pc.y, meta.mat, strat.factor, obs.alpha = 0.9, xlim = NA, ylim = NA){

    suppressMessages(library('ggplot2'))

    ##get PC basis
    if(class(pca) == 'prcomp'){
        pca.basis = pca[['x']]
    }else{
        pca.basis = pca
    }
    
    ##subset
    pca.basis = pca.basis[, c(pc.x, pc.y)]
    pca.basis = data.frame(obsnames = row.names(pca.basis), pca.basis)
    meta.mat = meta.mat[rownames(pca.basis), ]
    
    ##add strat.factor for color-purposes
    pca.basis = cbind(pca.basis, meta.mat[, strat.factor], stringsAsFactors = FALSE)
    colnames(pca.basis)[ncol(pca.basis)] = strat.factor
    
    ##Plot as points
    g.plot = ggplot(pca.basis, aes_string(x = pc.x, y = pc.y, colour = strat.factor)) + geom_point(alpha = obs.alpha) ##color = obs.colors    

    ##axis lims
    if(is.numeric(xlim)){
        g.plot = g.plot + coord_cartesian(xlim = xlim)
    }
    if(is.numeric(ylim)){
        g.plot = g.plot + coord_cartesian(ylim = ylim)
    }
    if(is.numeric(ylim) & is.numeric(xlim)){
        g.plot = g.plot + coord_cartesian(ylim = ylim, xlim = xlim)
    }
    
    ##Turn off bg + grid
    g.plot = g.plot + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black"))

    ##legend
    if(!is.character(meta.mat[, strat.factor]) && !is.factor(meta.mat[, strat.factor])){
        ##g.plot = g.plot + scale_colour_gradient(low = '#C6DBEF', high = '#08306B', name = 'ICMvsTE\nmarkers')
        ##g.plot = g.plot + scale_colour_gradientn(colours = c('#C6DBEF', '#C6DBEF', '#08306B', '#08306B'), values = c(0, 0.3, 0.7, 1), name = 'ICMvsTE\nmarkers')
        g.plot = g.plot + scale_colour_gradient2()
        ##custom: scale_colour_gradient(palette = )
    }else{
        strat.factor.col = paste(strat.factor, '.color', sep = '')
        strat2col.map = unique(meta.mat[, c(strat.factor, strat.factor.col)])
        strat2col.map = strat2col.map[order(strat2col.map[, strat.factor]), ]
        if(is.factor(meta.mat[, strat.factor])){
            strat2col.map = strat2col.map[levels(meta.mat[, strat.factor]), ]
        }
        
        g.plot = g.plot + scale_colour_manual(breaks = strat2col.map[, strat.factor], values = strat2col.map[, strat.factor.col])
    }
    
    return(g.plot)
}

pca.biplot <- function(pca, pc.x = "PC1", pc.y = "PC2", plot.vars, meta.mat, obs.color.col = NA, obs.alpha = 0.9, var.size = 1.5, var.alpha = 0.85, var.colors = '#DE2D26') {
###pca is an object from prcomp

    suppressMessages(library('ggplot2'))
    
    ##Observations
    pca.basis = pca[['x']]
    
    ##subset
    pca.basis = pca.basis[, c(pc.x, pc.y)]
    pca.basis = data.frame(obsnames = row.names(pca.basis), pca.basis)

    ##add color
    pca.basis = cbind(pca.basis, meta.mat[rownames(pca.basis), obs.color.col], stringsAsFactors = FALSE)
    colnames(pca.basis)[ncol(pca.basis)] = obs.color.col
    
    ##Variables
    pca.rot = pca[['rotation']]
    
    ##subset
    pca.rot = pca.rot[plot.vars, c(pc.x, pc.y)]
    pca.rot = data.frame(varnames = rownames(pca.rot), pca.rot)

    ##rescale pca.rot such that on the same scale as pca.basis
    mult.factor = min( (max(pca.basis[, pc.y]) - min(pca.basis[, pc.y]) / (max(pca.rot[, pc.y]) - min(pca.rot[, pc.y]))), (max(pca.basis[, pc.x]) - min(pca.basis[, pc.x]) / (max(pca.rot[, pc.x]) - min(pca.rot[, pc.x]))) )
    pc.x.rescaled = pca.rot[, pc.x] * mult.factor * 1.2
    pc.y.rescaled = pca.rot[, pc.y] * mult.factor * 1.2
    pca.rot = cbind(pca.rot, pc.x.rescaled, pc.y.rescaled)

    ##Observations as points
    g.plot = ggplot(pca.basis, aes_string(x = pc.x, y = pc.y, colour = obs.color.col)) + geom_point(alpha = obs.alpha) ##color = obs.colors    
    
    ##Variables as text
    g.plot = g.plot + coord_equal() + geom_text(data = pca.rot, aes(x = pc.x.rescaled, y = pc.y.rescaled, label = varnames), alpha = var.alpha, size = var.size, color = var.colors, vjust = 0.5, hjust = 0.5)

    ##Turn off bg + grid
    g.plot = g.plot + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black"))

    ##Legend
    if(!is.character(meta.mat[, obs.color.col]) && !is.factor(meta.mat[, obs.color.col])){
        ##g.plot = g.plot + scale_colour_gradient(low = '#C6DBEF', high = '#08306B', name = 'ICMvsTE\nmarkers')
        ##g.plot = g.plot + scale_colour_gradientn(colours = c('#C6DBEF', '#C6DBEF', '#08306B', '#08306B'), values = c(0, 0.3, 0.7, 1), name = 'ICMvsTE\nmarkers')
        g.plot = g.plot + scale_colour_gradient2()
        ##custom: scale_colour_gradient(palette = )
    }else{

        strat.factor.color.col = paste(obs.color.col, '.color', sep = '')
        strat2col.map = unique(meta.mat[, c(obs.color.col, strat.factor.color.col)])
        strat2col.map = strat2col.map[order(strat2col.map[, obs.color.col]), ]
        if(is.factor(meta.mat[, obs.color.col])){
            strat2col.map = strat2col.map[levels(meta.mat[, obs.color.col]), ]
        }

        g.plot = g.plot + scale_colour_manual(breaks = strat2col.map[, obs.color.col], values = strat2col.map[, strat.factor.color.col])
    }
            
    ##Arrows
    ##g.plot = g.plot + geom_segment(data = pca.rot, aes(x=0, y=0, xend = pc.x.rescaled, yend = pc.y.rescaled), arrow = arrow(length=unit(0.2,"cm")), alpha = 0.75, color = var.colors)
    
    return(g.plot)
}

rrnaseq.combat <- function(data.mat, meta.mat, batch.factor, rm.single = TRUE){

    suppressMessages(library('sva'))

    n.rows = nrow(data.mat)
    n.cols = ncol(data.mat)
    new.n.rows = 0
    new.n.cols = 0
    
    while(n.rows != new.n.rows | n.cols != new.n.cols){
        n.rows = nrow(data.mat)
        n.cols = ncol(data.mat)
        
        ##rm rows (genes) with constant variance
        const.rows = which(apply(data.mat, 1, function(x){var(x) == 0}))
        if(length(const.rows) != 0){
            warning('There are constant rows. These are removed.')
            data.mat = data.mat[setdiff(1:nrow(data.mat), const.rows), ]
        }

        ##rm cols (samples) with constant variance
        const.cols = which(apply(data.mat, 2, function(x){var(x) == 0}))
        if(length(const.cols) != 0){
            warning('There are constant cols. These are removed.')
            data.mat = data.mat[, setdiff(1:ncol(data.mat), const.cols)]
        }
        meta.mat = meta.mat[colnames(data.mat), ]
        batch = meta.mat[, batch.factor]

        ##Filter out singletons
        if(rm.single){
            batch = meta.mat[, batch.factor]
            embryo2ncells = table(batch)
            embryo.singles = names(embryo2ncells)[which(embryo2ncells == 1)]
            embryo.keep = setdiff(meta.mat[, 'stage.embryo'], embryo.singles)        
            pass.meta = merge(meta.mat, as.matrix(embryo.keep), by.x = 'stage.embryo', by.y = 1)
            pass.samples = pass.meta[, 'id']
            data.mat = data.mat[, pass.samples]
            meta.mat = meta.mat[pass.samples, ]
            batch = meta.mat[, batch.factor]
        }

        new.n.rows = nrow(data.mat)
        new.n.cols = ncol(data.mat)
    }
    
    ##model.matrix
    modcombat = model.matrix(~1, data = meta.mat)
    
    ##combat
    data.mat.adj = ComBat(dat = data.mat, batch = batch, mod = modcombat, par.prior = TRUE, prior.plots = FALSE) #numCovs = NULL

    return(data.mat.adj)    
}

plot.pam.boxplots <- function(data.mat, row.clust, bpcol = 'white') {
  
  unique.clusters = sort(unique(row.clust))
  n.clusters = length(unique.clusters)
  
  for (jcluster in 1:n.clusters){
    cluster = unique.clusters[jcluster]
    
    rowsincluster = names(row.clust)[which(row.clust==cluster)]
    n.objects = length(rowsincluster)
    if(n.objects == 1){
      data2plot = as.data.frame(as.matrix(t(data.mat[rowsincluster, ])))
    }
    else{
      data2plot = as.data.frame(data.mat[rowsincluster, ])
    }
    boxplot(data2plot, main = sprintf("Cluster %i, N=%i", cluster, n.objects), col=bpcol, pch=".", las=2)
    abline(h=0, lty = 'dashed', col = 'grey')
  }
}

filter.col <- function(sig.res, filt.cutoff, filt.col = 'log2.fc', abs.bool = TRUE){

    if(abs.bool){
        filt.res = sig.res[which(abs(sig.res[, filt.col]) >= filt.cutoff), ]
    }else{
        filt.res = sig.res[which(sig.res[, filt.col] >= filt.cutoff), ]
    }

    return(filt.res)
}

bigcor.par <- function(MAT, nblocks = 10, verbose = TRUE, ncores=20, yMAT = NA, ...){
#https://gist.github.com/bobthecat/5024079
    
  suppressMessages(library(ff, quietly = TRUE))
  suppressMessages(library(doMC))
  
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
    color.legend = FALSE
    if(!is.logical(meta.mat)){
        meta.mat = meta.mat[samples, ]        
        cell.color.col = paste(factor.color, 'color', sep = '.')
        cell.color = meta.mat[, cell.color.col]
        
        if(is.character(meta.mat[, factor.color]) || is.factor(meta.mat[, factor.color])){ #else numeric values on contiounous
            
            color.legend = TRUE            
            n.mixed = length(grep('; color: ', meta.mat[, factor.color]))
            if(n.mixed != 0){
                cell2color.map = t(as.data.frame(sapply(meta.mat[1:n.mixed, factor.color], strsplit, '; color: ')))
                colnames(cell2color.map) = c(factor.color, cell.color.col)
            }else{            
                cell2color.map = unique(meta.mat[, c(factor.color, cell.color.col)])
            }

            ##order and labels
            if(is.factor(meta.mat[, factor.color])){
                color.levels = levels(meta.mat[, factor.color])
                cell2color.map = cell2color.map[color.levels, ]
                cell2color.map.labels = color.levels
            }else{
                cell2color.map = cell2color.map[order(cell2color.map[, factor.color]), ]
                cell2color.map.labels = cell2color.map[, factor.color]
            }
            
            cell2color.map.col = cell2color.map[, cell.color.col]
        }
    }else{
        cell.color = par('col')
    }
        
    ##set up panels
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

        txt.y = 0.8
        if(color.legend == TRUE){
            txt.y = 0.3
            legend('topleft', legend = cell2color.map.labels, col = cell2color.map.col, lty = 1, cex = 0.7)
        }
        text(0.5, txt.y, sprintf('%s: %.1f%s', labels, e.var[grep(paste('^', labels, '$', sep = ''), colnames(pca.basis))] * 100, '%'));            
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
