

rseq.hclust <- function(rpkm, cor.meth, cor.res.list, store.col = 'col', ncores = 1, nblocks = 1){

    cor.store.col = paste(store.col, 'cor', sep = '.')
    hclust.store.col = paste(store.col, 'hclust', sep = '.')
    
    if(is.na(cor.res.list)){

        #rm cols with constant variance
        rpkm = rm.const.vec(rpkm, row.rm = FALSE)
        
        #Pairwise dist btw columns
        cor.res = bigcor.par(rpkm, use="pairwise.complete.obs", method = cor.meth, ncores = ncores, nblocks = nblocks)
    
        #get dist and hclust 
        col.dist.mat = as.dist(((cor.res * -1) + 1) / 2)
        clust.res = hclust(col.dist.mat)

        cor.res.list = list()
        cor.res.list[[cor.store.col]] = cor.res
        cor.res.list[[hclust.store.col]] = clust.res
    }else{

        if(length(which(names(cor.res.list) == cor.store.col)) != 0){
            cor.res = cor.res.list[[cor.store.col]]

            #get clust.res
            if(length(which(names(cor.res.list) == hclust.store.col)) != 0){
                clust.res = cor.res.list[[hclust.store.col]]
            }else{            
                col.dist.mat = as.dist(((cor.res * -1) + 1) / 2)
                clust.res = hclust(col.dist.mat)

                #store result
                cor.res.list[[hclust.store.col]] = clust.res
            }
        }else{

            #rm cols with constant variance
            rpkm = rm.const.vec(rpkm, row.rm = FALSE)
        
            #Pairwise dist btw columns
            cor.res = bigcor.par(rpkm, use="pairwise.complete.obs", method = cor.meth, ncores = ncores, nblocks = nblocks)

            #get dist and hclust 
            col.dist.mat = as.dist(((cor.res * -1) + 1) / 2)
            clust.res = hclust(col.dist.mat)

            #store result
            cor.res.list[[cor.store.col]] = cor.res
            cor.res.list[[hclust.store.col]] = clust.res
        }                
    }
    
    return(cor.res.list)
}


add.qc.col <- function(qc.meta.mat, qc.col, def.val = 0){
    if(length(which(colnames(qc.meta.mat) == qc.col)) == 0){
        new.qc.col = rep(def.val, nrow(qc.meta.mat))
        qc.meta.mat = cbind(qc.meta.mat, new.qc.col)
        colnames(qc.meta.mat)[ncol(qc.meta.mat)] = qc.col
    }
    return(qc.meta.mat)
}

get.qc.df <- function(qc.meta.file, meta.mat){
    
    if(!file.exists(qc.meta.file)){
        samples = meta.mat[, 'id']
        n.samples = length(samples)
        #qc.cols = c('n.reads', 'prc.reads.mapped', 'n.genes.expr', 'spear.cor') #pca
        qc.meta.mat = as.data.frame(matrix(0, nrow = n.samples, ncol = 0))
        rownames(qc.meta.mat) = samples
    }else{
        qc.meta.mat = readRDS(qc.meta.file)
    }

    return(qc.meta.mat)
}

plot.heatmap.sample <- function(cor.res, col.clust, meta.mat, strat.factor, ...){
    
    library('gplots')
    strat.factor.col = paste(strat.factor, '.color', sep = '')
    stratum.col = as.character(meta.mat[colnames(cor.res), strat.factor.col])
    strat.factor2color.map = unique(meta.mat[, c(strat.factor, strat.factor.col)])
    
    #colormaps
    if(0){
        col.pa = greenred(100) #alt: maPalette{marray)
        #col.pa = colorRampPalette(brewer.pal(9, "Blues")[4:9])(100)
    }
    
    ###
    #Plot
    ###
    heatmap.2(cor.res, Rowv = as.dendrogram(col.clust), Colv = as.dendrogram(col.clust), trace = 'none', symm = TRUE, ColSideColors = stratum.col, ...)
    legend('topright', legend = strat.factor2color.map[, strat.factor], col = strat.factor2color.map[, strat.factor.col], lty = 1)

}

gene2nsamples.expr.filter <- function(rpkm, ncells.cutoff, rpkm.cutoff, ncell.strat.factor, meta.mat, n.strata.cutoff = 1, nstrata.strat.factor = 'nostrat', nstrata.ncells.cutoff = 1, nstrata.nstrata.cutoff = 1){
        
    #number of samples (cells) with expression per stratum
    gene2nsamples.expr.list = gene2nsamples(rpkm, meta.mat, ncell.strat.factor, rpkm.cutoff)

    #ncells.cutoff per ncell.strat.factor
    fail.genes.list = lapply(gene2nsamples.expr.list, function(jstratum.gene2nsamples, ncells.cutoff){names(jstratum.gene2nsamples)[which(jstratum.gene2nsamples < ncells.cutoff)]}, ncells.cutoff = ncells.cutoff)

    #get genes that passed the ncells.cutoff in n.strata
    n.strata = length(gene2nsamples.expr.list)
    fail.genes2nstrata = table(unlist(fail.genes.list))
    pass.genes2nstrata = n.strata - fail.genes2nstrata
    fail.genes = names(pass.genes2nstrata)[which(pass.genes2nstrata < n.strata.cutoff)]
    pass.genes = setdiff(rownames(rpkm), fail.genes)

    print(sprintf('%i genes expressed in >=%i cells within >=%i strata, where stratification is based on %s, and expr cutoff %.1f', length(pass.genes), ncells.cutoff, n.strata.cutoff, ncell.strat.factor, rpkm.cutoff))
    
    #number of strata where expr is seen in at least nstrata.ncells.cutoff
    gene2nsamples.expr.list = gene2nsamples(rpkm, meta.mat, nstrata.strat.factor, rpkm.cutoff)
    gene2nstrata = table(unlist(lapply(gene2nsamples.expr.list, function(jstratum.gene2nsamples, nstrata.ncells.cutoff){names(jstratum.gene2nsamples)[which(jstratum.gene2nsamples >= nstrata.ncells.cutoff)]}, nstrata.ncells.cutoff = nstrata.ncells.cutoff)))
    nstrata.pass.genes = names(gene2nstrata)[which(gene2nstrata >= nstrata.nstrata.cutoff)]

    pass.genes = intersect(pass.genes, nstrata.pass.genes)
    print(sprintf('%i genes expressed as above and also expressed in >=%i cells within >=%i strata, stratified on %s', length(pass.genes), nstrata.ncells.cutoff, nstrata.nstrata.cutoff, nstrata.strat.factor))

    #subset data
    rpkm = rpkm[pass.genes, ]

    return(rpkm)
}

plot.gene2nsamples.expr <- function(rpkm, meta.mat, strat.factor, rpkm.cutoff, n.cells.cutoff){
    
    library('RColorBrewer')

    #number of cells with expression per stratum
    gene2nsamples.expr.list = gene2nsamples(rpkm, meta.mat, strat.factor, rpkm.cutoff)
    
    #"hist"
    tables = lapply(gene2nsamples.expr.list, table)
    max.cells = max(as.integer(unlist(lapply(tables, names))))
    breaks = seq(0, max.cells, by = 1)
    hist.list = lapply(gene2nsamples.expr.list, hist, breaks = breaks, plot = FALSE)
  
    #colormap
    strat.factor.col = paste(strat.factor, '.color', sep = '')
    stratum2color.map = unique(meta.mat[, c(strat.factor, strat.factor.col)])
    rownames(stratum2color.map) = stratum2color.map[, strat.factor]

    #axis limits
    xlim = c(0, max.cells)
    ylim = range(unlist(lapply(hist.list, '[[', 'density')))
    ylim[2] = ylim[2] + ylim[2] * 0.1
  
    #plot
    strata = names(hist.list)
    j.stratum = 1
    stratum = strata[j.stratum]
    j.hist = hist.list[[stratum]]
    plot(x = j.hist[['mids']] - 0.5, j.hist[['density']], col = stratum2color.map[stratum, strat.factor.col], xlab = '# of cells', xlim = xlim, main = '', type = 'p', pch = 16, ylim = ylim, ylab = 'Fraction of all genes')
    for(j.stratum in 1:length(strata)){
      stratum = strata[j.stratum]
      j.hist = hist.list[[stratum]]
      points(x = j.hist[['mids']] - 0.5, j.hist[['density']], col = stratum2color.map[stratum, strat.factor.col], xlim = xlim, pch = 16)
    }
    legend('topright', legend = stratum2color.map[, strat.factor], col = stratum2color.map[, strat.factor.col], lty = 1)
    abline(v = n.cells.cutoff - 0.1)

}
    
gene2nsamples <- function(rpkm, meta.mat, strat.factor, rpkm.cutoff){

    #Ensure that the same set of samples in rpkm and meta.mat
    shared.samples = intersect(colnames(rpkm), rownames(meta.mat))
    rpkm = rpkm[, shared.samples]
    meta.mat = meta.mat[shared.samples, ]
    
    #number of cells with expression per stratum
    if(strat.factor == 'unstrat'){
        stratum2id = list(unstrat = colnames(rpkm))
    }else{
        stratum2id = tapply(meta.mat[, 'id'], meta.mat[, strat.factor], unique)
    }
    gene2nsamples.expr.list = lapply(stratum2id, function(cells, rpkm, rpkm.cutoff){rpkm = rpkm[, cells, drop = FALSE]; gene2nsamples.expr = apply(rpkm, 1, function(jgene.rpkm, rpkm.cutoff){length(which(jgene.rpkm > rpkm.cutoff));}, rpkm.cutoff = rpkm.cutoff)}, rpkm = rpkm, rpkm.cutoff = rpkm.cutoff)

    return(gene2nsamples.expr.list)
}

plot.strata.expr.dist <- function(expr.mat, meta.mat, strat.factor){

    #rm columns with < 2 observations != inf (logged 0's)
    few.datapoints.samples = colnames(expr.mat)[which(apply(expr.mat, 2, function(x){x = setdiff(x, c(Inf, -Inf)); length(unique(x))}) < 2)]
    pass.samples = setdiff(colnames(expr.mat), few.datapoints.samples)
    expr.mat = expr.mat[, pass.samples]
    if(length(few.datapoints.samples) >0){
        warning(sprintf('There were %i columns with less than two datapoints (after Inf, -Inf removal). These columns were excluded.', length(few.datapoints.samples)))
    }

    meta.mat = meta.mat[pass.samples, ]
    
    #Get dists
    meta.list = split(meta.mat, meta.mat[, strat.factor])
    strata = names(meta.list)
    density.list = list()
    for(j.strat in strata){
        j.meta = meta.list[[j.strat]]

        #subset
        j.samples = j.meta[, 'id']
        j.rpkm = expr.mat[, j.samples]

        density.list[[j.strat]] = density(j.rpkm)
    }
        
    #colormap
    strat.factor.col = paste(strat.factor, '.color', sep = '')
    strat2color.map = unique(meta.mat[, c(strat.factor, strat.factor.col)])
    rownames(strat2color.map) = strat2color.map[, strat.factor]
    
    #get ranges
    xlim = range(unlist(lapply(density.list, '[[', 'x'))) #xlim = c(-100, 100)
    ylim = range(unlist(lapply(density.list, '[[', 'y')))
    xlim = c(floor(xlim[1]), ceiling(xlim[2]))

    #plot
    n.strata = length(strata)
    j.strat.it = 1
    j.strat = strata[j.strat.it]
    plot(density.list[[j.strat]], col = strat2color.map[j.strat, strat.factor.col], xlim = xlim, ylim = ylim, xaxt = 'n', xlab = 'log10(RPKM)', ylab = 'Density: n.genes', main = '')
    for(j.strat.it in 2:n.strata){
        j.strat = strata[j.strat.it]
        lines(density.list[[j.strat]], col = strat2color.map[j.strat, strat.factor.col])
    }
    x.vec = seq(xlim[1], xlim[2], by = 1)
    axis(1, at = x.vec, labels = x.vec, las = 2, cex.axis = 0.5)
    legend('topright', legend = strat2color.map[, strat.factor], col = strat2color.map[, strat.factor.col], lty = 1)

}

plot.sample.expr.dist <- function(expr.mat, meta.mat, strat.factor){

    #rm columns with < 2 observations != inf (logged 0's)
    few.datapoints.samples = colnames(expr.mat)[which(apply(expr.mat, 2, function(x){x = setdiff(x, c(Inf, -Inf)); length(unique(x))}) < 2)]
    expr.mat = expr.mat[, setdiff(colnames(expr.mat), few.datapoints.samples)]
    if(length(few.datapoints.samples) >0){
        warning(sprintf('There were %i columns with less than two datapoints (after Inf, -Inf removal). These columns were excluded.', length(few.datapoints.samples)))
    }
    
    #Get dists
    density.list = apply(expr.mat, 2, density)
    
    #colormap
    strat.factor.col = paste(strat.factor, '.color', sep = '')
    strat2color.map = unique(meta.mat[, c(strat.factor, strat.factor.col)])

    #get ranges
    xlim = range(unlist(lapply(density.list, '[[', 'x'))) #xlim = c(-100, 100)
    ylim = range(unlist(lapply(density.list, '[[', 'y')))
    xlim = c(floor(xlim[1]), ceiling(xlim[2]))

    #plot
    n.samples = length(density.list)
    j.sample = 1
    plot(density.list[[j.sample]], xlim = xlim, ylim = ylim, col = meta.mat[j.sample, strat.factor.col], main = '', xlab = 'log10(RPKM)', xaxt = 'n')
    for(j.sample in 2:n.samples){
        lines(density.list[[j.sample]], col = meta.mat[j.sample, strat.factor.col])
    }
    x.vec = seq(xlim[1], xlim[2], by = 1)
    axis(1, at = x.vec, labels = x.vec, las = 2, cex.axis = 0.5)    
    legend('topright', legend = strat2color.map[, strat.factor], col = strat2color.map[, strat.factor.col], lty = 1)

}

get.rpkm.files <- function(rpkm.cloud.dir, rpkm.raw.dir){
    rpkm.file = file.path(rpkm.cloud.dir, 'rpkm.rds')
    counts.file = file.path(rpkm.raw.dir, 'counts.rds')
    ercc.counts.file = file.path(rpkm.cloud.dir, 'ercc.counts.rds')
    ercc.rpkm.file = file.path(rpkm.cloud.dir, 'ercc.rpkm.rds')
    rds.files = c(rpkm.file, counts.file, ercc.counts.file, ercc.rpkm.file)
    names(rds.files) = gsub('\\.rds$', '', unlist(lapply(rds.files, basename)))

    return(rds.files)
}

get.tab.files <- function(rpkm.raw.dir){
    rpkm.tab.file = file.path(rpkm.raw.dir, 'rpkm.tab')
    counts.tab.file = file.path(rpkm.raw.dir, 'counts.tab')
    tab.files = c(rpkm.tab.file, counts.tab.file)
    names(tab.files) = gsub('\\.tab$', '', unlist(lapply(tab.files, basename)))

    return(tab.files)
}

read.rpkm <- function(rpkm.files){
    
    #the first comment line #samples should be kept, ignore the other comment lines
    base.files = basename(rpkm.files)
    n.files = length(base.files)
    n.uniq.files = length(unique(base.files))
    if(n.files != n.uniq.files){
    	 base.files = paste(base.files, 1:n.files, sep = '.')
    }	 
    names(base.files) = rpkm.files

    counts.list = list()
    length(counts.list) = n.files
    names(counts.list) = base.files
    
    rpkm.list = list()
    length(rpkm.list) = n.files
    names(rpkm.list) = base.files
    for(j.file in rpkm.files){
        print(j.file)
        
        header = read.table(j.file, sep = '\t', stringsAsFactors = FALSE, nrows = 1, comment.char = '')
        header = header[2:length(header)]
        header = gsub('_unique.bam', '', header)
        data.df = read.table(j.file, sep = '\t', header = FALSE, stringsAsFactors = FALSE)
        
        gene.id = data.df[, 1]
        nm.id = data.df[, 2]
        gene2nm.mat = cbind(gene.id, nm.id)
        data.df = data.df[, setdiff(1:ncol(data.df), 1:2)]
	n.data.samples = ncol(data.df)
	n.samples = n.data.samples / 2
	header = header[1:n.samples]
        rpkm = as.matrix(data.df[, 1:n.samples])
        counts = as.matrix(data.df[, (n.samples+1):ncol(data.df)])

        #Use second ID column for the ERCC ids
        ercc.ind = grep('ERCC-', nm.id)
        gene.id[ercc.ind] = nm.id[ercc.ind]

        #Set gene.id as rowname
        colnames(rpkm) = header
        colnames(counts) = header
        rownames(rpkm) = gene.id
        rownames(counts) = gene.id
        rpkm.list[[base.files[j.file]]] = rpkm
        counts.list[[base.files[j.file]]] = counts
    }
    
    #handle non-unique gene-names (change to uniq.genes names)
    rpkm.list = lapply(rpkm.list, set.unique.rownames)
    counts.list = lapply(counts.list, set.unique.rownames)
    
    #bind rpkm
    n.files = length(rpkm.list)
    gene2freq = table(unlist(lapply(rpkm.list, rownames)))
    shared.genes = names(gene2freq)[which(gene2freq == n.files)]
    setdiff.genes = setdiff(names(gene2freq), shared.genes)
    n.diff.genes = length(setdiff.genes)
    if(n.diff.genes != 0){
      warning(sprintf('The gene annotation differs between your files! There were %i differing genes, among these are: %s', n.diff.genes, paste(head(setdiff.genes, 30), collapse = ', ')))
    }
    n.genes = length(shared.genes)
    all.rpkm = matrix(nrow = n.genes, ncol = 0)
    for(j.rpkm in rpkm.list){
        all.rpkm = cbind(all.rpkm, j.rpkm[shared.genes, ])
    }
    rpkm = all.rpkm
    
    #Exclude ERCC ids
    ercc.id = rownames(rpkm)[grep('^ERCC-', rownames(rpkm))]
    rpkm = rpkm[setdiff(rownames(rpkm), ercc.id), ]

    #bind counts
    ##shared.genes = rownames(counts.list[[1]])
    ##n.genes = length(shared.genes)
    all.counts = matrix(nrow = n.genes, ncol = 0)
    for(j.counts in counts.list){
        all.counts = cbind(all.counts, j.counts[shared.genes, ])
    }
    counts = all.counts
    
    #Exclude ERCC ids
    ercc.id = rownames(counts)[grep('^ERCC-', rownames(counts))]
    counts = counts[setdiff(rownames(counts), ercc.id), ]
    
    
    
    #####
    #Get ERCC
    #####
    ercc.id = unique(unlist(lapply(rpkm.list, function(jsample.rpkm){rownames(jsample.rpkm)[grep('^ERCC-', rownames(jsample.rpkm))]})))    
    #ercc.id = rownames(rpkm.list[[1]])[grep('^ERCC-', rownames(rpkm.list[[1]]))]
    n.ercc = length(ercc.id)
    if(n.ercc != 92){
        warning('Number of ERCC ids identified differs from 92!')
    }
    
    #bind rpkm
    all.rpkm = matrix(nrow = n.ercc, ncol = 0)
    for(j.rpkm in rpkm.list){
        j.ercc.id = intersect(rownames(j.rpkm), ercc.id)
        if(length(j.ercc.id) != n.ercc){
            warning("All ERCC not present in rpkm file.")
            j.ercc = matrix(NA, nrow = n.ercc, ncol = ncol(j.rpkm), dimnames = list(ercc.id, colnames(j.rpkm)))
        }else{
            j.ercc = j.rpkm[ercc.id, ]
        }
        all.rpkm = cbind(all.rpkm, j.ercc)
    }
    rownames(all.rpkm) = ercc.id
    ercc.rpkm = all.rpkm
    
    #bind counts
    all.counts = matrix(nrow = n.ercc, ncol = 0)
    for(j.counts in counts.list){
        j.ercc.id = intersect(rownames(j.counts), ercc.id)
        if(length(j.ercc.id) != n.ercc){
            warning("All ERCC not present in counts file.")
            j.ercc = matrix(NA, nrow = n.ercc, ncol = ncol(j.counts), dimnames = list(ercc.id, colnames(j.counts)))
        }else{
            j.ercc = j.counts[ercc.id, ]
        }
        all.counts = cbind(all.counts, j.ercc)        
    }
    rownames(all.counts) = ercc.id
    ercc.counts = all.counts
    
    
    return(list(rpkm = rpkm, counts = counts, ercc.rpkm = ercc.rpkm, ercc.counts = ercc.counts))
}

read.haplo.rpkm <- function(rpkm.files, id = 'gene'){

    #the first comment line #samples should be kept, ignore the other comment lines
    base.files = basename(rpkm.files)
    n.files = length(base.files)
    
    rpkm.list = list()
    length(rpkm.list) = n.files
    names(rpkm.list) = base.files
    for(j.file in rpkm.files){
        print(j.file)
        
        header = read.table(j.file, sep = '\t', stringsAsFactors = FALSE, nrows = 1, comment.char = '')
        header = header[2:length(header)]
        data.df = read.table(j.file, sep = '\t', header = FALSE, stringsAsFactors = FALSE)

        if(id == 'gene'){
            gene.id = data.df[, 1]
        }
        if(id == 'nm'){
            gene.id = data.df[, 2]
        }
        data.df = data.df[, setdiff(1:ncol(data.df), 1:2)]
        n.samples = length(header)
        rpkm = as.matrix(data.df[, 1:n.samples])
        
        colnames(rpkm) = header
        rownames(rpkm) = gene.id
        rpkm.list[[basename(j.file)]] = rpkm
    }

    #rm ERCC
    ercc.id = gene.id[grep('^ERCC_', gene.id)]
    rpkm.list = lapply(rpkm.list, function(rpkm, ercc.id){rpkm[setdiff(rownames(rpkm), ercc.id), ];}, ercc.id = ercc.id)
    
    #handle non-unique gene-names (change to uniq.genes names)
    rpkm.list = lapply(rpkm.list, set.unique.rownames)
    
    return(rpkm.list)
}

set.unique.rownames <- function(rpkm){
       
    gene.names = rownames(rpkm)

    #add uid and gene.uid col
    n.genes = length(gene.names)
    uid = 1:n.genes
    gene2uid = cbind(gene.names, uid)
    gene.uid = apply(gene2uid, 1, function(jrow){paste(jrow[1], '.uid_', jrow[2], sep = '')})
    gene2uid = cbind(gene2uid, gene.uid)

    #set rpkm rownames to gene.uid (no reordering has taken place)
    rownames(rpkm) = gene2uid[, 'gene.uid']

    #strip uid suffix from unique genes    
    gene2n = sort(table(gene.names), decreasing = TRUE)
    gene2uid = merge(as.data.frame(gene2uid, stringsAsFactors = FALSE), as.matrix(gene2n), by.x = 'gene.names', by.y = 'row.names')
    colnames(gene2uid)[ncol(gene2uid)] = 'n.rep'

    strip.ind = which(gene2uid[, 'n.rep'] == 1)
    gene.uid.stripped = gene2uid[, 'gene.uid']
    names(gene.uid.stripped) = gene2uid[, 'gene.uid']
    gene.uid.stripped[strip.ind] = gsub('\\.uid_[0-9]+$', '', gene.uid.stripped[strip.ind])

    #set rpkm rownames to gene.uid.stripped
    rownames(rpkm) = gene.uid.stripped[rownames(rpkm)]
            
    return(rpkm)
}

meta.add.colormaps <- function(meta.mat){

    library('RColorBrewer')
    
    #stage2color
    stages = sort(unique(meta.mat[, 'stage']))
    n.stages = length(stages)
    stage2color.map = 2:(n.stages + 1)
    names(stage2color.map) = stages
    stage.color = stage2color.map[meta.mat[, 'stage']]
    meta.mat = cbind(meta.mat, stage.color)

    #stage.ind2color
    stage.ind = gsub(' ', '', apply(meta.mat[, c('stage', 'individual')], 1, paste, collapse = '.'))
    meta.mat = cbind(meta.mat, stage.ind)
    factor = 'stage.ind'
    factors = sort(unique(meta.mat[, factor]))
    n.factors = length(factors)
    factor2color.map = brewer.pal(n.factors, 'Set3') #Possibly better, because Set3 is so bright: pal8 = brewer.pal(8, 'Dark2'); colorRampPalette(pal8)(10)
    names(factor2color.map) = factors
    factor.color = factor2color.map[meta.mat[, factor]]    
    meta.mat = cbind(meta.mat, factor.color)
    colnames(meta.mat)[ncol(meta.mat)] = paste(factor, 'color', sep = '.')
    
    #Sort by id
    meta.mat = meta.mat[order(meta.mat[, 'id']), ]

    #coerce all factor cols to chars
    factor.cols = which(unlist(lapply(meta.mat, is.factor)))
    meta.mat[, factor.cols] = lapply(meta.mat[, factor.cols], as.character)

    return(meta.mat)
}

mixcolor.mat <- function(colors.mat, mix.factors, mix.factors.pal){
###mix colors (additive mixing)
    library('colorspace')
    library('RColorBrewer')

    n.factors = length(mix.factors)

    ##get colors
    factor1 = mix.factors[1]
    factor1.color.col = paste(factor1, 'color', sep = '.')
    factor1.color = colors.mat[, factor1.color.col]

    factor2 = mix.factors[2]
    factor2.color.col = paste(factor2, 'color', sep = '.')
    factor2.color = colors.mat[, factor2.color.col]

    ##mix
    mixed.color = hex(mixcolor(0.5, hex2RGB(factor1.color), hex2RGB(factor2.color), where = 'RGB'))

    ##if additional factors, mix them NON-commutatively
    if(n.factors >= 3){
        for(j.factor in 3:n.factors){

            ##get colors
            factor1.color = mixed.color
            
            factor2 = mix.factors[j.factor]
            factor2.color.col = paste(factor2, 'color', sep = '.')
            factor2.color = colors.mat[, factor2.color.col]

            ##mix
            mixed.color = hex(mixcolor(0.5, hex2RGB(factor1.color), hex2RGB(factor2.color), where = 'RGB'))
        }    
    }
    
    ##color legend
    mixed.factor.val = rep(NA, nrow(colors.mat))
    for(j.factor.it in 1:n.factors){
        j.factor = mix.factors[j.factor.it]
        j.factor.col = brewer.pal(9, mix.factors.pal[j.factor.it])[9]
        mixed.factor.val[j.factor.it] = paste(j.factor, j.factor.col, sep = '; color: ')
    }
    
    ##bind to colors.mat    
    colors.mat = cbind(colors.mat, mixed.factor.val, mixed.color, stringsAsFactors = FALSE)

    ##colnames
    mixed.factor = paste(mix.factors, collapse = '.')
    mixed.color.col = paste(mixed.factor, 'color', sep = '.')
    colnames(colors.mat)[(ncol(colors.mat) - 1):ncol(colors.mat)] = c(mixed.factor, mixed.color.col)

    return(colors.mat)
}

get.gene.colormap <- function(j.rpkm, genes2color.list, pad.frac = 0, sum.stat = 'median'){    
###gene-based colormap: column-wise median of Z-scaled normalized rows (gene expr scaled relative to other cells in the set)

    ##rm rows (genes) with zero variance
    const.rows = which(apply(j.rpkm, 1, function(x){var(x) == 0}))
    if(length(const.rows) != 0){
        warning('There are constant rows. These are removed.')
        j.rpkm = j.rpkm[setdiff(1:nrow(j.rpkm), const.rows), ]
    }

    ##scale rows
    j.rpkm = t(scale(t(j.rpkm)))

    ##loop gene sets
    gene.sets = names(genes2color.list)
    n.sets = length(gene.sets)
    for(j.set in 1:n.sets){

        j.set.name = gene.sets[j.set]
        
        ##get genes
        j.genes = genes2color.list[[j.set.name]][['genes']]
        j.palette = genes2color.list[[j.set.name]][['palette.name']]
        j.genes = intersect(j.genes, rownames(j.rpkm))
        
        ##color by median or mean expression
        if(sum.stat == 'mean'){
            j.weights = genes2color.list[[j.set.name]][['weights']]
            j.weights = j.weights[j.genes]
            j.stats = as.data.frame(as.matrix(apply(j.rpkm[j.genes, , drop = FALSE], 2, function(jcell.expr, w){sum(jcell.expr * w);}, w = j.weights)))
        }
        if(sum.stat == 'median'){
            j.stats = as.data.frame(as.matrix(apply(j.rpkm[j.genes, , drop = FALSE], 2, median)))
        }
        colnames(j.stats) = j.set.name

        
        ##get color
        j.genes2colormap = add.factor.color(j.stats, factor = j.set.name, discrete = FALSE, palette.name = j.palette, pad.frac = pad.frac)


        ##bind
        if(j.set == 1){
            genes2colormap = j.genes2colormap
        }else{
            genes2colormap = cbind(genes2colormap, j.genes2colormap)
        }

    }
    return(genes2colormap)
}

add.factor.color <- function(meta.mat, factor, discrete = TRUE, palette.name = 'Dark2', pad.frac = 0){
    library('RColorBrewer')
    
    if(discrete){
        factors = sort(unique(meta.mat[, factor]))
        n.factors = length(factors)
        if(any(factors == 'NA')){
            na.ind = which(factors == 'NA')
            factors = factors[c(setdiff(1:n.factors, na.ind), na.ind)]
        }
        if(n.factors < 3){
            factor2color.map = brewer.pal(3, palette.name)[1:n.factors]
        }else{
            if(n.factors > 8){
                factor2color.map = colorRampPalette(brewer.pal(8, palette.name))(n.factors)
            }else{
                factor2color.map = brewer.pal(n.factors, palette.name)[1:n.factors]
            }
        }
        names(factor2color.map) = factors
        factor.color = factor2color.map[as.character(meta.mat[, factor])]
    }else{ #continuous

        ##color palette (e.g. 'Reds', 'Blues', ...)
        n.cols = 30
        col.pal = colorRampPalette(brewer.pal(9, palette.name)[3:9])(n.cols)

        ##pad end with the end color to make the bright portion smaller
        if(pad.frac > 0){
            ##pad ends of palette with the "end" values
            pad.len = (pad.frac * n.cols) / (1 - pad.frac)
            col.pal = c(col.pal, rep(col.pal[length(col.pal)], pad.len))
        }
        n.cols = length(col.pal)
        
        ##scale the range to 1:n.cols as to get the index in the col.pal
        if(is.vector(meta.mat)){
            vals = meta.mat
        }else{
            vals = meta.mat[, factor]
        }
        min.val = min(vals, na.rm = TRUE)
        max.val = max(vals, na.rm = TRUE)
        range.val = max.val - min.val
        rescaled = (vals - min.val) / range.val
        col.pal.ind = round((rescaled * (n.cols - 1)) + 1)

        #get colors
        factor.color = col.pal[col.pal.ind]
    }

    #add column to meta.mat if it doesn't already exist
    factor.color.name = paste(factor, 'color', sep = '.')
    if(length(which(colnames(meta.mat) == factor.color.name)) != 0){
        meta.mat[, factor.color.name] = factor.color
    }else{
        meta.mat = cbind(meta.mat, factor.color, stringsAsFactors = FALSE)
        colnames(meta.mat)[ncol(meta.mat)] = factor.color.name
    }

    return(meta.mat)
}

rm.zero.genes <- function(data.mat){
#Rm all zero rows
    
    zero.genes = which(apply(data.mat, 1, function(jrow){all(jrow == 0)}))
    data.mat = data.mat[setdiff(1:nrow(data.mat), zero.genes), ]

    return(data.mat)
}
