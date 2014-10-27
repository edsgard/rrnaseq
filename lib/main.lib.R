
get.meta <- function(mapstats.tab.file, meta.file, meta.tab.file){

    
    #####
    #Read mapstats
    #####
    mapstats = read.table(mapstats.tab.file, sep = '\t', header = TRUE, stringsAsFactors = FALSE, fill = TRUE)
    colnames(mapstats) = gsub('\\.*$', '', colnames(mapstats))
    
    colnames(mapstats) = sub('Reads', 'n.reads', colnames(mapstats))
    colnames(mapstats) = sub('Readlength', 'read.len', colnames(mapstats))
    colnames(mapstats) = sub('Uniqely.mapping', 'uniq.mapped.prc', colnames(mapstats))
    colnames(mapstats) = sub('Unmapped.reads', 'unmapped.prc', colnames(mapstats))
    colnames(mapstats) = sub('Multimapping.reads', 'multimapping.prc', colnames(mapstats))    
    colnames(mapstats) = tolower(colnames(mapstats))

    #set rownames to sample
    nodata.ind = grep('no data for', mapstats[, 'n.reads'])
    if(length(nodata.ind) != 0){
        nodata.samples = gsub('no data for ', '', mapstats[nodata.ind, 'n.reads'])
        mapstats[nodata.ind, 'sample'] = nodata.samples
        mapstats[nodata.ind, 'n.reads'] = NA
        mapstats[, 'n.reads'] = as.integer(mapstats[, 'n.reads'])
    }
    rownames(mapstats) = mapstats[, 'sample']
        
    meta.mat = mapstats

    #ensure that coverage cols are numeric
    meta.mat[, 'coverage.first.5'] = as.numeric(meta.mat[, 'coverage.first.5'])
    meta.mat[, 'mid.5'] = as.numeric(meta.mat[, 'mid.5'])
    meta.mat[, 'last.5'] = as.numeric(meta.mat[, 'last.5'])

    
    ########
    #Add cols
    ########

    #group = gsub('[0-9]*', '', meta.mat[, 'sample'])
    #nr = gsub('[A-Z]*', '', meta.mat[, 'sample'])
    #id = paste(group, nr, sep = '.')
        
    #id
    #meta.mat = cbind(meta.mat, group, id, stringsAsFactors = FALSE)
    id = meta.mat[, 'sample']
    meta.mat = cbind(meta.mat, id, stringsAsFactors = FALSE)
    rownames(meta.mat) = meta.mat[, 'id']
           
    #Add colormap cols
    meta.mat = add.factor.color(meta.mat, 'id')
    #meta.mat = add.factor.color(meta.mat, 'group')
        
    nostrat = rep('nostrat', nrow(meta.mat))
    nostrat.color = rep('blue', nrow(meta.mat))
    meta.mat = cbind(meta.mat, nostrat, nostrat.color, stringsAsFactors = FALSE)
    
    #order on id
    meta.mat = meta.mat[order(meta.mat[, 'id']), ]
    
    #Dump
    saveRDS(meta.mat, file = meta.file)
    write.table(meta.mat, quote = FALSE, col.names = TRUE, row.names = FALSE, sep = '\t', file = meta.tab.file)

}

get.expr <- function(rpkm.files, rpkm.rsync.out.dir, rpkm.cloud.out.dir, meta.file = NA, meta.sample.col = 'sample', rpkm.files.listed = FALSE){


    #Get rpkmforgenes filenames
    rpkm.files = strsplit(rpkm.files, ',')[[1]]

    #If rpkm.files.listed set to TRUE, then the rpkm.files is a file containing a list of rpkm files to be read
    if(rpkm.files.listed){
        rpkm.files = read.table(rpkm.files, stringsAsFactors = FALSE)[, 1]
    }
    
    #Read RPKM (and count) files
    rpkm.counts.list = read.rpkm(rpkm.files)

    #Subset/crosscheck against cells present in meta.mat
    if(!is.na(meta.file)){
        
        meta.mat = readRDS(meta.file)    
        file.samples = colnames(rpkm.counts.list[['rpkm']])
        sample.subset = intersect(file.samples, meta.mat[, meta.sample.col])
        print(length(sample.subset)) #
        rpkm.counts.list = lapply(rpkm.counts.list, function(jdata, sample.subset){jdata[, sample.subset]}, sample.subset = sample.subset)

        #Change names: sample -> id    
        rpkm.counts.list = lapply(rpkm.counts.list, function(jdata, meta.mat){rownames(meta.mat) = meta.mat[, meta.sample.col]; colnames(jdata) = meta.mat[colnames(jdata), 'id']; return(jdata);}, meta.mat = meta.mat)
    }
    
    #Check
    lapply(rpkm.counts.list, dim)

    
    ###
    #Dump
    ###
        
    #Get all filenames
    rds.files = get.rpkm.files(rpkm.cloud.out.dir, rpkm.rsync.out.dir)
    tab.files = get.tab.files(rpkm.rsync.out.dir)

    #Write
    for(j.rds in names(rds.files)){
        saveRDS(rpkm.counts.list[[j.rds]], file = rds.files[j.rds])
    }

    for(j.tab in names(tab.files)){
        write.table(rpkm.counts.list[[j.tab]], quote = F, sep = '\t', col.names = NA, file = tab.files[j.tab])
    }
        
}

mapstats <- function(meta.file, counts.file, res.dir, strat.factor = 'nostrat', cex = 0.8, qc.meta.file, mapped.prc.cutoff, nreads.cutoff, filter.bool, plot.bool){

    library('RColorBrewer')

    #Hardcoded params
    pdf.w = 20
    pdf.h = 20
    qc.meta.tab.file = paste(sub('\\.rds$', '', qc.meta.file), 'tab', sep = '.')
    
    #Read data
    meta.mat = readRDS(meta.file)
    counts.mat = readRDS(counts.file)
    print(dim(meta.mat)) #

    #Error-check
    shared.samples = intersect(colnames(counts.mat), meta.mat[, 'id'])
    n.shared.samples = length(shared.samples)
    n.count.samples = ncol(counts.mat)
    n.meta.samples = length(unique(meta.mat[, 'id']))
    if(n.shared.samples != n.count.samples){
        warning(sprintf("The sample names in the count matrix (column names) do not fully agree with the sample id column in the matrix with the sample meta-information.\n n.data.samples: %i\n n.meta.samples: %i\n n.shared.samples: %i", n.count.samples, n.meta.samples, n.shared.samples))
    }
    
    #Subset counts.mat
    counts.mat = counts.mat[, shared.samples]
    
    #Colormaps
    strat.factor.color = paste(strat.factor, 'color', sep = '.')
    factor2color.map = unique(meta.mat[, c(strat.factor, strat.factor.color)])

    #Create dirs
    dir.create(res.dir, recursive = TRUE)

    
    if(plot.bool){
        
        ##################
        #Number of reads
        ##################
    
        #histogram
        pdf(file = file.path(res.dir, 'n.reads.dens.pdf'))
        plot(density(meta.mat[, 'n.reads']), xlab = 'n.reads', main= '')
        dev.off()

        #histogram log10
        pdf(file = file.path(res.dir, 'n.reads.dens.pdf'))
        plot(density(log10(meta.mat[, 'n.reads'])), xlab = 'log10(n.reads)', main= '', col = meta.mat[, strat.factor.color])
        #abline(v = log10(nreads.cutoff))
        dev.off()
        

        #barplot    
        pdf(file = file.path(res.dir, 'n.reads.bar.pdf'), width = pdf.w, height = pdf.h)
        bp = barplot(meta.mat[, 'n.reads'], axes = FALSE, axisnames = FALSE, ylab = '# of reads', main = '# of reads', col = meta.mat[, strat.factor.color], border = NA)
        axis(2)
        axis(1, at = bp, labels = meta.mat[, 'id'], cex.axis = cex, las = 2)
        #abline(h = nreads.cutoff)
        legend('topleft', legend = factor2color.map[, strat.factor], col = factor2color.map[, strat.factor.color], lty = 1)
        dev.off()


    
        ###
        #Percent mapped reads
        ###
        pdf(file = file.path(res.dir, 'percent.mapped.bar.pdf'), width = pdf.w, height = pdf.h)
        bp = barplot(as.matrix(t(meta.mat[, c('uniq.mapped.prc', 'multimapping.prc', 'unmapped.prc')])), beside=FALSE, axes = FALSE, axisnames = FALSE, ylab = '%', col = c('green', 'blue', 'red'), legend.text = c('Uniq', 'Multi', 'Unmapped'), args.legend = list(x = 37, y = 115), border = NA)
        axis(2, at = seq(0, 100, 10), labels = seq(0, 100, 10))
        axis(1, at = bp, labels = meta.mat[, 'id'], cex.axis = cex, las = 2)
        dev.off()    

        pdf(file = file.path(res.dir, 'percent.mapped.dens.pdf'))
        plot(density(meta.mat[, 'uniq.mapped.prc']), main = '', xlab = '% uniqely mapped reads', ylab = 'Sample density')
        dev.off()    
    

        ###
        #Number of uniquely mapped reads
        ###
        n.reads.uniq.mapped = meta.mat[, 'n.reads'] * (meta.mat[, 'uniq.mapped.prc'] / 100)
        meta.mat = cbind(meta.mat, n.reads.uniq.mapped)

    
        #barplot
        pdf(file = file.path(res.dir, 'n.reads.uniq.mapped.bar.pdf'), width = pdf.w, height = pdf.h)
        bp = barplot(meta.mat[, 'n.reads.uniq.mapped'], axes = FALSE, axisnames = FALSE, ylab = '# of reads', main = '# of reads', col = meta.mat[, strat.factor.color], border = NA)
        axis(2)
        axis(1, at = bp, labels = meta.mat[, 'id'], cex.axis = cex, las = 2)
        abline(h = nreads.cutoff)
        legend('topleft', legend = factor2color.map[, strat.factor], col = factor2color.map[, strat.factor.color], lty = 1)
        dev.off()

        #density
        pdf(file = file.path(res.dir, 'n.reads.uniq.mapped.dens.pdf'))
        plot(density(log10(meta.mat[, 'n.reads.uniq.mapped'])), xlab = 'log10(n.reads)', main= '', col = meta.mat[, strat.factor.color])
        dev.off()
    
    
    
        ####
        #Gene body coverage
        ####
    
        #boxplot for each stratum
        cov.cols = c('first.5', 'mid.5', 'last.5')

        #exclude samples where cov.cols missing
        rm.ind = which(is.na(meta.mat[, cov.cols[2]]))
        meta.mat.filt = meta.mat[setdiff(1:nrow(meta.mat), rm.ind), ]

        #split by stratum
        meta.factor = split(meta.mat.filt, meta.mat.filt[, strat.factor])
        
        factors = names(meta.factor)
        n.factors = length(factors)
        pdf(file = file.path(res.dir, 'genebody.cov.pdf'))
        
        if(n.factors >=2){      
            n.rows = ceiling(n.factors / 2)
            par(mfrow = c(n.rows, 2))
        }else{
            par(mfrow = c(1, 1))
        }
        
        par(mar = c(3,3,2,2)) #c(bottom, left, top, right), c(5, 4, 4, 2)
        for(j.factor in factors){
            j.meta.factor = meta.factor[[j.factor]]
            colnames(j.meta.factor) = sub('coverage.', '', colnames(j.meta.factor))
            if(nrow(j.meta.factor) > 0){
                boxplot(j.meta.factor[, cov.cols], main = j.factor, ylab = 'coverage')
            }
        }
        dev.off()

    

        #######
        #Number of reads mapped to gene annotation
        ######
        counts.mat.sample = colSums(counts.mat)

        #barplot
        pdf(file = file.path(res.dir, 'n.reads.annotmapped.barplot.pdf'), width = pdf.w, height = pdf.h)
        bp = barplot(counts.mat.sample, axes = FALSE, axisnames = FALSE, ylab = '# of reads', main = '', col = meta.mat[, strat.factor.color], border = NA)
        axis(2)
        axis(1, at = bp, labels = meta.mat[, 'id'], cex.axis = cex, las = 2)
        abline(h = nreads.cutoff)
        #legend('topleft', legend = factor2color.map[, strat.factor], col = factor2color.map[, strat.factor.color], lty = 1)
        dev.off()

        #density
        pdf(file = file.path(res.dir, 'n.reads.annotmapped.dens.pdf'))
        plot(density(log10(counts.mat.sample)), xlab = 'log10(# of reads)', ylab = 'Density: # of samples', main= '', col = meta.mat[, strat.factor.color])
        dev.off()


    
        ####
        #Fraction of reads mapped to annotation
        ####
        n.reads = meta.mat[, 'n.reads.uniq.mapped']

        #TBD: this should be wrong!?: frac.refseq = counts.mat.sample / (counts.mat.sample + n.reads)
        frac.refseq = counts.mat.sample / n.reads
    
        #barplot
        pdf(file = file.path(res.dir, 'fraction.annotmapped.barplot.pdf'), width = pdf.w, height = pdf.h)
        bp = barplot(frac.refseq, axes = FALSE, axisnames = FALSE, ylab = 'Fraction of uniquely mapped reads mapping to annotation', main = '', col = meta.mat[, strat.factor.color], border = NA)
        axis(2)
        axis(1, at = bp, labels = meta.mat[, 'id'], cex.axis = cex, las = 2)
        abline(h = nreads.cutoff)
        #legend('topleft', legend = factor2color.map[, strat.factor], col = factor2color.map[, strat.factor.color], lty = 1)
        dev.off()

        #density
        pdf(file = file.path(res.dir, 'fraction.refseq.dens.pdf'))
        plot(density(frac.refseq), xlab = 'Fraction of uniqely mapped reads mapping to annotation', main= '', col = meta.mat[, strat.factor.color])
        dev.off()
    }

    if(filter.bool){

        
        #Get qc.meta.mat, create new data-structure if file doesn't exist
        qc.meta.mat = get.qc.df(qc.meta.file, meta.mat)

        ###
        #Prc mapped reads
        ###
        qc.col = 'prc.mapped.reads'
        
        #Get failing samples
        fail.samples = meta.mat[which(meta.mat[, 'uniq.mapped.prc'] < mapped.prc.cutoff), 'id']
        length(fail.samples)

        #Create new column if it doesn't exist
        qc.meta.mat = add.qc.col(qc.meta.mat, qc.col)

        #Set filter
        qc.meta.mat[fail.samples, qc.col] = 1

        ###
        #Number of reads
        ###
        qc.col = 'n.reads'

        #Get failing samples
        fail.samples = meta.mat[which(meta.mat[, 'n.reads'] < nreads.cutoff), 'id']
        length(fail.samples) #

        #Create new column if it doesn't exist
        qc.meta.mat = add.qc.col(qc.meta.mat, qc.col)

        #Set filter
        qc.meta.mat[, qc.col] = 0
        qc.meta.mat[fail.samples, qc.col] = 1
        
        #Dump
        saveRDS(qc.meta.mat, file = qc.meta.file)
        write.table(qc.meta.mat, quote = F, col.names = NA, sep = '\t', file = qc.meta.tab.file)
    }
}

expr.dhist <- function(meta.mat, rpkm, exprdist.pdf, strat.factor, log.fcn = get('log10'), strat.plot = FALSE){

    library('RColorBrewer')
    
    #Create output dir    
    dir.create(dirname(exprdist.pdf), showWarnings = FALSE, recursive = TRUE)

    #Log
    if(is.function(log.fcn)){
        rpkm = log.fcn(rpkm)
        print('Logging')
    }

    #Plot
    pdf(file = exprdist.pdf)
    if(strat.plot){
        plot.strata.expr.dist(rpkm, meta.mat, strat.factor)
    }else{
        plot.sample.expr.dist(rpkm, meta.mat, strat.factor)
    }
    dev.off()    
}

gene.filter <- function(meta.file, rpkm.file, rpkm.postqc.file, gene2nsamples.pdf, strat.factor, rpkm.cutoff, n.cells.cutoff, n.strata.cutoff, nstrata.strat.factor, nstrata.ncells.cutoff, nstrata.nstrata.cutoff, filter.bool, plot.bool){
    
    #################
    #Gene filter: Number of cells a gene is expressed in
    #################
    
    #plot
    if(plot.bool){

        #Create output dir    
        dir.create(dirname(gene2nsamples.pdf), showWarnings = FALSE, recursive = TRUE)

        pdf(file = gene2nsamples.pdf)
        plot.gene2nsamples.expr(rpkm, meta.mat, strat.factor, rpkm.cutoff, n.cells.cutoff)
        dev.off()
    }
    
    #Filter genes
    if(filter.bool){
        print(nrow(rpkm))
        rpkm = gene2nsamples.expr.filter(rpkm, n.cells.cutoff, rpkm.cutoff, strat.factor, meta.mat, n.strata.cutoff, nstrata.strat.factor, nstrata.ncells.cutoff, nstrata.nstrata.cutoff)
        print(nrow(rpkm))

        #Dump
        saveRDS(rpkm, file = rpkm.postqc.file)
    }
}

sample2ngenes.expr <- function(rpkm, meta.mat, sample2ngenes.pdf, log.rpkm.cutoff, n.genes.cutoff, qc.meta.mat, filter.bool, plot.bool, meta.add, strat.factor = 'nostrat'){


    library('RColorBrewer')

    #Params
    rpkm.cutoff = 10^(log.rpkm.cutoff)
    
    
    ####################
    #Sample QC: Sample2n.genes.expr
    ####################

    #Get number of expr genes per sample
    n.genes.expr = apply(rpkm, 2, function(jsample, rpkm.cutoff){length(which(jsample >= rpkm.cutoff))}, rpkm.cutoff)
    
    #Plot
    if(plot.bool){

        #subset meta.mat
        meta.mat.new = meta.mat[colnames(rpkm), ]
        
        #add column if it doesn't already exist
        meta.mat.new = add.qc.col(meta.mat.new, 'n.genes.expr', def.val = NA)
        
        #set values
        meta.mat.new[names(n.genes.expr), 'n.genes.expr'] = n.genes.expr        

        #densities per stratum
        ngenes.dens.list = tapply(meta.mat.new[, 'n.genes.expr'], meta.mat.new[, strat.factor], density)

        #colormap
        strat.factor.col = paste(strat.factor, '.color', sep = '')
        strat.factor2color.map = unique(meta.mat.new[, c(strat.factor, strat.factor.col)])
        rownames(strat.factor2color.map) = strat.factor2color.map[, strat.factor]
        strat.factor2color.map = strat.factor2color.map[order(strat.factor2color.map[, strat.factor]), ]
        
        #Create output dir
        dir.create(dirname(sample2ngenes.pdf), showWarnings = FALSE, recursive = TRUE)

        #get ranges
        xlim = range(unlist(lapply(ngenes.dens.list, '[[', 'x'))) #xlim = c(-100, 100)
        ylim = range(unlist(lapply(ngenes.dens.list, '[[', 'y')))

        #strata
        strata = names(ngenes.dens.list)
        n.strata = length(strata)
        
        #plot
        pdf(file = sample2ngenes.pdf)
        j.stratum.it = 1
        j.stratum = strata[j.stratum.it]
        plot(ngenes.dens.list[[j.stratum]], xlab = 'n.genes', main = '', ylab = 'Density: # of samples', col = strat.factor2color.map[j.stratum, strat.factor.col], xlim = xlim, ylim = ylim)
        for(j.stratum.it in 2:n.strata){
            j.stratum = strata[j.stratum.it]
            lines(ngenes.dens.list[[j.stratum]], col = strat.factor2color.map[j.stratum, strat.factor.col])
        }
        legend('topright', legend = strat.factor2color.map[, strat.factor], col = strat.factor2color.map[, strat.factor.col], lty = 1)
        abline(v = n.genes.cutoff, col = 'grey', lty = 2)
        dev.off()
    }

    
    #Set filter
    if(filter.bool){

        qc.col = 'n.genes.expr'
        
        #Get failing samples
        fail.ind = which(n.genes.expr < n.genes.cutoff)
        fail.samples = names(fail.ind)    
        length(fail.samples) #
        #cat(fail.samples)
        
        #Create new column if it doesn't exist
        qc.meta.mat = add.qc.col(qc.meta.mat, qc.col)

        #Set filter
        qc.meta.mat[, qc.col] = 0
        qc.meta.mat[fail.samples, qc.col] = 1
    }

    #add n.gene.expr factor to meta.mat
    if(meta.add){
        
        #add column if it doesn't already exist
        meta.mat = add.qc.col(meta.mat, 'n.genes.expr', def.val = NA)
        
        #set values
        meta.mat[names(n.genes.expr), 'n.genes.expr'] = n.genes.expr        
        
        #Add colormap cols
        meta.mat = add.factor.color(meta.mat, 'n.genes.expr', discrete = FALSE)
    }
    
    res.list = list(n.genes.expr = n.genes.expr, qc.meta.mat = qc.meta.mat, meta.mat = meta.mat)
    return(res.list)
}

sampledist.heatmap <- function(meta.mat, rpkm, sample.heatmap.pdf = 'sample.heatmap.pdf', strat.factor = 'nostrat', cor.meth = 'spearman', cex.sample = 0.7, cor.res.list = NA, ncores = 1, nblocks = 1, rm.diag = TRUE, log.fcn = get('log10'), ...){

    library('RColorBrewer')

    #Create output dir    
    dir.create(dirname(sample.heatmap.pdf), showWarnings = FALSE, recursive = TRUE)

    #Log
    if(is.function(log.fcn)){
        rpkm = log.fcn(rpkm + 1)
        print('Logging')
    }
    
    #######################
    #Sample-distance heatmap
    #######################

    if(is.na(cor.res.list)){

        #rm cols with constant variance
        const.cols = which(apply(rpkm, 2, function(x){var(x) == 0}))
        if(length(const.cols) != 0){
            warning('There were constant columns. These were removed.')
            rpkm = rpkm[, setdiff(1:ncol(rpkm), const.cols)]
        }
        
        #Pairwise dist btw samples
        col.cor = bigcor.par(rpkm, use="pairwise.complete.obs", method = cor.meth, ncores = ncores, nblocks = nblocks)
    
        #get dist and hclust 
        col.dist.mat = as.dist(((col.cor * -1) + 1) / 2)
        col.clust = hclust(col.dist.mat)
        
        cor.res.list = list(col.cor = col.cor, col.clust = col.clust)
    }else{

        if(length(which(names(cor.res.list) == 'col.cor')) != 0){
            col.cor = cor.res.list[['col.cor']]

            #get col.clust
            if(length(which(names(cor.res.list) == 'col.clust')) != 0){
                col.clust = cor.res.list[['col.clust']]
            }else{            
                col.dist.mat = as.dist(((col.cor * -1) + 1) / 2)
                col.clust = hclust(col.dist.mat)

                #store result
                cor.res.list[['col.clust']] = col.clust
            }
        }else{
            
            #rm cols with constant variance
            const.cols = which(apply(rpkm, 2, function(x){var(x) == 0}))
            if(length(const.cols) != 0){
                warning('There were constant columns. These were removed.')
                rpkm = rpkm[, setdiff(1:ncol(rpkm), const.cols)]
            }
        
            #Pairwise dist btw samples
            col.cor = bigcor.par(rpkm, use="pairwise.complete.obs", method = cor.meth, ncores = ncores, nblocks = nblocks)

            #get dist and hclust 
            col.dist.mat = as.dist(((col.cor * -1) + 1) / 2)
            col.clust = hclust(col.dist.mat)

            #add to list
            cor.res.list = c(cor.res.list, list(col.cor = col.cor, col.clust = col.clust))            
        }                
    }
    
    if(rm.diag){
    	diag(col.cor) = median(col.cor)
    }
    
    #Plot
    heatmap.pdf = file.path(dirname(sample.heatmap.pdf), sub('heatmap', paste('heatmap.', cor.meth, sep = ''), basename(sample.heatmap.pdf)))
    pdf(file = heatmap.pdf, width = 10, height = 10)
    plot.heatmap.sample(col.cor, col.clust, meta.mat, strat.factor, cexRow = cex.sample, cexCol = cex.sample, ...)
    dev.off()


    return(cor.res.list)
}

sampledist.boxplot <- function(meta.mat, rpkm, sample.cor.pdf, qc.meta.file, cor.meth, max.cor.cutoff, strat.factor = 'nostrat', cex.axis, filter.bool, plot.bool, cor.res.list = NA, ...){


    library('RColorBrewer')

    ###
    #Params
    ###
    qc.meta.tab.file = paste(sub('\\.rds$', '', qc.meta.file), 'tab', sep = '.')
    sample.cor.pdf.dir = dirname(sample.cor.pdf)
    sample.cor.pdf.base = basename(sample.cor.pdf)
    cor.boxplot.pdf = file.path(sample.cor.pdf.dir, sub('cor', paste('cor', cor.meth, 'boxplot', sep = '.'), sample.cor.pdf.base))
    cor.max.dens.pdf = file.path(sample.cor.pdf.dir, sub('cor', paste('cor', cor.meth, 'max.dens', sep = '.'), sample.cor.pdf.base))

    
    ###
    #Filter data
    ###
    meta.mat.filt = meta.mat[colnames(rpkm), ]
    
    
    ######################
    #Sample-to-sample correlation boxplot
    ######################
    #Plot boxplot of pairwise correlations for each sample

    if(is.na(cor.res.list)){

        #rm cols with constant variance
        const.cols = which(apply(rpkm, 2, function(x){var(x) == 0}))
        if(length(const.cols) != 0){
            warning('There were constant columns. These were removed.')
            rpkm = rpkm[, setdiff(1:ncol(rpkm), const.cols)]
        }

        #Pairwise dist btw samples
        col.cor = bigcor.par(rpkm, use="pairwise.complete.obs", method = cor.meth, ...)

        #store result
        cor.res.list = list(col.cor = col.cor)

    }else{
        if(length(which(names(cor.res.list) == 'col.cor')) != 0){
            col.cor = cor.res.list[['col.cor']]
        }else{
            
            #rm cols with constant variance
            const.cols = which(apply(rpkm, 2, function(x){var(x) == 0}))
            if(length(const.cols) != 0){
                warning('There were constant columns. These were removed.')
                rpkm = rpkm[, setdiff(1:ncol(rpkm), const.cols)]
            }

            #Pairwise dist btw samples
            col.cor = bigcor.par(rpkm, use="pairwise.complete.obs", method = cor.meth, ...)

            #store result
            cor.res.list = c(cor.res.list, list(col.cor = col.cor))
        }
    }
    
    #Max correlation to another sample
    col.cor.noone = col.cor
    col.cor.noone[which(col.cor.noone > 0.99)] = NA
    max.cor = apply(col.cor.noone, 2, max, na.rm = TRUE)
    max.cor = sort(max.cor)

    #filter meta.mat
    meta.mat.filt = meta.mat[names(max.cor), ]

    
    if(plot.bool){

        #Create output dir    
        dir.create(dirname(sample.cor.pdf), showWarnings = FALSE, recursive = TRUE)

        ###
        #Boxplots of all pairwise cors per sample
        ###

        #Order by median correlation
        median.cor = apply(col.cor, 2, median)
        median.cor = sort(median.cor)

        #plot
        pdf(cor.boxplot.pdf)
        boxplot(col.cor[, names(median.cor)], las = 2, cex.axis = cex.axis)
        dev.off()

        ###
        #Plot histogram of max.cor
        ###

        #split by strata
        strat2sample.list = tapply(meta.mat.filt[, 'id'], meta.mat.filt[, strat.factor], unique)
        strata = names(strat2sample.list)
        n.strata = length(strata)

        #densities per stratum
        strat2maxcor.dens.list = lapply(strat2sample.list, function(jstrat.samples, max.cor){density(max.cor[jstrat.samples])}, max.cor = max.cor)

        #get ranges
        xlim = range(unlist(lapply(strat2maxcor.dens.list, '[[', 'x')))
        ylim = range(unlist(lapply(strat2maxcor.dens.list, '[[', 'y')))

        #colormap
        strat.factor.col = paste(strat.factor, '.color', sep = '')
        strat.factor2color.map = unique(meta.mat.filt[, c(strat.factor, strat.factor.col)])
        strat.factor2color.map = strat.factor2color.map[order(strat.factor2color.map[, strat.factor]), ]
        rownames(strat.factor2color.map) = strat.factor2color.map[, strat.factor]

        #plot
        pdf(cor.max.dens.pdf)
        j.stratum.it = 1
        j.stratum = strata[j.stratum.it]
        plot(strat2maxcor.dens.list[[j.stratum]], main = '', xlab = paste('max pairwise corr: ', cor.meth, sep = ''), ylab = 'Density: # of samples', col = strat.factor2color.map[j.stratum, strat.factor.col], xlim = xlim, ylim = ylim)
        for(j.stratum.it in 2:n.strata){
            j.stratum = strata[j.stratum.it]
            lines(strat2maxcor.dens.list[[j.stratum]], col = strat.factor2color.map[j.stratum, strat.factor.col])
        }
        legend('topright', legend = strat.factor2color.map[, strat.factor], col = strat.factor2color.map[, strat.factor.col], lty = 1)
        abline(v = max.cor.cutoff, col = 'grey', lty = 2)        
        dev.off()
    }    
    
    #Get fail.samples
    if(filter.bool){

        qc.col = paste(cor.meth, 'cor', sep = '.')
        
        #Get failing samples
        fail.samples = names(max.cor)[which(max.cor < max.cor.cutoff)]
        length(fail.samples) #16

        #Get qc.meta.mat, create new data-structure if file doesn't exist
        qc.meta.mat = get.qc.df(qc.meta.file, meta.mat)

        #Create new column if it doesn't exist
        qc.meta.mat = add.qc.col(qc.meta.mat, qc.col)

        #Set filter
        qc.meta.mat[, qc.col] = 0
        qc.meta.mat[fail.samples, qc.col] = 1

        #Dump
        saveRDS(qc.meta.mat, file = qc.meta.file)
        write.table(qc.meta.mat, quote = F, col.names = NA, sep = '\t', file = qc.meta.tab.file)
    }
    
    return(cor.res.list)
}

sample.hclust <- function(meta.file, rpkm.file, sample.hclust.pdf, cor.meth, n.boot, strat.factor, ind.factor, cex, cor.res.list = NA, ...){


    library('RColorBrewer')

    ###
    #Params
    ###
    sample.hclust.pdf.dir = dirname(sample.hclust.pdf)
    sample.hclust.pdf.base = basename(sample.hclust.pdf)
    clust.pdf = file.path(sample.hclust.pdf.dir, sub('hclust.', paste('hclust.', cor.meth, '.', sep = ''), sample.hclust.pdf.base))
    col.clust.pdf = file.path(sample.hclust.pdf.dir, sub('hclust.', paste('hclust.colors.', cor.meth, '.', sep = ''), sample.hclust.pdf.base))
    strat.factor.col = paste(strat.factor, 'color', sep = '.')
    ind.factor.col = paste(ind.factor, 'color', sep = '.')
    pdf.w = 20
    pdf.h = 20
    
    ###
    #Load data
    ###
    meta.mat = readRDS(meta.file)
    rpkm = readRDS(rpkm.file)
    dim(rpkm) #Refseq: 24249 x n.samples
    dim(meta.mat) #n.samples x n.cols

    #Create output dir    
    dir.create(dirname(sample.hclust.pdf), showWarnings = FALSE, recursive = TRUE)

    
    if(n.boot == 0){
        library('moduleColor') #depends on lib 'impute' which is not available in R>=3. Install impute via biocLite.
        library('graphics') #dendrapply

        #########################
        #Sample hierarchical clustering w/o bootstrap
        #########################
        #colored labels in dendrogram: plotColoredClusters {ClassDiscovery} (#not available for R>=3. also e.g. (PerturbationClusterTest)); dendextend; dendrapply{graphics}    


        if(is.na(cor.res.list)){
            
            #rm cols with constant variance
            const.cols = which(apply(rpkm, 2, function(x){var(x) == 0}))
            if(length(const.cols) != 0){
                warning('There were constant columns. These were removed.')
                rpkm = rpkm[, setdiff(1:ncol(rpkm), const.cols)]
            }

            #Pairwise dist btw samples
            col.cor = bigcor.par(rpkm, use="pairwise.complete.obs", method = cor.meth, ...)
    		
            #get dist and hclust 
            col.dist.mat = as.dist(((col.cor * -1) + 1) / 2)
            col.clust = hclust(col.dist.mat)
            
            cor.res.list = list(col.cor = col.cor, col.clust = col.clust)
        }else{
            if(length(which(names(cor.res.list) == 'col.cor')) != 0){

                #get cor.res
                col.cor = cor.res.list[['col.cor']]

                #get col.clust
                if(length(which(names(cor.res.list) == 'col.clust')) != 0){
                    col.clust = cor.res.list[['col.clust']]
                }else{            
                    col.dist.mat = as.dist(((col.cor * -1) + 1) / 2)
                    col.clust = hclust(col.dist.mat)

                    #store result
                    cor.res.list = c(cor.res.list, list(col.clust = col.clust))
                }
            }else{
            
                #rm cols with constant variance
                const.cols = which(apply(rpkm, 2, function(x){var(x) == 0}))
                if(length(const.cols) != 0){
                    warning('There were constant columns. These were removed.')
                    rpkm = rpkm[, setdiff(1:ncol(rpkm), const.cols)]
                }
        
                #Pairwise dist btw samples
                col.cor = bigcor.par(rpkm, use="pairwise.complete.obs", method = cor.meth, ...)
                
                #get dist and hclust 
                col.dist.mat = as.dist(((col.cor * -1) + 1) / 2)
                col.clust = hclust(col.dist.mat)

                #add to list
                cor.res.list = c(cor.res.list, list(col.cor = col.cor, col.clust = col.clust))            
            }
        }
        
        #Colormap
        d = dendrapply(as.dendrogram(col.clust), labelCol, meta.mat = meta.mat, strat.factor = strat.factor, cex = cex)

        #Plot cluster
        pdf(clust.pdf, width = pdf.w, height = pdf.h)
        plot(d, cex = cex)
        dev.off()

        #Plot labels
        stratum.color = as.character(meta.mat[colnames(col.cor), strat.factor.col])
        ind.color = as.character(meta.mat[colnames(col.cor), ind.factor.col])
        pdf(col.clust.pdf, width = 100, height = 7)
        plotHclustColors(col.clust, cbind(ind.color, stratum.color), cex.rowLabels = 10, rowLabels = '')
        dev.off()

    }else{
        
        ################
        #Sample hclust with pvclust (bootstrapped)
        ################

        if(cor.meth == 'spearman'){
            cor.dist.fcn = spear.cor.col.dist
        }

        #rm cols with constant variance
        const.cols = which(apply(rpkm, 2, function(x){var(x) == 0}))
        if(length(const.cols) != 0){
            warning('There were constant columns. These were removed.')
            rpkm = rpkm[, setdiff(1:ncol(rpkm), const.cols)]
        }
        
        #cluster
        pvclust.res = pvclust(rpkm, method.dist = cor.dist.fcn, nboot = n.boot)
        
        #dump
        if(is.na(cor.res.list)){
            cor.res.list = list(pvclust.res = pvclust.res)
        }else{
            cor.res.list = c(cor.res.list, list(pvclust.res = pvclust.res))
        }
        #saveRDS(pvclust.res, file = '/Volumes/Data/cloud/gdrive/work/rspd/code/my/rrnaseq/test/rqc/data/tmp.rds')


        #plot
        pdf(file = clust.pdf)
        plot.pvclust(pvclust.res)
        dev.off()
        #pvrect(pvclust.res, alpha = 0.95)
    }

    return(cor.res.list)
}

pca <- function(meta.mat, rpkm, pca.pdf = 'pca.pdf', plot.comp = c(1, 2), log.fcn = get('log10'), pc.cutoff = 0, filter.bool = FALSE, gt.bool = FALSE, qc.meta.file = 'qc.meta.rds', strat.factor = 'nostrat', point.cex = 0.4){
    
    library('RColorBrewer')
    
    ###
    #Params
    ###   
    pca.evar.pdf = file.path(dirname(pca.pdf), sub('pca\\.', 'pca.evar.', basename(pca.pdf)))
    qc.meta.tab.file = paste(sub('\\.rds$', '', qc.meta.file), 'tab', sep = '.')
    strat.factor.col = paste(strat.factor, 'color', sep = '.')
    

    #Create output dir    
    dir.create(dirname(pca.pdf), showWarnings = FALSE, recursive = TRUE)

    
    ######################
    #PCA
    ######################


    ###
    #Filter rows (typically genes) with constant values
    ###
    rm.rows = which(apply(rpkm, 1, function(jrow){if(var(jrow) == 0){const = TRUE}else{const = FALSE}; return(const);}))
    rpkm = rpkm[setdiff(1:nrow(rpkm), rm.rows), ]
    
    ###
    #Log-transform
    ###
    if(is.function(log.fcn)){
        rpkm = rpkm + 1e-10
        rpkm = log.fcn(rpkm)
    }
    
    #PCA
    pca = prcomp(t(rpkm), scale=T, center=T)
    pca.basis = pca[['x']]

    #Explained var
    e.var = pca[['sdev']]^2 / sum(pca[['sdev']]^2)
    names(e.var) = colnames(pca.basis)
    
    #Plot explained var
    pdf(file = pca.evar.pdf)
    barplot(e.var[1:30] * 100, las = 2, ylab = 'Explained variance (%)', col = 'blue')
    dev.off()

    n.comp = length(plot.comp)
    
    if(n.comp > 2){
        
        ###
        #Pairs scatter-plot
        ###
        pdf(pca.pdf)
        plot.pairs.pca(pca.basis[, plot.comp], e.var, meta.mat, factor.color = strat.factor, point.text.cex = point.cex)
        dev.off()
        
    }
    if(n.comp == 2){
                
        ###
        #Single scatter plot
        ###

        #Params
        pc.x = paste('PC', plot.comp[1], sep = '')
        pc.y = paste('PC', plot.comp[2], sep = '')
        pca.pdf = file.path(dirname(pca.pdf), sub('pca\\.', paste('pca.', pc.x, '_', pc.y, '.', sep = ''), basename(pca.pdf)))
        
        #Colormaps
        meta.mat = meta.mat[rownames(pca.basis), ]
        strat2col.map = unique(meta.mat[, c(strat.factor, strat.factor.col)])
        strat2col.map = strat2col.map[order(strat2col.map[, strat.factor]), ]
        strat.color = meta.mat[, strat.factor.col]
    
        #Plot
        pdf(file = pca.pdf)
        plot(pca.basis[, pc.x], pca.basis[, pc.y], col = strat.color, xlab = pc.x, ylab = pc.y, pch = 16)
        legend('topleft', legend = strat2col.map[, strat.factor], lty = 1, col = strat2col.map[, strat.factor.col])
        dev.off()

        #celltype2col.map = unique(meta.mat[, c('celltype', 'celltype.pch')])
        #celltype.pch = meta.mat[, 'celltype.pch']

        #plot(pca.basis[, pc.x], pca.basis[, pc.y], pch = celltype.pch, col = strat.color, xlab = pc.x, ylab = pc.y)
        #legend('topleft', legend = strat2col.map[, strat.factor], lty = 1, col = strat2col.map[, strat.factor.col])
        #legend('bottomleft', legend = celltype2col.map[, 'celltype'], pch = celltype2col.map[, 'celltype.pch'])
    }
    if(n.comp == 1){
            
        ###
        #Histogram, projected on a single PC-component
        ###

        #Params
        pc = paste('PC', plot.comp, sep = '')
        pca.pdf = file.path(dirname(pca.pdf), sub('pca\\.', paste('pca', pc, 'dhist.', sep = '.'), basename(pca.pdf)))

        #plot
        pdf(file = pca.pdf)
        hist(pca.basis[, pc])
        plot(density(pca.basis[, pc]), xlab = pc, main = '')
        abline(v = pc.cutoff)
        dev.off()

        ###
        #Filter based on a single PC
        ###
        if(filter.bool){

            qc.col = pc
            
            #Get fail.samples
            if(gt.bool == TRUE){
                pc.fail.samples = rownames(pca.basis)[which(pca.basis[, pc] >= pc.cutoff)]
            }else{
                pc.fail.samples = rownames(pca.basis)[which(pca.basis[, pc] < pc.cutoff)]
            }
            #cat(pc.fail.samples)
            length(pc.fail.samples)            
            
            #Get qc.meta.mat, create new data-structure if file doesn't exist
            qc.meta.mat = get.qc.df(qc.meta.file, meta.mat)

            #Create new column if it doesn't exist
            qc.meta.mat = add.qc.col(qc.meta.mat, pc)

            #Set filter
            qc.meta.mat[, qc.col] = 0
            qc.meta.mat[pc.fail.samples, pc] = 1

            #Dump
            saveRDS(qc.meta.mat, file = qc.meta.file)
            write.table(qc.meta.mat, quote = F, col.names = NA, sep = '\t', file = qc.meta.tab.file)
        }
    }
    
}

scatter <- function(data.file, scatter.pdf, samples, transpose, log.fcn){
    
    #Read data
    data.mat = readRDS(data.file)

    #Set default if samples not set
    if(is.na(samples)){samples = colnames(data.mat);}

    #Create output dir    
    dir.create(dirname(scatter.pdf), showWarnings = FALSE, recursive = TRUE)

    
    #############
    #Sample-sample scatter plot
    ##############

    #Transpose
    if(transpose){
        data.mat = t(data.mat)
    }

    #Log
    if(is.function(log.fcn)){
        data.mat = log.fcn(data.mat)
    }

    n.samples = length(samples)
    pdf(file = scatter.pdf, width = 10, height = 10)
    par(mfrow = c(n.samples, n.samples))    
    for(j.sample in samples){
        for(i.sample in samples){
            smoothScatter(data.mat[, i.sample], data.mat[, j.sample], xlab = i.sample, ylab = j.sample)
        }
    }
    dev.off()
}

rseq.heatmap <- function(data.mat, cor.meth = 'spearman', meta.mat = NA, strat.heat.factor = NA, log.fcn = log10, cor.res.list = NA, colorpad = FALSE, ColSideColors = NA, col.pal = NA, ncores = 1, nblocks = 1, cluster.col = TRUE, ...){ 
# cexRow = 0.2 #500 rows: 0.2
# cexCol = 0.2 #500 cols: 0.2
#Colorbars
# strat.heat.factors = c('stage', 'stage.ind', 'enrichment', 'tri.lin', 'bi.lin') 
# genes = rownames(rpkm)
# RowSideColors = feat.annot[genes, 'marker.color']        
# marker2color.map = unique(feat.annot[, c('marker', 'marker.color')])

    
    library('marray') #maPalette
    library('gplots')


    ###
    #Log
    ###
    if(is.function(log.fcn)){
        data.mat = data.mat + 1e-10
        data.mat = log.fcn(data.mat)
        print('Logging')
    }

    
    ##***
    ##Cluster cols
    ##***

    
    if(cluster.col == TRUE){
        ##rm constant cols
        data.mat = rm.const.vec(data.mat, row.rm = FALSE)
    
        ##pairwise dist and cluster
        cor.res.list = rseq.hclust(data.mat, cor.meth, cor.res.list, store.col = 'col', ncores, nblocks)

        ##get dendrogram
        col.dendro = as.dendrogram(cor.res.list[['col.hclust']])
        col.dendro = reorder(col.dendro, colMeans(data.mat, na.rm = TRUE))
    }else{
        col.dendro = NULL
    }

    
    ###
    #Cluster rows
    ###
    #rm constant rows
    data.mat = rm.const.vec(data.mat, col.rm = FALSE)
    
    #pairwise dist and cluster
    cor.res.list = rseq.hclust(t(data.mat), cor.meth, cor.res.list, store.col = 'row', ncores, nblocks)
    
    #get dendrogram
    row.dendro = as.dendrogram(cor.res.list[['row.hclust']])
    row.dendro = reorder(row.dendro, rowMeans(data.mat, na.rm = TRUE))


    ###
    #Set default colormap if col.pal == NA
    ###
    if(is.na(col.pal)){
        if(colorpad){
            pad.len = 20
    
            #pad red and blue to the ends of the pal    
            col.pal = maPalette(low = "blue", high = 'red', mid='white', k = 10)
            col.pal = c(rep(col.pal[1], pad.len), col.pal, rep(col.pal[length(col.pal)], pad.len))
        }else{
            #no padding    
            col.pal = maPalette(low = "blue", high = 'red', mid='white', k = 50)
            #col.pal = greenred(100) #alt: maPalette{marray)
        }
    }

    
    ###
    #Colorbars
    ###
    meta.col.bar = !is.na(meta.mat) & !is.na(strat.heat.factor) & is.na(ColSideColors)
    if(meta.col.bar){
        color.col = paste(strat.heat.factor, 'color', sep = '.')
        samples = colnames(data.mat)
        meta.mat = meta.mat[samples, ]
        ColSideColors = as.character(meta.mat[, color.col])
        sample2color.map = unique(meta.mat[, c(strat.heat.factor, color.col)])
        sample2color.map = sample2color.map[order(sample2color.map[, strat.heat.factor]), ]
    }


    ###
    #Plot
    ###
    if(is.na(ColSideColors)){
        heatmap.2(data.mat, Colv = col.dendro, Rowv = row.dendro, col = col.pal, trace = 'none', ...)
    }else{
        heatmap.2(data.mat, Colv = col.dendro, Rowv = row.dendro, col = col.pal, trace = 'none', ColSideColors = ColSideColors, ...)
    }
    if(meta.col.bar){
        legend('topright', legend = sample2color.map[, strat.heat.factor], col = sample2color.map[, color.col], lty = 1)
    }

    return(cor.res.list)
}

sample.filter <- function(data.file, qc.meta.file, out.file, filter.cols){

    
    #Read data
    data.mat = readRDS(data.file)
    qc.mat = readRDS(qc.meta.file)
    
    #Parse args
    filter.cols = strsplit(filter.cols, ',')[[1]]

    #Error-check
    qc.cols = colnames(qc.mat)
    if(length(setdiff(filter.cols, qc.cols)) != 0){
        stop('ERROR: A filter column was provided that is not present in qc.mat')
    }

    #Create output dir    
    dir.create(dirname(out.file), showWarnings = FALSE, recursive = TRUE)
    
    #Get samples to filter
    if(filter.cols == 'all'){
        filter.cols = colnames(qc.mat)
    }
    fail.samples = rownames(qc.mat)[which(apply(qc.mat[, filter.cols], 1, function(jrow){any(jrow == 1)}))]
    pass.samples = setdiff(colnames(data.mat), fail.samples)

    #Filter samples
    print(ncol(data.mat))
    data.mat = data.mat[, pass.samples]
    print(ncol(data.mat)) #

    #Dump
    saveRDS(data.mat, file = out.file)
}

ercc <- function(meta.mat, real.rpkm, real.counts, ercc.rpkm, ercc.counts, annot, out.dir, strat.factor, n.ercc.min, min.prop, min.rpkm, log.size.cutoff, ercc.rpkm.frac.adjust, cex, spike.vol){

    
    #Params
    strat.factor.color = paste(strat.factor, 'color', sep = '.')
    pdf.w = 20
    pdf.h = 20
    
    #Create output dir    
    dir.create(out.dir, showWarnings = FALSE, recursive = TRUE)

    #PDFs
    sampledist.pdf = file.path(out.dir, 'ercc.sampledist.pdf')
    ercc.pdf = file.path(out.dir, 'ercc.pdf')
    mean2cv.pdf = file.path(out.dir, 'mean2cv.pdf')
    real.mean2cv.pdf = file.path(out.dir, 'real.mean2cv.pdf')
    size.factor.hist.pdf = file.path(out.dir, 'size.factor.hist.pdf')
    fraction.ercc.barplot.pdf = file.path(out.dir, 'fraction.ercc.barplot.pdf')
    fraction.ercc.dens.pdf = file.path(out.dir, 'fraction.ercc.dens.pdf')    
    
    #Ensure that same set samples: real.prkm, ercc.rpkm and meta.mat (for example if real.rpkm has been sample-filtered)
    pass.samples = colnames(real.rpkm)
    ercc.rpkm = ercc.rpkm[, pass.samples]
    meta.mat = meta.mat[pass.samples, ]

    #Subset real.counts genes
    real.counts = real.counts[intersect(rownames(real.counts), rownames(real.rpkm)), ]

    
    ###
    #Convert concentration to number of added transcripts
    ###
    c = annot[, 'mix1.conc'] #attomol/uL
    n.attomol = spike.vol * c #attomol
    n.molecules = round(n.attomol * 1e-18 * 6.022 * 1e23)
    annot = cbind(annot, n.attomol, n.molecules)
    #print(annot[order(annot[, 'n.molecules']), ])
    rownames(annot) = annot[, 'ERCC.ID']

    
    ###
    #Filter samples
    ###
    print(dim(ercc.rpkm)) #92 n.samples
    
    #Rm samples with no ERCC in the rpkm annot (they didn't have ERCC spiked in to begin with either)
    ercc.ind = which(apply(ercc.rpkm, 2, function(j.ercc){length(which(!is.na(j.ercc))) > 0}))
    pass.samples = names(ercc.ind) 
    ercc.rpkm = ercc.rpkm[, pass.samples]
    print(ncol(ercc.rpkm)) #

    #Rm samples with reads mapped to less than n ERCC transcripts
    ercc.ind = which(apply(ercc.rpkm, 2, function(j.ercc){length(which(j.ercc != 0)) > n.ercc.min}))
    pass.samples = names(ercc.ind)
    ercc.rpkm = ercc.rpkm[, pass.samples]
    print(length(pass.samples)) #
        
    #For each cell, proportion of expressed ERCC transcripts > min.rpkm should be greater than min.prop
    ercc.ind = which(apply(ercc.rpkm, 2, function(j.ercc){j.ercc = j.ercc[which(j.ercc != 0)]; prop.expr = length(which(j.ercc > min.rpkm)) / length(j.ercc); return(prop.expr >= min.prop)}))
    pass.samples = names(ercc.ind)    
    print(length(pass.samples)) #

    #apply filter
    ercc.rpkm = ercc.rpkm[, pass.samples]    
    meta.mat = meta.mat[pass.samples, ]
    real.rpkm = real.rpkm[, pass.samples]

    
    ###
    #Filter samples on fraction of spike-in RNA
    ###
    ercc.sum = apply(ercc.rpkm, 2, sum)
    real.sum = apply(real.rpkm, 2, sum)
    ercc.frac = ercc.sum / (ercc.sum + real.sum)
    
    #calculate a size factor relative the median frac
    size.factor = ercc.frac / median(ercc.frac)
    log.size.factor = log10(size.factor)

    #Plot hist of size.factor
    #plot(density(size.factor))
    pdf(file = size.factor.hist.pdf)
    hist(log.size.factor, breaks = seq(min(log.size.factor) - 0.1 * abs(min(log.size.factor)), max(log.size.factor) + 0.1, 0.01), xlab = 'log10(size.factor)')
    dev.off()

    #Rm samples with more than ten-fold difference to median
    pass.samples = names(log.size.factor)[which(abs(log.size.factor) <= log.size.cutoff)]
    print(length(pass.samples)) #

    #apply filter
    ercc.rpkm = ercc.rpkm[, pass.samples]    
    meta.mat = meta.mat[pass.samples, ]
    real.counts = real.counts[, pass.samples]
    ercc.counts = ercc.counts[, pass.samples]        
    size.factor = size.factor[pass.samples]

    
    ###
    #Correct on fraction of spike-in RNA
    ###
    if(ercc.rpkm.frac.adjust){
        samples = names(size.factor)
        ercc.corr = as.data.frame(lapply(samples, function(j.sample){ercc.rpkm[, j.sample] / size.factor[j.sample]}))
        colnames(ercc.corr) = samples
    
        ercc.rpkm = ercc.corr
    }

    
    ###
    #Plot RPKM histogram
    ###
    
    ercc.log = log2(ercc.rpkm)
    n.samples = ncol(ercc.rpkm)
    densities = apply(ercc.log, 2, density)    
    xlim = range(unlist(lapply(densities, '[[', 'x')))
    ylim = range(unlist(lapply(densities, '[[', 'y')))
    
    batch2color.map = unique(meta.mat[, c(strat.factor, strat.factor.color)])
    batch2color.map = batch2color.map[order(batch2color.map[, strat.factor]), ]
    pdf(file = sampledist.pdf)
    j.sample = 1
    plot(densities[[j.sample]], col = meta.mat[j.sample, strat.factor.color], ylim = ylim, xlim = xlim, xlab = 'log2(RPKM)', main = '', xaxt = 'n')
    for(j.sample in 2:n.samples){
        lines(densities[[j.sample]], col = meta.mat[j.sample, strat.factor.color], ylim = ylim, xlim = xlim, xaxt = 'n')
    }
    legend('topright', legend = batch2color.map[, strat.factor], col = batch2color.map[, strat.factor.color], lty = 1)
    min.x = floor(min(xlim))
    max.x = ceiling(max(xlim))
    tick.at = seq(min.x, max.x, 1)
    tick.labels = tick.at
    axis(side = 1, at = tick.at, labels = tick.labels, las = 2, cex.axis = 0.5)
    dev.off()
        
    
    ####################
    #Plot
    ####################

    annot.filt = annot
    conc.col = 'mix1.conc'
    annot.filt = annot.filt[order(annot.filt[, conc.col]), ]
    ercc.id = annot.filt[, 'ERCC.ID']

    
    #plot each batch separately
    batch.cells = tapply(meta.mat[, 'id'], meta.mat[, strat.factor], unique)
    batches = names(batch.cells)
    for(j.batch in batches){

        #subset cells
        jbatch.cells = batch.cells[[j.batch]]
        ercc.log.filt = log2(ercc.rpkm[, jbatch.cells])

        #set filenames
        j.ercc.pdf = sub('ercc.pdf', paste(strat.factor, '_', j.batch, '.ercc.pdf', sep = ''), ercc.pdf)
        j.mean2cv.pdf = sub('mean2cv.pdf', paste(strat.factor, '_', j.batch, '.mean2cv.pdf', sep = ''), mean2cv.pdf)
        j.real.mean2cv.pdf = sub('real.mean2cv.pdf', paste(strat.factor, '_', j.batch, '.real.mean2cv.pdf', sep = ''), real.mean2cv.pdf)

        ########
        #RPKM vs conc
        ########
    
        #order by conc
        ercc.log.filt = ercc.log.filt[ercc.id, ]

        #set -inf to -10
        ercc.log.filt.noinf = t(as.matrix(apply(ercc.log.filt, 1, function(j.ercc){j.ercc[which(j.ercc == '-Inf')] = -10; return(j.ercc)})))
        dim(ercc.log.filt.noinf) = dim(ercc.log.filt)
        
        #plot
        pdf(file = j.ercc.pdf)
        plot.ercc(ercc.log.filt.noinf, log2(annot.filt[, conc.col]))
        dev.off()

    
        #####
        #CV^2 vs mean
        #####
        #variance for each transcript, or variance for each conc
        j.ercc = ercc.rpkm[, jbatch.cells]
        
        #get stats
        ts2stats = mean.sd.cv(j.ercc)
        pass.ts = rownames(ts2stats)[which(log10(ts2stats[, 'm']) > 0)]
        mean.log = log10(ts2stats[pass.ts, 'm'])
        cv2 = ts2stats[pass.ts, 'cv']^2
        
        #plot
        pdf(file = j.mean2cv.pdf)
        plot(mean.log, log10(cv2), xlab = 'log10(mean(RPKM))', ylab = 'log10(CV^2)', pch = 19, col = 'red')
        dev.off()

        #add real data
        rpkm = real.rpkm[, jbatch.cells]
        
        ts2stats = mean.sd.cv(rpkm)
        pass.ts = rownames(ts2stats)[which(log10(ts2stats[, 'm']) > 0)]
        real.mean.log = log10(ts2stats[pass.ts, 'm'])
        real.cv2 = ts2stats[pass.ts, 'cv']^2

        #plot
        pdf(file = j.real.mean2cv.pdf)
        smoothScatter(real.mean.log, log10(real.cv2), xlab = 'log10(mean(RPKM))', ylab = 'log10(CV^2)')
        points(mean.log, log10(cv2), pch = 19, col = 'red')
        dev.off()                
    }

    
    ###
    #Fraction of counts ERCC vs Real
    ###
    
    #calculate fractions
    sample2sum = colSums(real.counts)
    sample2sum.ercc = colSums(ercc.counts)
    ercc.frac = sample2sum.ercc / (sample2sum + sample2sum.ercc)    

    #barplot
    pdf(file = fraction.ercc.barplot.pdf, width = pdf.w, height = pdf.h)
    bp = barplot(ercc.frac, axes = FALSE, axisnames = FALSE, ylab = 'Fraction of ERCC reads', main = '', col = meta.mat[, strat.factor.color], border = NA)
    axis(2)
    axis(1, at = bp, labels = meta.mat[, 'id'], cex.axis = cex, las = 2)
#    abline(h = nreads.cutoff)
#    legend('topleft', legend = factor2color.map[, strat.factor], col = factor2color.map[, strat.factor.color], lty = 1)
    dev.off()

    #density
    pdf(file = fraction.ercc.dens.pdf)
    plot(density(ercc.frac), xlab = 'Fraction of ERCC reads', main= '', col = meta.mat[, strat.factor.color])
    dev.off()


    ###
    #Filter: Exclude ERCC transcripts with 0 added transcripts
    ###
    added.ercc.ids = annot[which(annot[, 'n.molecules'] != 0), 'ERCC.ID']
    ercc.counts = ercc.counts[added.ercc.ids, ]
    annot = annot[added.ercc.ids, ]
    print(length(added.ercc.ids))    

    
    ###
    #RPKM vs n.mol
    ###

    annot.filt = annot
    conc.col = 'n.molecules'
    annot.filt = annot.filt[order(annot.filt[, conc.col]), ]
    ercc.id = annot.filt[, 'ERCC.ID']

    batch.cells = tapply(meta.mat[, 'id'], meta.mat[, strat.factor], unique)
    batches = names(batch.cells)
    for(j.batch in batches){

        #subset cells
        jbatch.cells = batch.cells[[j.batch]]
        ercc.log.filt = log2(ercc.rpkm[, jbatch.cells])

        #set filenames
        j.ercc.pdf = sub('ercc.pdf', paste(strat.factor, '_', j.batch, '.nmol.ercc.pdf', sep = ''), ercc.pdf)
        j.mean2cv.pdf = sub('mean2cv.pdf', paste(strat.factor, '_', j.batch, '.mean2cv.pdf', sep = ''), mean2cv.pdf)
        j.real.mean2cv.pdf = sub('real.mean2cv.pdf', paste(strat.factor, '_', j.batch, '.real.mean2cv.pdf', sep = ''), real.mean2cv.pdf)
    
        #order by conc
        ercc.log.filt = ercc.log.filt[ercc.id, ]

        #set -inf to -10
        ercc.log.filt.noinf = t(as.matrix(apply(ercc.log.filt, 1, function(j.ercc){j.ercc[which(j.ercc == '-Inf')] = -10; return(j.ercc)})))
        dim(ercc.log.filt.noinf) = dim(ercc.log.filt)
        
        #plot
        pdf(file = j.ercc.pdf)
        plot.ercc(ercc.log.filt.noinf, log2(annot.filt[, conc.col]), xlab = 'log2(n molecules)')
        dev.off()
    }

    ###
    #linear fit
    ###
    
}
