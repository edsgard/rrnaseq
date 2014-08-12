
plot.cv2vsmean <- function(real.stats, ercc.stats, tech.fit.res, sel.genes, sel.highlight = TRUE){
#plot.cv2vsmean(real.stats, ercc.stats, ercc.fit.res, sel.genes, sel.highlight = sel.highlight)
    
    #get gene summary stats
    gene.means = real.stats[, 'mean']
    gene.cv2 = real.stats[, 'cv2']
    min.biol.cv2 = unique(real.stats[, 'min.biol.cv2'])
    
    ercc.gene.means = ercc.stats[, 'mean']
    ercc.gene.cv2 = ercc.stats[, 'cv2']

    #get tech fit res
    a0 = tech.fit.res['a0']
    a1.tilde = tech.fit.res['a1.tilde']
    psia1theta = tech.fit.res['psia1theta']
    
    #Limits
    xlim = range(log10(gene.means))
    xlim[1] = floor(xlim[1])
    xlim[2] = ceiling(xlim[2])
    
    ylim = range(log10(gene.cv2))
    ylim[1] = floor(ylim[1])
    ylim[2] = ceiling(ylim[2])
    
    #labels
    xlab = "average normalized read count"
    ylab = "squared coefficient of variation (CV^2)"

    #color map
    n.genes = length(gene.means)
    gene2color = rep('blue', n.genes)
    names(gene2color) = names(gene.means)
    if(sel.highlight){
        gene2color[sel.genes] = '#C0007090'
    }
    
    #Plot the real data
    if(sel.highlight){
        plot(log10(gene.means), log10(gene.cv2), pch=20, cex=.2, col = gene2color, xlab = xlab, ylab = ylab, xaxt = 'n', yaxt = 'n')
    }else{
        smoothScatter(log10(gene.means), log10(gene.cv2), pch=20, cex=.2, col = gene2color, xlab = xlab, ylab = ylab, xaxt = 'n', yaxt = 'n')        
    }    
    x = xlim[1]:xlim[2]
    y = ylim[1]:ylim[2]
    axis(1, x, as.character(10^x))
    axis(2, y, as.character(10^y), las = 2)
        
    #Add the technical noise fit, as before
    xg = 10^seq(xlim[1], xlim[2], length.out=1000)
    lines(log10(xg), log10((a1.tilde)/xg + a0), col="#FF000080", lwd=3)
    
    #Add a curve showing the expectation for the chosen biological CV^2 threshold
    lines(log10(xg), log10(psia1theta/xg + a0 + min.biol.cv2), lty="dashed", col="#C0007090", lwd=3 )
    
    #Add the normalised ERCC points
    points(log10(ercc.gene.means), log10(ercc.gene.cv2), pch=20, cex=1, col="black")

}    

get.real.stats <- function(expr, tech.fit.res, ercc.size.factors, real.size.factors, min.biol.cv2){

    #basic summary stats
    gene.means = rowMeans(expr)
    gene.vars = rowVars(expr)
    gene.cv2 = gene.vars / gene.means^2

    #tech fit res
    a0 = tech.fit.res['a0']
    a1.tilde = tech.fit.res['a1.tilde']
    psia1theta = tech.fit.res['psia1theta']
    
    #other test-statistica
    m = ncol(expr)
    cv2 = a0 + min.biol.cv2 + a0 * min.biol.cv2
    testDenom = (gene.means * psia1theta + gene.means^2 * cv2) / ( 1 + cv2/m )
    gene.var.quants = gene.vars * (m-1) / testDenom
    
    #get p-value from a chi-sq dist
    p = 1 - pchisq(gene.var.quants, m-1)

    #mult-test correction
    p.adj = p.adjust(p, "BH")

    #bind res
    test.res = cbind(gene.means, gene.vars, gene.cv2, min.biol.cv2, gene.var.quants, p, p.adj)
    colnames(test.res)[1:3] = c('mean', 'var', 'cv2')
    
    #order by test var
    test.res = test.res[order(test.res[, 'gene.var.quants'], decreasing = TRUE), ]
 
    return(test.res)
}

plot.ercc.fit <- function(ercc.stats, fit.res, df, quant.min = 0.025, add.quant = TRUE, cex = 1){

    gene.means = ercc.stats[, 'mean']
    gene.cv2 = ercc.stats[, 'cv2']

    a0 = fit.res['a0']
    a1.tilde = fit.res['a1.tilde']
    
    #xlim
    xlim = range(log10(gene.means))
    xlim[1] = floor(xlim[1])
    xlim[2] = ceiling(xlim[2])

    #ylim
    ylim = range(log10(gene.cv2))
    ylim[1] = floor(ylim[1])
    ylim[2] = ceiling(ylim[2])
    
    #labels
    xlab = "average normalized read count"
    ylab = "squared coefficient of variation (CV^2)"
        
    #Plot the observed data
    plot(log10(gene.means), log10(gene.cv2), pch=20, cex = cex, col = 'blue', xlab = xlab, ylab = ylab, xaxt = 'n', yaxt = 'n')
    x = xlim[1]:xlim[2]
    y = ylim[1]:ylim[2]
    axis(1, x, as.character(10^x))
    axis(2, y, as.character(10^y), las = 2)
        
    #Add the technical noise fit, as before
    xg = 10^seq(xlim[1], xlim[2], length.out=1000)
    lines(log10(xg), log10((a1.tilde)/xg + a0), col="#FF000080", lwd=3)
                
    #Plot quantile lines around the fit
    if(add.quant){
        quant.max = 1 - quant.min
        lines(log10(xg), log10(( (a1.tilde)/xg + a0  ) * qchisq( quant.max, df ) / df), col="#FF000080", lwd=2, lty="dashed")
        lines(log10(xg), log10(( (a1.tilde)/xg + a0  ) * qchisq( quant.min, df ) / df), col="#FF000080", lwd=2, lty="dashed")
    }
}    

ercc.fit <- function(ercc.stats, ercc.size.factors, real.size.factors, pass.genes){
                
    #Fit
    fit = glmgam.fit(cbind(a0 = 1, a1.tilde = 1 / ercc.stats[pass.genes, 'mean']), ercc.stats[pass.genes, 'cv2'])

    #get coefs
    xi = mean(1 / ercc.size.factors)
    a0 = unname(fit$coefficients["a0"])
    a1.tilde = fit$coefficients["a1.tilde"]
    a1 = unname(a1.tilde - xi)

    #psi + a1*theta
    pass.samples = names(ercc.size.factors)
    psia1theta = mean(1 / real.size.factors[pass.samples]) + a1 * mean(ercc.size.factors / real.size.factors[pass.samples])

    #Explained variance
    residual.var = var(log(fitted.values(fit)) - log(ercc.stats[pass.genes, 'cv2']) )
    total.var = var(log(ercc.stats[pass.genes, 'cv2']))
    expl.var = 1 - residual.var / total.var

    #bind
    fit.res = c(xi, a0, a1.tilde, a1, psia1theta, residual.var, total.var, expl.var)
    names(fit.res) = as.character(expression(xi, a0, a1.tilde, a1, psia1theta, residual.var, total.var, expl.var))

    return(fit.res)
}

get.ercc.stats <- function(ercc){
    
    mean = rowMeans(ercc)
    var = rowVars(ercc)
    cv2 = var / mean^2

    #bind
    ercc.stats = cbind(mean, var, cv2)

    return(ercc.stats)
}

ercc.filter <- function(ercc.raw, n.ercc.min, min.count, min.prop){
    
    #Rm samples with no ERCC in the rpkm annot (they didn't have ERCC spiked in to begin with either)
    ercc.ind = which(apply(ercc.raw, 2, function(j.ercc){length(which(!is.na(j.ercc))) > 0}))
    pass.samples = names(ercc.ind) 
    ercc.raw = ercc.raw[, pass.samples]
    ncol(ercc.raw) #456

    #Rm samples with reads mapped to less than n ERCC transcripts
    ercc.ind = which(apply(ercc.raw, 2, function(j.ercc){length(which(j.ercc != 0)) > n.ercc.min}))
    pass.samples = names(ercc.ind)
    length(pass.samples) #446
    ercc.raw = ercc.raw[, pass.samples]
        
    #For each cell, proportion of expressed transcripts > min.rpkm should be greater than min.prop
    ercc.ind = which(apply(ercc.raw, 2, function(j.ercc){j.ercc = j.ercc[which(j.ercc != 0)]; prop.expr = length(which(j.ercc > min.count)) / length(j.ercc); return(prop.expr >= min.prop)}))
    pass.samples = names(ercc.ind)    
    length(pass.samples) #352

    #apply filter
    ercc.raw = ercc.raw[, pass.samples]    

    ###
    #Filter ercc transcripts that are all zero
    ###
    ercc.ind = which(apply(ercc.raw, 1, function(j.ercc){length(which(j.ercc > min.count)) >= n.ercc.min}))
    pass.genes = names(ercc.ind)
    ercc.raw = ercc.raw[pass.genes, ]

    return(ercc.raw)
}

var.genes.brennecke <- function(count.file, ercc.count.file, meta.file, out.dir, params){

    
    ###
    #Libs
    ###
    library('DESeq2')
    library('genefilter')
    library('statmod')

    
    ###
    #Set output files
    ###
    
    #PDFs
    ercc.fit.pdf = file.path(out.dir, 'pdf', 'ercc.fit.pdf')
    cv2vsmean.pdf = file.path(out.dir, 'pdf', 'cv2vsmean.pdf')

    #RDS out
    stats.list.rds = file.path(out.dir, 'rds', 'stats.list.rds')
    stats.mat.rds = file.path(out.dir, 'rds', 'stats.mat.rds')
    stats.mat.tab = file.path(out.dir, 'rds', 'stats.mat.tab')


    #Create output dirs
    dir.create(dirname(ercc.fit.pdf), recursive = TRUE, showWarnings = FALSE)
    dir.create(dirname(stats.list.rds), recursive = TRUE, showWarnings = FALSE)

    
    ###
    #Read data
    ###
    count = readRDS(count.file)
    ercc.raw = readRDS(ercc.count.file)
    meta.mat = readRDS(meta.file)

    
    ###
    #Filter samples such that all agree with samples in count
    ###

    pass.samples = colnames(count)
    count = count[, pass.samples]
    ercc.raw = ercc.raw[, pass.samples]
    meta.mat = meta.mat[pass.samples, ]

    dim(ercc.raw) #456

    
    ###
    #Filter ERCC samples on presence of spike-in
    ###

    #Params
    n.ercc.min = params[['n.ercc.min']]
    min.count = params[['min.count']]
    min.prop = params[['min.prop']]

    #Filter
    ercc.filt = ercc.filter(ercc.raw, n.ercc.min, min.count, min.prop)
    dim(ercc.filt) #84, 352    

    
    ###
    #Normalize counts
    ###
    #Real data
    real.size.factors = estimateSizeFactorsForMatrix(count)
    expr = t( t(count) / real.size.factors )

    #ERCC
    ercc.size.factors = estimateSizeFactorsForMatrix(ercc.filt)
    ercc = t( t(ercc.filt) / ercc.size.factors )

    
    #######
    #Estimate tech noise
    #######

    #Params
    min.cv2 = params[['min.cv2']]
    quant.cutoff = params[['quant.cutoff']]
    
    #Get stats
    ercc.stats = get.ercc.stats(ercc)    
    
    #Filter genes on min mean expression
    min.mean = unname( quantile( ercc.stats[which(ercc.stats[, 'cv2'] > min.cv2 ), 'mean'], quant.cutoff ) )
    pass.genes = rownames(ercc.stats)[which(ercc.stats[, 'mean'] >= min.mean)]
    length(pass.genes) #25

    #Fit ERCC
    ercc.fit.res = ercc.fit(ercc.stats, ercc.size.factors, real.size.factors, pass.genes)
    print(ercc.fit.res['expl.var'])
    #all: 0.91
            
    
    
    ######
    #Test for high var
    ######

    #Params
    strat.factor = params[['strat.factor']]
    min.biol.cv2 = params[['min.biol.cv2']]
    min.expr = params[['min.expr']]
    min.cells.expr = params[['min.cells.expr']]
        
    #Get stats and test
    real.stats = get.real.stats(expr, ercc.fit.res, ercc.size.factors, real.size.factors, min.biol.cv2)

    #Split cells by strat.factor
    meta2stage.list = split(meta.mat, meta.mat[, strat.factor])
    stages = names(meta2stage.list)
    
    n.stages = length(stages)
    real.stats.list = list()
    length(real.stats.list) = n.stages
    names(real.stats.list) = stages
    for(stage in stages){

        print(stage)
        
        #subset cells
        samples.subset = meta2stage.list[[stage]][, 'id']
        expr.subset = expr[ , samples.subset]
        real.size.factors.subset = real.size.factors[samples.subset]

        #Rm genes with expr in <= min.cells.expr
        gene2nsamples.expr = apply(expr.subset, 1, function(jgene.rpkm, rpkm.cutoff){length(which(jgene.rpkm > rpkm.cutoff));}, rpkm.cutoff = min.expr)
        rm.genes = names(gene2nsamples.expr)[which(gene2nsamples.expr < min.cells.expr)]
        expr.subset = expr.subset[setdiff(rownames(expr.subset), rm.genes), ]
        
        #get stats
        j.real.stats = get.real.stats(expr.subset, ercc.fit.res, ercc.size.factors, real.size.factors.subset, min.biol.cv2)

        #set removed genes to NA
        rm.genes.mat = matrix(NA, nrow = length(rm.genes), ncol = ncol(j.real.stats), dimnames = list(rm.genes, colnames(j.real.stats)))
        
        j.real.stats = rbind(j.real.stats, rm.genes.mat)
        

        #save
        real.stats.list[[stage]] = j.real.stats
    }

    #add stats using all cells
    real.stats.list[['nostrat']] = real.stats
    

    #########
    #Dump
    #########
    
    #Make a table with all genes merging the real.stats.list    
    stages = names(real.stats.list)
    all.genes = unique(unlist(lapply(real.stats.list, rownames)))
    all.stats.mat = as.data.frame(matrix(all.genes, ncol = 1, nrow = length(all.genes)))    
    rownames(all.stats.mat) = all.genes
    colnames(all.stats.mat) = 'gene'
    all.stats.mat[, 'gene'] = all.genes
        
    for(stage in stages){

        #add stage to colname
        j.real.stats = real.stats.list[[stage]]
        colnames(j.real.stats) = paste(stage, colnames(j.real.stats), sep = '.')
        
        all.stats.mat = merge(all.stats.mat, j.real.stats, by.x = 'gene', by.y = 'row.names', stringsAsFactors = FALSE)
    }

    #Dump
    saveRDS(real.stats.list, file = stats.list.rds)
    saveRDS(all.stats.mat, file = stats.mat.rds)    
    write.table(all.stats.mat, quote = FALSE, row.names = FALSE, sep = '\t', file = stats.mat.tab)

    
    ############
    #Plot
    ############    
    
    ###
    #Plot fit of ERCC
    ###
    
    quant.min = params[['quant.min']]
    
    n.samples = ncol(expr)
    df = n.samples - 1
    pdf(ercc.fit.pdf)
    plot.ercc.fit(ercc.stats, ercc.fit.res, df, quant.min = quant.min, add.quant = FALSE)
    dev.off()

    
    ###
    #Plot real data
    ###

    #Params
    alpha = params[['alpha']]
    n.top = params[['n.top']]
    
    #Plot every stage
    stages = names(real.stats.list)
    for(stage in stages){
        
        real.stats = real.stats.list[[stage]]

        #Filter NAs
        real.stats = real.stats[which(!is.na(real.stats[, 'cv2'])), ]
        
        #Select subset of sig/top genes
        sig.genes = rownames(real.stats)[which(real.stats[, 'p.adj'] <= alpha)]
        print(length(sig.genes))
        #length(which(p == 0))
        
        sel.genes = rownames(real.stats)[1:n.top]

        #set dir
        j.dir = file.path(dirname(cv2vsmean.pdf), 'stage', stage)
        dir.create(j.dir, recursive = TRUE, showWarnings = FALSE)
        stage.cv2vsmean.pdf = file.path(j.dir, basename(cv2vsmean.pdf))
        
        #plot
        sel.highlight = FALSE
        j.cv2vsmean.pdf = sub('cv2vsmean', paste('cv2vsmean.selhighlight_', sel.highlight, sep = ''), stage.cv2vsmean.pdf)        
        pdf(j.cv2vsmean.pdf)
        plot.cv2vsmean(real.stats, ercc.stats, ercc.fit.res, sel.genes, sel.highlight = sel.highlight)
        dev.off()

        sel.highlight = TRUE
        j.cv2vsmean.pdf = sub('cv2vsmean', paste('cv2vsmean.selhighlight_', sel.highlight, '.top_', n.top, sep = ''), stage.cv2vsmean.pdf)
        pdf(j.cv2vsmean.pdf)
        plot.cv2vsmean(real.stats, ercc.stats, ercc.fit.res, sel.genes, sel.highlight = sel.highlight)
        dev.off()

        sel.highlight = TRUE
        j.cv2vsmean.pdf = sub('cv2vsmean', paste('cv2vsmean.selhighlight_', sel.highlight, '.alpha_', alpha, sep = ''), stage.cv2vsmean.pdf)
        pdf(j.cv2vsmean.pdf)
        plot.cv2vsmean(real.stats, ercc.stats, ercc.fit.res, sig.genes, sel.highlight = sel.highlight)
        dev.off()
    }

            
}
