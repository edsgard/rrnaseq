#Differential expression analysis


run.samseq <- function(meta.mat, counts, class.factor, factor.contrasts, alpha = 0.05){

    #get response
    data.samples = colnames(counts)
    group.samples = get.group.samples(meta.mat, class.factor, factor.contrasts, data.samples)
    
    c1.samples = group.samples[['c1.samples']]
    c2.samples = group.samples[['c2.samples']]
    annot.samples = c(c1.samples, c2.samples)
        
    if(length(annot.samples) > 0){

        cat(sprintf('n.samples in group "%s": %i\n', factor.contrasts[1], length(c1.samples)))
        cat(sprintf('n.samples in group "%s": %i\n', factor.contrasts[2], length(c2.samples)))

        #get response
        y = rep(1, length(annot.samples))
        names(y) = annot.samples
        y[c2.samples] = 2
        
        #subset on annotated samples present in counts
        j.counts = counts[, annot.samples]
    
        #Test
        sam.res = SAMseq(j.counts, y, resp.type = "Two class unpaired", fdr.output = alpha)

        #Get genes with FDR <= alpha (q-val < 100 * alpha)
        sig.res = get.samseq.sig(sam.res, alpha, fix.genes = TRUE, j.counts = j.counts)
        nrow(sig.res) #

        #Order by Score (test-stat) and log2.fc
        sig.res = sig.res[order(abs(sig.res[, 'Score(d)']), abs(sig.res[, 'log2.fc']), decreasing = TRUE), ]

        #Add rank column
        rank = 1:nrow(sig.res)
        sig.res = cbind(sig.res, rank, stringsAsFactors = FALSE)

    }else{
        sig.res = NA
    }
    return(sig.res)
}

get.group.samples <- function(meta.mat, class.factor, factor.contrasts, data.samples){
    
    c1.samples = rownames(meta.mat)[which(meta.mat[, class.factor] == factor.contrasts[1])]
    c2.samples = rownames(meta.mat)[which(meta.mat[, class.factor] == factor.contrasts[2])]

    #get annotated samples present in data
    c1.samples = intersect(c1.samples, data.samples)
    c2.samples = intersect(c2.samples, data.samples)

    return(list(c1.samples = c1.samples, c2.samples = c2.samples))
}

get.fc.stats <- function(group.samples, rpkm, eps = 1e-10){
    
    c1.samples = group.samples[[1]]
    c1.rpkm = rpkm[, c1.samples]
    c1.stats = get.row.stats(c1.rpkm)
    colnames(c1.stats) = paste('c1', colnames(c1.stats), sep = '.')
    
    c2.samples = group.samples[[2]]
    c2.rpkm = rpkm[, c2.samples]
    c2.stats = get.row.stats(c2.rpkm)
    colnames(c2.stats) = paste('c2', colnames(c2.stats), sep = '.')
    
    fc = (c2.stats[, 'c2.mean'] + eps) / (c1.stats[, 'c1.mean'] + eps)
    log2.fc = log2(fc)

    fc.stats = cbind(c1.stats, c2.stats, fc, log2.fc)

    return(fc.stats)
}

get.row.stats <- function(data.mat){
    median = apply(data.mat, 1, median)
    mean = apply(data.mat, 1, mean)
    sd = apply(data.mat, 1, sd)
    cv = sd / mean

    row.stats = cbind(median, mean, sd, cv)
    return(row.stats)
}

get.samseq.sig <- function(sam.res, alpha, fix.genes = FALSE, j.counts = NA){

    #Get up-reg genes
    genes.up = sam.res[['siggenes.table']][['genes.up']]

    #Fix data types
    genes.up = as.data.frame(genes.up, stringsAsFactors = FALSE)
    num.cols = c('Score(d)', 'Fold Change', 'q-value(%)')
    genes.up[, num.cols] = lapply(genes.up[, num.cols], as.numeric)

    #get sig
    up.sig = genes.up[which(genes.up[, 'q-value(%)'] <= alpha * 100), ]
    n.sig = nrow(up.sig)
    direction = rep('up', n.sig)
    up.sig = cbind(up.sig, direction, stringsAsFactors = FALSE)

    #Get down-reg genes
    genes.lo = sam.res[['siggenes.table']][['genes.lo']]

    #Fix data types
    genes.lo = as.data.frame(genes.lo, stringsAsFactors = FALSE)
    num.cols = c('Score(d)', 'Fold Change', 'q-value(%)')
    genes.lo[, num.cols] = lapply(genes.lo[, num.cols], as.numeric)

    #get sig
    lo.sig = genes.lo[which(genes.lo[, 'q-value(%)'] <= alpha * 100), ]
    n.sig = nrow(lo.sig)
    direction = rep('down', n.sig)
    lo.sig = cbind(lo.sig, direction, stringsAsFactors = FALSE)

    #bind up- and down-reg
    sig.res = rbind(up.sig, lo.sig)

    #add log2(fold-change) column
    log2.fc = log2(sig.res[, 'Fold Change'])
    sig.res = cbind(sig.res, log2.fc, stringsAsFactors = FALSE)
    
    #Fix gene names
    if(fix.genes){
        gene.ind = as.numeric(sig.res[, 'Gene Name'])
        gene = rownames(j.counts)[gene.ind]
        sig.res[, 'Gene Name'] = gene
    }

    rownames(sig.res) = sig.res[, 'Gene Name']
    
    #sort by score
    sig.res = sig.res[order(abs(sig.res[, 'Score(d)']), decreasing = TRUE), ]
    
    return(sig.res)
}
