

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

plot.pairs.pca <- function(pca.basis, e.var, meta.mat = NA, factor.color = 'stage', point.text.cex = 0.7){

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
    point.panel <- function(x, y){points(x, y, col = cell.color, pch = 16)}
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
