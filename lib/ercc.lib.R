
mean.sd.cv <- function(ercc){

    #sd, mean
    ts2sd = apply(ercc, 1, sd, na.rm = TRUE)
    ts2mean = apply(ercc, 1, mean, na.rm = TRUE)
    
    #filter NA
    pass.ts = names(ts2sd)[which(!is.na(ts2sd))]
    ts2sd = ts2sd[pass.ts]
    ts2mean = ts2mean[pass.ts]

    #cv
    cv = ts2sd / ts2mean

    ts2stats = cbind(ts2mean, ts2sd, cv)
    colnames(ts2stats) = c('m', 'sd', 'cv')
    rownames(ts2stats) = pass.ts
    
    return(ts2stats)
}

plot.ercc <- function(ercc, control.conc, xlab = 'log2(mix1 conc)', ylab = 'log2(RPKM)'){
    
    ylim = range(ercc)
    n.samples = ncol(ercc)
    j.sample = 1
    plot(control.conc, ercc[, j.sample], col = j.sample, ylim = ylim, type = 'p', ylab = ylab, xlab = xlab)
    for(j.sample in 2:n.samples){
        points(control.conc, ercc[, j.sample], col = j.sample, ylim = ylim)
    }
}
