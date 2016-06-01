

 test.mat = matrix(rnorm(100), ncol = 10, nrow = 10, byrow = TRUE)

#Libs
cloud.dir = '/Volumes/Data/cloud/gdrive'
lib.dir = file.path(cloud.dir, 'work/rspd/code/my/git/rrnaseq/lib')
lib.file = file.path(lib.dir, 'pvclust.R')
source(lib.file)
lib.file = file.path(lib.dir, 'pvclust-internal.R')
source(lib.file)
lib.file = file.path(lib.dir, 'basic.lib.R')
source(lib.file)

main <- function(){

    #params
    n.boot = 2
    cor.dist.fcn = spear.cor.col.dist
    mc.cores = 2
    
    #orig
    orig.res = pvclust(test.mat, method.dist = cor.dist.fcn, nboot = n.boot)

    #par
    par.res = pvclust.par(test.mat, method.dist = cor.dist.fcn, nboot = 2, mc.cores = mc.cores)    

    
    
}
