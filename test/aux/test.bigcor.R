
#Libs
cloud.dir = '/Volumes/Data/cloud/gdrive'
lib.dir = file.path(cloud.dir, 'work/rspd/code/my/git/rrnaseq/lib')
lib.file = file.path(lib.dir, 'basic.lib.R')
source(lib.file)

cor.meth = 'spearman'
nblocks = 3
ncore = 2

#Execute
#main(cor.meth, nblocks, ncore)

main <- function(cor.meth, nblocks, ncore){
    
    ###
    #Test with one MAT
    ###
    print('One mat')
    
    #make data
    MAT = matrix(rnorm(8 * 10), nrow = 10)

    #big cor
    big.res = bigcor.par(MAT, nblocks = nblocks, ncore = ncore, method = cor.meth)

    #original cor
    orig.res = cor(MAT, method = cor.meth)

    #diff
    res.diff = big.res - orig.res
    summary(as.numeric(res.diff))
    max.diff = max(res.diff)
    if(max.diff < 1e-6){
        print('PASSED')
    }


    ###
    #Test with two identical MATs
    ###
    print('Two identical mat')
    
    #make data
    MAT = matrix(rnorm(8 * 10), nrow = 10)

    #big cor
    big.res = bigcor.par(MAT, nblocks = nblocks, ncore = ncore, yMAT = MAT)

    #original cor
    orig.res = cor(x = MAT, y = MAT)

    #diff
    res.diff = big.res - orig.res
    summary(as.numeric(res.diff))
    max.diff = max(res.diff)
    if(max.diff < 1e-6){
        print('PASSED')
    }

    
    ###
    #Test with two different MATs
    ###
    print('Two different mat')
    
    #make data
    MAT = matrix(rnorm(8 * 10), nrow = 10)
    y = matrix(rnorm(8 * 10), nrow = 10)
    
    #big cor
    big.res = bigcor.par(MAT, nblocks = nblocks, ncore = ncore, yMAT = y)

    #original cor
    orig.res = cor(x = MAT, y = y)

    #diff
    res.diff = big.res - orig.res
    summary(as.numeric(res.diff))
    max.diff = max(res.diff)
    if(max.diff < 1e-6){
        print('PASSED')
    }
    
    
}

