#########
#Change working dir to output dir
#########
cd /Volumes/Data/cloud/gdrive/work/rspd/code/my/git/rrnaseq/test/rqc/data


#########
#Get data
#########
#Get sample meta-information
get_meta -i /Volumes/Data/cloud/gdrive/work/rspd/code/my/git/rrnaseq/test/rqc/data/mapstats.tab -o meta.rds

#Get expression data-structures
get_expr -i /Volumes/Data/cloud/gdrive/work/rspd/code/my/git/rrnaseq/test/rpkmforgenes/refseq_rpkms.tab -o /Volumes/Data/cloud/gdrive/work/rspd/code/my/git/rrnaseq/test/rqc/data


#########
#Plot basic expression QC-metrics
#########
#Plot mapping stats
mapstats -m meta.rds -c counts.rds -o /Volumes/Data/cloud/gdrive/work/rspd/code/my/git/rrnaseq/test/rqc/pdf/seq_mapping -q qc.rds

#Plot sample expression histogram
sample_expr_dhist -m meta.rds -r rpkm.rds -o /Volumes/Data/cloud/gdrive/work/rspd/code/my/git/rrnaseq/test/rqc/pdf/basic/sample.expr.dens.pdf

#Filter genes based on the number of samples they are expressed in and plot histogram of it
gene_filter -m meta.rds -r rpkm.rds -o rpkm.postqc.rds -d /Volumes/Data/cloud/gdrive/work/rspd/code/my/git/rrnaseq/test/rqc/pdf/basic/gene2nsamples.pdf

#Plot number of expressed genes per sample
sample2ngenes_expr -m meta.rds -r rpkm.rds -o /Volumes/Data/cloud/gdrive/work/rspd/code/my/git/rrnaseq/test/rqc/pdf/basic/sample2ngenes.dens.pdf -q qc.rds

#Sample distance heatmap
sampledist_heatmap -m meta.rds -r rpkm.rds -o /Volumes/Data/cloud/gdrive/work/rspd/code/my/git/rrnaseq/test/rqc/pdf/sampledist/sample.heatmap.pdf -e sampledist.cor.res.rds -n 10 -k 4

#Sample distance boxplot
sampledist_boxplot -m meta.rds -r rpkm.rds -o /Volumes/Data/cloud/gdrive/work/rspd/code/my/git/rrnaseq/test/rqc/pdf/sampledist/sample.cor.pdf -q qc.rds -e sampledist.cor.res.rds

#Sample hierarchical clustering
#no boostrap
sample_hclust -m meta.rds -r rpkm.rds -o /Volumes/Data/cloud/gdrive/work/rspd/code/my/git/rrnaseq/test/rqc/pdf/sampledist/sample.hclust.pdf -b 0 -e sampledist.cor.res.rds
#with bootstrap
#sample_hclust -m meta.rds -r rpkm.rds -o /Volumes/Data/cloud/gdrive/work/rspd/code/my/git/rrnaseq/test/rqc/pdf/sampledist/sample.hclust.b_10.pdf -b 10

#PCA
pca -m meta.rds -r rpkm.rds -o /Volumes/Data/cloud/gdrive/work/rspd/code/my/git/rrnaseq/test/rqc/pdf/pca/pca.pdf -c 1,2,3
pca -m meta.rds -r rpkm.rds -o /Volumes/Data/cloud/gdrive/work/rspd/code/my/git/rrnaseq/test/rqc/pdf/pca/pca.pdf -c 1,2
pca -m meta.rds -r rpkm.rds -o /Volumes/Data/cloud/gdrive/work/rspd/code/my/git/rrnaseq/test/rqc/pdf/pca/pca.pdf -c 1 -f -q qc.rds

#Scatter plot
scatter -d rpkm.rds -o /Volumes/Data/cloud/gdrive/work/rspd/code/my/git/rrnaseq/test/rqc/pdf/scatter/sample.scatter.pdf -s 1,2,3


#########
#Apply sample filters (and gene filter on the new set of samples)
#########

#Filter samples based on qc-measures
sample_filter -d rpkm.rds -q qc.rds -o rpkm.postqc.rds -f n.genes.expr,spearman.cor
sample_filter -d counts.rds -q qc.rds -o counts.postqc.rds -f n.genes.expr,spearman.cor

#Filter genes using the new (filtered) set of samples
gene_filter -m meta.rds -r rpkm.postqc.rds -o rpkm.postqc.rds -p FALSE
gene_filter -m meta.rds -r counts.postqc.rds -o counts.postqc.rds -p FALSE


#########
#Plot ERCC stats
#########
ercc -m meta.rds -r rpkm.rds -c counts.rds -e ercc.rpkm.rds -d ercc.counts.rds -o /Volumes/Data/cloud/gdrive/work/rspd/code/my/git/rrnaseq/test/rqc/pdf/ercc/fracadjusted_false
ercc -m meta.rds -r rpkm.rds -c counts.rds -e ercc.rpkm.rds -d ercc.counts.rds -f -o /Volumes/Data/cloud/gdrive/work/rspd/code/my/git/rrnaseq/test/rqc/pdf/ercc/fracadjusted_true


#########
#Differential expression analysis
#########
#Find genes with higher variance than the technical variance (ERCC)
var_genes_brennecke -m meta.rds -c counts.postqc.rds -e ercc.counts.rds -o /Volumes/Data/cloud/gdrive/work/rspd/code/my/git/rrnaseq/test/rqc/diffexp/brennecke

#Heatmap of most variable genes found by var_genes_brennecke
head -n100 /Volumes/Data/cloud/gdrive/work/rspd/code/my/git/rrnaseq/test/rqc/diffexp/brennecke/rds/genes.ntop_500.txt >/Volumes/Data/cloud/gdrive/work/rspd/code/my/git/rrnaseq/test/rqc/diffexp/brennecke/rds/genes.ntop_100.txt
heatmap -i rpkm.rds -r /Volumes/Data/cloud/gdrive/work/rspd/code/my/git/rrnaseq/test/rqc/diffexp/brennecke/rds/genes.ntop_100.txt -o /Volumes/Data/cloud/gdrive/work/rspd/code/my/git/rrnaseq/test/rqc/diffexp/brennecke/pdf/topgenes.heatmap.pdf
