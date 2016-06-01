#########
#Change working dir to output dir
#########
cd /mnt/kauffman/edsgard/cloud/btsync/work/rspd/code/my/git/rrnaseq/test/rqc/data


#########
#Get data
#########
#Get sample meta-information
get_meta -i /mnt/kauffman/edsgard/cloud/btsync/work/rspd/code/my/git/rrnaseq/test/rqc/data/mapstats.tab -o meta.rds

#Get expression data-structures
get_expr -i /mnt/kauffman/edsgard/cloud/btsync/work/rspd/code/my/git/rrnaseq/test/rpkmforgenes/refseq_rpkms.tab -o /mnt/kauffman/edsgard/cloud/btsync/work/rspd/code/my/git/rrnaseq/test/rqc/data


#########
#Plot basic expression QC-metrics
#########
#Plot mapping stats
mapstats -m meta.rds -c counts.rds -o /mnt/kauffman/edsgard/cloud/btsync/work/rspd/code/my/git/rrnaseq/test/rqc/pdf/seq_mapping -q qc.rds

#Plot sample expression histogram
expr_dhist -m meta.rds -r rpkm.rds -o /mnt/kauffman/edsgard/cloud/btsync/work/rspd/code/my/git/rrnaseq/test/rqc/pdf/basic/sample.expr.dens.pdf

#Filter genes based on the number of samples they are expressed in and plot histogram of it
gene_filter -m meta.rds -r rpkm.rds -o rpkm.postqc.rds -d /mnt/kauffman/edsgard/cloud/btsync/work/rspd/code/my/git/rrnaseq/test/rqc/pdf/basic/gene2nsamples.pdf

#Plot number of expressed genes per sample
sample2ngenes_expr -m meta.rds -r rpkm.rds -o /mnt/kauffman/edsgard/cloud/btsync/work/rspd/code/my/git/rrnaseq/test/rqc/pdf/basic/sample2ngenes.dens.pdf -q qc.rds

#Sample distance heatmap
sampledist_heatmap -m meta.rds -r rpkm.rds -o /mnt/kauffman/edsgard/cloud/btsync/work/rspd/code/my/git/rrnaseq/test/rqc/pdf/sampledist/sample.heatmap.pdf -e sampledist.cor.res.rds -n 10 -k 4
sampledist_heatmap -m meta.rds -r rpkm.rds -o /mnt/kauffman/edsgard/cloud/btsync/work/rspd/code/my/git/rrnaseq/test/rqc/pdf/sampledist/sample.heatmap.strat_ngenesexpr.pdf -e sampledist.cor.res.rds -s n.genes.expr

#Sample distance boxplot
sampledist_boxplot -m meta.rds -r rpkm.rds -o /mnt/kauffman/edsgard/cloud/btsync/work/rspd/code/my/git/rrnaseq/test/rqc/pdf/sampledist/sample.cor.pdf -q qc.rds -e sampledist.cor.res.rds

#Sample hierarchical clustering
#no boostrap
sample_hclust -m meta.rds -r rpkm.rds -o /mnt/kauffman/edsgard/cloud/btsync/work/rspd/code/my/git/rrnaseq/test/rqc/pdf/sampledist/sample.hclust.pdf -b 0 -e sampledist.cor.res.rds
#with bootstrap
#sample_hclust -m meta.rds -r rpkm.rds -o /mnt/kauffman/edsgard/cloud/btsync/work/rspd/code/my/git/rrnaseq/test/rqc/pdf/sampledist/sample.hclust.b_10.pdf -b 10

#PCA
pca -m meta.rds -r rpkm.rds -o /mnt/kauffman/edsgard/cloud/btsync/work/rspd/code/my/git/rrnaseq/test/rqc/pdf/pca/pca.pdf -c 1,2,3
pca -m meta.rds -r rpkm.rds -o /mnt/kauffman/edsgard/cloud/btsync/work/rspd/code/my/git/rrnaseq/test/rqc/pdf/pca/pca.strat_ngenesexpr.pdf -c 1,2,3 -s n.genes.expr
pca -m meta.rds -r rpkm.rds -o /mnt/kauffman/edsgard/cloud/btsync/work/rspd/code/my/git/rrnaseq/test/rqc/pdf/pca/pca.pdf -c 1,2
pca -m meta.rds -r rpkm.rds -o /mnt/kauffman/edsgard/cloud/btsync/work/rspd/code/my/git/rrnaseq/test/rqc/pdf/pca/pca.pdf -c 1 -f -q qc.rds

#Scatter plot
scatter -d rpkm.rds -o /mnt/kauffman/edsgard/cloud/btsync/work/rspd/code/my/git/rrnaseq/test/rqc/pdf/scatter/sample.scatter.pdf -s 1,2,3


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
#Save as tab-format
#########
rds2tab -i meta.rds -o meta.tab
rds2tab -i rpkm.postqc.rds -o rpkm.postqc.tab
rds2tab -i counts.postqc.rds -o counts.postqc.tab


#########
#Plot ERCC stats
#########
ercc -m meta.rds -r rpkm.rds -c counts.rds -e ercc.rpkm.rds -d ercc.counts.rds -o /mnt/kauffman/edsgard/cloud/btsync/work/rspd/code/my/git/rrnaseq/test/rqc/pdf/ercc/fracadjusted_false
ercc -m meta.rds -r rpkm.rds -c counts.rds -e ercc.rpkm.rds -d ercc.counts.rds -f -o /mnt/kauffman/edsgard/cloud/btsync/work/rspd/code/my/git/rrnaseq/test/rqc/pdf/ercc/fracadjusted_true


#########
#Differential expression analysis
#########
#Find genes with higher variance than the technical variance (ERCC)
#NB: Suggested to set -w to at least 1, here set to 0 since this test dataset only has 3 cells.
var_genes_brennecke -w 0 -m meta.rds -c counts.postqc.rds -e ercc.counts.rds -o /mnt/kauffman/edsgard/cloud/btsync/work/rspd/code/my/git/rrnaseq/test/rqc/diffexp/brennecke

#Subset on genes that are the most variable across all cells
ntop=200
head -${ntop} /mnt/kauffman/edsgard/cloud/btsync/work/rspd/code/my/git/rrnaseq/test/rqc/diffexp/brennecke/rds/genes.ranked.txt >/mnt/kauffman/edsgard/cloud/btsync/work/rspd/code/my/git/rrnaseq/test/rqc/diffexp/brennecke/rds/genes.ntop_${ntop}.txt
subset -i rpkm.postqc.rds -o rpkm.postqc.topvargenes.rds -r /mnt/kauffman/edsgard/cloud/btsync/work/rspd/code/my/git/rrnaseq/test/rqc/diffexp/brennecke/rds/genes.ntop_${ntop}.txt


#########
#Plotting of most variable genes
#########

#Heatmap of most variable genes
heatmap -x 0.5 -i rpkm.postqc.topvargenes.rds -o /mnt/kauffman/edsgard/cloud/btsync/work/rspd/code/my/git/rrnaseq/test/rqc/diffexp/brennecke/pdf/heatmap.topvargenes.pdf

#PCA of most variable genes
pca -m meta.rds -r rpkm.postqc.topvargenes.rds -o /mnt/kauffman/edsgard/cloud/btsync/work/rspd/code/my/git/rrnaseq/test/rqc/diffexp/brennecke/pdf/pca.topvargenes.pdf -c 1,2,3
