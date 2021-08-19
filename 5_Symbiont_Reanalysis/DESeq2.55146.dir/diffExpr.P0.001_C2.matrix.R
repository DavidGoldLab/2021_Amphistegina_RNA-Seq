library(cluster)
library(Biobase)
library(qvalue)
library(fastcluster)
options(stringsAsFactors = FALSE)
NO_REUSE = F

# try to reuse earlier-loaded data if possible
if (file.exists("diffExpr.P0.001_C2.matrix.RData") && ! NO_REUSE) {
    print('RESTORING DATA FROM EARLIER ANALYSIS')
    load("diffExpr.P0.001_C2.matrix.RData")
} else {
    print('Reading matrix file.')
    primary_data = read.table("diffExpr.P0.001_C2.matrix", header=T, com='', row.names=1, check.names=F, sep='\t')
    primary_data = as.matrix(primary_data)
}
source("/usr/local/Cellar/trinity/2.11.0/libexec/Analysis/DifferentialExpression/R/heatmap.3.R")
source("/usr/local/Cellar/trinity/2.11.0/libexec/Analysis/DifferentialExpression/R/misc_rnaseq_funcs.R")
source("/usr/local/Cellar/trinity/2.11.0/libexec/Analysis/DifferentialExpression/R/pairs3.R")
source("/usr/local/Cellar/trinity/2.11.0/libexec/Analysis/DifferentialExpression/R/vioplot2.R")
data = primary_data
myheatcol = colorpanel(75, 'purple','black','yellow')
samples_data = read.table("../samples_file.txt", header=F, check.names=F, fill=T)
samples_data = samples_data[samples_data[,2] != '',]
colnames(samples_data) = c('sample_name', 'replicate_name')
sample_types = as.character(unique(samples_data[,1]))
rep_names = as.character(samples_data[,2])
data = data[, colnames(data) %in% rep_names, drop=F ]
nsamples = length(sample_types)
sample_colors = rainbow(nsamples)
names(sample_colors) = sample_types
sample_type_list = list()
for (i in 1:nsamples) {
    samples_want = samples_data[samples_data[,1]==sample_types[i], 2]
    sample_type_list[[sample_types[i]]] = as.vector(samples_want)
}
sample_factoring = colnames(data)
for (i in 1:nsamples) {
    sample_type = sample_types[i]
    replicates_want = sample_type_list[[sample_type]]
    sample_factoring[ colnames(data) %in% replicates_want ] = sample_type
}
initial_matrix = data # store before doing various data transformations
data = log2(data+1)
sample_factoring = colnames(data)
for (i in 1:nsamples) {
    sample_type = sample_types[i]
    replicates_want = sample_type_list[[sample_type]]
    sample_factoring[ colnames(data) %in% replicates_want ] = sample_type
}
sampleAnnotations = matrix(ncol=ncol(data),nrow=nsamples)
for (i in 1:nsamples) {
  sampleAnnotations[i,] = colnames(data) %in% sample_type_list[[sample_types[i]]]
}
sampleAnnotations = apply(sampleAnnotations, 1:2, function(x) as.logical(x))
sampleAnnotations = sample_matrix_to_color_assignments(sampleAnnotations, col=sample_colors)
rownames(sampleAnnotations) = as.vector(sample_types)
colnames(sampleAnnotations) = colnames(data)
data = as.matrix(data) # convert to matrix

# Centering rows
data = t(scale(t(data), scale=F))

write.table(data, file="diffExpr.P0.001_C2.matrix.log2.centered.dat", quote=F, sep='	');
if (nrow(data) < 2) { stop("

**** Sorry, at least two rows are required for this matrix.

");}
if (ncol(data) < 2) { stop("

**** Sorry, at least two columns are required for this matrix.

");}
sample_cor = cor(data, method='pearson', use='pairwise.complete.obs')
write.table(sample_cor, file="diffExpr.P0.001_C2.matrix.log2.centered.sample_cor.dat", quote=F, sep='	')
sample_dist = dist(t(data), method='euclidean')
hc_samples = hclust(sample_dist, method='complete')
pdf("diffExpr.P0.001_C2.matrix.log2.centered.sample_cor_matrix.pdf")
sample_cor_for_plot = sample_cor
if (is.null(hc_samples)) { RowV=NULL; ColV=NULL} else { RowV=as.dendrogram(hc_samples); ColV=RowV }
heatmap.3(sample_cor_for_plot, dendrogram='both', Rowv=RowV, Colv=ColV, col = myheatcol, scale='none', symm=TRUE, key=TRUE,density.info='none', trace='none', symkey=FALSE, symbreaks=F, margins=c(10,10), cexCol=1, cexRow=1, cex.main=0.75, main=paste("sample correlation matrix
", "diffExpr.P0.001_C2.matrix.log2.centered") , ColSideColors=sampleAnnotations, RowSideColors=t(sampleAnnotations))
dev.off()
gene_cor = NULL
gene_dist = dist(data, method='euclidean')
if (nrow(data) <= 1) { message('Too few genes to generate heatmap'); quit(status=0); }
hc_genes = hclust(gene_dist, method='complete')
heatmap_data = data
pdf("diffExpr.P0.001_C2.matrix.log2.centered.genes_vs_samples_heatmap.pdf")
heatmap.3(heatmap_data, dendrogram='both', Rowv=as.dendrogram(hc_genes), Colv=as.dendrogram(hc_samples), col=myheatcol, scale="none", density.info="none", trace="none", key=TRUE, keysize=1.2, cexCol=1, margins=c(10,10), cex.main=0.75, main=paste("samples vs. features
", "diffExpr.P0.001_C2.matrix.log2.centered" ) , ColSideColors=sampleAnnotations)
dev.off()
save(list=ls(all=TRUE), file="diffExpr.P0.001_C2.matrix.RData")
