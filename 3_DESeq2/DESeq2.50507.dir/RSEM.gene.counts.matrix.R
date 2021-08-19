library(cluster)
library(Biobase)
library(qvalue)
library(fastcluster)
options(stringsAsFactors = FALSE)
NO_REUSE = F

# try to reuse earlier-loaded data if possible
if (file.exists("RSEM.gene.counts.matrix.RData") && ! NO_REUSE) {
    print('RESTORING DATA FROM EARLIER ANALYSIS')
    load("RSEM.gene.counts.matrix.RData")
} else {
    print('Reading matrix file.')
    primary_data = read.table("../RSEM.gene.counts.matrix", header=T, com='', row.names=1, check.names=F, sep='\t')
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
data = data[rowSums(data)>=10,]
initial_matrix = data # store before doing various data transformations
cs = colSums(data)
data = t( t(data)/cs) * 1e6;
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

MA_plot = function(x, y, ...) {
    M = log( (exp(x) + exp(y)) / 2)
    A = x - y;
    res = list(x=M, y=A)
    return(res)
}
MA_color_fun = function(x,y) {
    col = sapply(y, function(y) ifelse(abs(y) >= 1, 'red', 'black')) # color 2-fold diffs
    return(col)
}
Scatter_color_fun = function(x,y) {
    col = sapply(abs(x-y), function(z) ifelse(z >= 1, 'red', 'black')) # color 2-fold diffs
    return(col)
}
for (i in 1:nsamples) {
    sample_name = sample_types[[i]]
    cat('Processing replicate QC analysis for sample: ', sample_name, "
")
    samples_want = sample_type_list[[sample_name]]
    samples_want = colnames(data) %in% samples_want
    if (sum(samples_want) > 1) {
        pdf(file=paste(sample_name, '.rep_compare.pdf', sep=''))
        d = data[,samples_want]
        initial_matrix_samples_want = initial_matrix[,samples_want]
        op <- par(mar = c(10,10,10,10))
        barplot(colSums(initial_matrix_samples_want), las=2, main=paste("Sum of Frags for replicates of:", sample_name), ylab='', cex.names=0.7)
        par(op)
        pairs3(d, pch='.', CustomColorFun=Scatter_color_fun, main=paste('Replicate Scatter:', sample_name)) # scatter plots
        pairs3(d, XY_convert_fun=MA_plot, CustomColorFun=MA_color_fun, pch='.', main=paste('Replicate MA:', sample_name)); # MA plots
        reps_cor = cor(d, method="pearson", use='pairwise.complete.obs')
        hc_samples = hclust(as.dist(1-reps_cor), method="complete")
        heatmap.3(reps_cor, dendrogram='both', Rowv=as.dendrogram(hc_samples), Colv=as.dendrogram(hc_samples), col = cm.colors(256), scale='none', symm=TRUE, key=TRUE,density.info='none', trace='none', symbreaks=F, margins=c(10,10), cexCol=1, cexRow=1, main=paste('Replicate Correlations:', sample_name) )
        dev.off()
    }
}
write.table(data, file="RSEM.gene.counts.matrix.minRow10.CPM.log2.centered.dat", quote=F, sep='	');
pdf("RSEM.gene.counts.matrix.minRow10.CPM.log2.centered.prcomp.principal_components.pdf")
prin_comp_data = data
pca = prcomp(prin_comp_data, center = FALSE, scale. = FALSE)
pc_pct_variance = (pca$sdev^2)/sum(pca$sdev^2)
def.par <- par(no.readonly = TRUE) # save default, for resetting...
gridlayout = matrix(c(1:4),nrow=2,ncol=2, byrow=TRUE);
layout(gridlayout, widths=c(1,1));
write.table(pca$rotation, file="RSEM.gene.counts.matrix.minRow10.CPM.log2.centered.PCA.prcomp.scores", quote=F, sep="	")
write.table(pca$x, file="RSEM.gene.counts.matrix.minRow10.CPM.log2.centered.PCA.prcomp.loadings", quote=F, sep="	")
PCA.loadings=pca$x
PCA.scores = pca$rotation
for (i in 1:(max(2,2)-1)) {
    xrange = range(PCA.scores[,i])
    yrange = range(PCA.scores[,i+1])
    samples_want = rownames(PCA.scores) %in% sample_type_list[[sample_types[1]]]
    pc_i_pct_var = sprintf("(%.2f%%)", pc_pct_variance[i]*100)
    pc_i_1_pct_var = sprintf("(%.2f%%)", pc_pct_variance[i+1]*100)
    plot(PCA.scores[samples_want,i], PCA.scores[samples_want,i+1], xlab=paste('PC',i, pc_i_pct_var), ylab=paste('PC',i+1, pc_i_1_pct_var), xlim=xrange, ylim=yrange, col=sample_colors[1])
    for (j in 2:nsamples) {
        samples_want = rownames(PCA.scores) %in% sample_type_list[[sample_types[j]]]
        points(PCA.scores[samples_want,i], PCA.scores[samples_want,i+1], col=sample_colors[j], pch=j)
    }
    plot.new()
    legend('topleft', as.vector(sample_types), col=sample_colors, pch=1:nsamples, ncol=2)
}

par(def.par)
pcloadings_mat_vals = PCA.loadings[,1:2]
print(dim(pcloadings_mat_vals))
pcloadings_mat = matrix_to_color_assignments(pcloadings_mat_vals, col=colorpanel(256,'purple','black','yellow'), by='col')
print(dim(pcloadings_mat))
colnames(pcloadings_mat) = paste('PC', 1:ncol(pcloadings_mat))
dev.off()
gene_cor = NULL
