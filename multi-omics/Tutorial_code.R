# See : https://compgenomr.github.io/book/use-case-multi-omics-data-from-colorectal-cancer.html
# see : https://nbisweden.github.io/workshop_omics_integration/session_nmf/NMF_main.html
# see : https://divingintogeneticsandgenomics.com/post/matrix-factorization-for-single-cell-rnaseq-data/

setwd("D:/CoBi/semester 1/Computational Biomedicine/Tutorial/Tutorial_13")

# increase the timeout for internet connection
options(timeout = 400)
#devtools::install_github("compgenomr/compGenomRData")
#install.packages("FactoMineR")
#install.packages("NMF")
library(compGenomRData)
library(ggplot2)
library(ggfortify)
library(FactoMineR)
library(NMF)
library(enrichR)
############ 1 ############ 
# read data and modify if needed
data_gex = readRDS("COREAD_CMS13_gex.rds")

data_wes = readRDS("COREAD_CMS13_muts.rds")
data_wes[data_wes > 0] = 1

data_cnv = readRDS("COREAD_CMS13_cnv.rds")

# read sample annotation and correct row names
metadata = read.csv("COREAD_CMS13_subtypes.csv", row.names = 1)
rownames(metadata) = gsub(pattern = "-", replacement = ".", rownames(metadata))

# create a dataframe of sample name and cms type for later annotation
metadata_cms = data.frame(cms = as.factor(metadata$cms_label))
rownames(metadata_cms) = rownames(metadata)

# PCA for each data
pca_gex = prcomp(t(data_gex), scale. = TRUE)
autoplot(pca_gex, data = metadata, color = "cms_label") + ggtitle("Scaled PCA for Gene Expression (GEX)")

pca_wes = prcomp(t(data_wes), scale. = TRUE)
autoplot(pca_wes, data = metadata, color = "cms_label") + ggtitle("Scaled PCA for Somatic Mutations (WES)")

pca_cnv = prcomp(t(data_cnv), scale. = TRUE)
autoplot(pca_cnv, data = metadata, color = "cms_label") + ggtitle("Scaled PCA for Copy Number Variation (CNV)")

# trick: you can see which features contribute to PC with the following
loadings = pca_gex$rotation[,2]
top_features = sort(abs(loadings), decreasing = T)

# heatmaps for each data "euclidean" distance
pheatmap(data_gex, main = "Heatmap for Gene Expression (GEX)",annotation_col = metadata_cms,
         clustering_distance_rows = "euclidean",
         clustering_distance_cols = "euclidean")

pheatmap(data_wes, main = "Heatmap for Somatic Mutatioms (WES)",annotation_col = metadata_cms,
         clustering_distance_rows = "euclidean",
         clustering_distance_cols = "euclidean")

pheatmap(data_cnv, main = "Heatmap for Copy Number Variations (CNV)",annotation_col = metadata_cms,
         clustering_distance_rows = "euclidean",
         clustering_distance_cols = "euclidean")

# heatmap for each data "pearson correlation" distance
pheatmap(data_gex, main = "Heatmap for Gene Expression (GEX)",annotation_col = metadata_cms,
         clustering_distance_rows = "correlation",
         clustering_distance_cols = "correlation")

pheatmap(data_wes, main = "Heatmap for Somatic Mutatioms (WES)",annotation_col = metadata_cms,
         clustering_distance_rows = "correlation",
         clustering_distance_cols = "correlation")

############ 2 ############
# join all data in one and transpose matrix
all_data = t(rbind(data_gex, data_wes, data_cnv))

# run MFA and define features
mfa = MFA(all_data, group = c(dim(data_gex)[1], dim(data_wes)[1], dim(data_cnv)[1]), 
          graph = FALSE, )

# extract decomposed matrices H and W
H_matrix = mfa$global.pca$ind$coord
W_matrix = mfa$quanti.var$coord

# add labelling to the dfs, be sure that names labels match the names in original frame
H_matrix_df = as.data.frame(H_matrix)
H_matrix_df$subtype = factor(metadata[rownames(H_matrix_df),]$cms_label)

# plot scatter of the two main dimensions/factors, 
ggplot2::ggplot(H_matrix_df, ggplot2::aes(x=Dim.1, y=Dim.2, color=subtype)) +
  ggplot2::geom_point() + ggplot2::ggtitle("Scatter plot of MFA")

# plot heatmap of the two main dimensions
pheatmap(t(H_matrix)[1:5,], annotation_col = metadata_cms, show_colnames = TRUE, 
         main = "MFA for multi-omics integration")

############ 3 ############
# first feature-normalize matrices, each row/rowsum
gex_featnorm = data_gex/rowSums(data_gex)
wes_featnorm = data_wes/rowSums(data_wes)
cnv_featnorm = data_cnv/rowSums(data_cnv)

# then normalize by omics type's Frobenius norm
gex_featnorm_frob = gex_featnorm/norm(as.matrix(gex_featnorm), type = "F")
wes_featnorm_frob = wes_featnorm/norm(as.matrix(wes_featnorm), type = "F")
cnv_featnorm_frob = cnv_featnorm/norm(as.matrix(cnv_featnorm), type = "F")

# run the Non-negative matrix factorization, rank 2 means
nmf_obj = nmf(t(rbind(gex_featnorm_frob,
                      wes_featnorm_frob,
                      cnv_featnorm_frob)),
              rank = 2, method = "Frobenius")

# extract W and H
nmf_H = NMF::basis(nmf_obj)
nmf_W = NMF::coef(nmf_obj)

# make a data frame of H and edit colnames 
nmf_H_df = as.data.frame(nmf_H)
colnames(nmf_H_df) = c("dim1","dim2")
nmf_H_df$subtype = factor(metadata[rownames(nmf_H_df),]$cms_label)

# plot how will NMF separates the two cancer subtypes
ggplot(nmf_H_df, aes(x=dim2, y=dim1, color = subtype)) +
geom_point() + ggtitle("Scatter plot of 2-factor NMF")

# plot heatmap for the same data, better than scatter in case of >2 factors
pheatmap(t(nmf_H_df[,1:2]), annotation_col = metadata_cms, show_colnames = TRUE, main = "Heatmap of 2-factor NMF")

############ 4 ############
# combine all row names of 3 datasets and their omics type
data_anno = data.frame(omics = 
  c(rep("gex", dim(data_gex)[1]),
  rep("wes", dim(data_wes)[1]),
  rep("cnv", dim(cnv_featnorm_frob)[1])))

# use that to annotate the features in the W matrix
nmf_W_t = t(nmf_W)
rownames(data_anno) = c(rownames(data_gex), paste0("mut:", rownames(data_wes)), rownames(cnv_featnorm_frob))
rownames(nmf_W_t) = rownames(data_anno)

# plot heatmap of feature weights in the W matrix
# there is a problem, i cant visualize features with highest contribution to factors
pheatmap(nmf_W_t[,], cluster_cols = FALSE, clustering_distance_rows = "manhattan",
         annotation_row = data_anno, main = "NMF coefficients", 
         show_rownames = TRUE, fontsize_row = 2)

# this part includes only data_gex, can be extended 
genes_factor_1 = names(which(nmf_W_t[1:dim(data_gex)[1],1] > nmf_W_t[1:dim(data_gex)[1],2]))
genes_factor_2 = names(which(nmf_W_t[1:dim(data_gex)[1],1] < nmf_W_t[1:dim(data_gex)[1],2]))

go_factor_1 = enrichr(genes_factor_1,databases = c("GO_Biological_Process_2023"))$GO_Biological_Process_2023

go_factor_2 = enrichr(genes_factor_2,
                                databases = c("GO_Biological_Process_2023")
)$GO_Biological_Process_2023
