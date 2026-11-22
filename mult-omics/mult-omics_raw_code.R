library(compGenomRData)

# read in the csv from the companion package as a data frame
csvfile <- system.file("extdata", "multi-omics", "COREAD_CMS13_gex.csv", 
                       package="compGenomRData")
x1 <- read.csv(csvfile, row.names=1)
# Fix the gene names in the data frame
rownames(x1) <- sapply(strsplit(rownames(x1), "\\|"), function(x) x[1])
# Output a table
knitr::kable(head(t(head(x1))), caption="Example gene expression data (head)")

# read in the csv from the companion package as a data frame
csvfile <- system.file("extdata", "multi-omics", "COREAD_CMS13_muts.csv", 
                       package="compGenomRData")
x2 <- read.csv(csvfile, row.names=1)
# Set mutation data to be binary (so if a gene has more than 1 mutation,
# we only count one)
x2[x2>0]=1
# output a table
knitr::kable(head(t(head(x2))), caption="Example mutation data (head)")

# read in the csv from the companion package as a data frame
csvfile <- system.file("extdata", "multi-omics", "COREAD_CMS13_cnv.csv", 
                       package="compGenomRData")
x3 <- read.csv(csvfile, row.names=1)
head(x3)
# output a table
knitr::kable(head(t(head(x3))), 
             caption="Example copy number data for CRC samples")

# read in the csv from the companion package as a data frame
csvfile <- system.file("extdata", "multi-omics", "COREAD_CMS13_subtypes.csv",
                       package="compGenomRData")
covariates <- read.csv(csvfile, row.names=1)
# Fix the TCGA identifiers so they match up with the omics data
rownames(covariates) <- gsub(pattern = '-', replacement = '\\.',
                             rownames(covariates))
covariates <- covariates[colnames(x1),]
# create a dataframe which will be used to annotate later graphs
anno_col <- data.frame(cms=as.factor(covariates$cms_label))
rownames(anno_col) <- rownames(covariates)

# Load necessary packages
library(ggplot2)

# Perform PCA for gene expression data
pca_gex <- prcomp(t(x1), scale. = TRUE)
pca_muts <- prcomp(t(x2), scale. = TRUE)
pca_cnv <- prcomp(t(x3), scale. = TRUE)

# Convert PCA results into data frames for visualization
pca_gex_df <- data.frame(PC1 = pca_gex$x[,1], PC2 = pca_gex$x[,2], Group = anno_col$cms)
pca_muts_df <- data.frame(PC1 = pca_muts$x[,1], PC2 = pca_muts$x[,2], Group = anno_col$cms)
pca_cnv_df <- data.frame(PC1 = pca_cnv$x[,1], PC2 = pca_cnv$x[,2], Group = anno_col$cms)

# Function to plot PCA results
plot_pca <- function(pca_df, title) {
  ggplot(pca_df, aes(x = PC1, y = PC2, color = Group)) +
    geom_point(size = 3, alpha = 0.7) +
    labs(title = title, x = "PC1", y = "PC2") +
    theme_minimal()
}

# Plot PCA results
plot_pca(pca_gex_df, "PCA of Gene Expression Data")
plot_pca(pca_muts_df, "PCA of Mutation Data")
plot_pca(pca_cnv_df, "PCA of Copy Number Variation Data")


pheatmap::pheatmap(x1,
                   annotation_col = anno_col,
                   show_colnames = FALSE,
                   show_rownames = FALSE,
                   main="Gene expression data")

pheatmap::pheatmap(x2,
                   annotation_col = anno_col,
                   show_colnames = FALSE,
                   show_rownames = FALSE,
                   main="Mutation data")

pheatmap::pheatmap(x3,
                   annotation_col = anno_col,
                   show_colnames = FALSE,
                   show_rownames = FALSE,
                   main="Copy number data")

library(FactoMineR)
# run the MFA function from the FactoMineR package
r.mfa <- FactoMineR::MFA(
  t(rbind(x1,x2,x3)), # binding the omics types together
  c(dim(x1)[1], dim(x2)[1], dim(x3)[1]), # specifying the dimensions of each
  graph=FALSE, ncp = 10)

# first, extract the H and W matrices from the MFA run result
mfa.h <- r.mfa$global.pca$ind$coord
mfa.w <- r.mfa$quanti.var$coord

# create a dataframe with the H matrix and the CMS label
mfa_df <- as.data.frame(mfa.h)
mfa_df$subtype <- factor(covariates[rownames(mfa_df),]$cms_label)

# create the plot
ggplot2::ggplot(mfa_df, ggplot2::aes(x=Dim.1, y=Dim.2, color=subtype)) +
  ggplot2::geom_point() + ggplot2::ggtitle("Scatter plot of MFA")

pheatmap::pheatmap(t(mfa.h)[1:2,], annotation_col = anno_col,
                   show_colnames = FALSE,
                   main="MFA for multi-omics integration")

pheatmap::pheatmap(t(mfa.h)[1:10,], annotation_col = anno_col,
                   show_colnames = FALSE,
                   main="MFA: Heatmap of Leading 10 Factors")

# NMF
library(NMF)
split_neg_columns <- function(matrix_data) {
  # Get original column names
  original_colnames <- colnames(matrix_data)
  
  # Initialize matrices for positive and negative values
  pos_matrix <- matrix(0, nrow=nrow(matrix_data), ncol=ncol(matrix_data))
  neg_matrix <- matrix(0, nrow=nrow(matrix_data), ncol=ncol(matrix_data))
  
  # Preserve row names
  rownames(pos_matrix) <- rownames(matrix_data)
  rownames(neg_matrix) <- rownames(matrix_data)
  
  # Iterate over the matrix to separate positives and negatives
  for (j in 1:ncol(matrix_data)) {
    pos_matrix[, j] <- pmax(0, matrix_data[, j])  # Keep positives
    neg_matrix[, j] <- pmax(0, -matrix_data[, j]) # Convert negatives to positive
  }
  
  # Assign modified column names
  colnames(pos_matrix) <- paste0(original_colnames, "_pos")
  colnames(neg_matrix) <- paste0(original_colnames, "_neg")
  
  # Combine positive and negative matrices
  new_matrix <- cbind(pos_matrix, neg_matrix)
  
  return(new_matrix)
}


# Feature-normalize the data
x1.featnorm <- x1 / rowSums(x1)
x2.featnorm <- x2 / rowSums(x2)
x3.featnorm <- x3 / rowSums(x3)

# Normalize by each omics type's frobenius norm
x1.featnorm.frobnorm <- x1.featnorm / norm(as.matrix(x1.featnorm), type="F")
x2.featnorm.frobnorm <- x2.featnorm / norm(as.matrix(x2.featnorm), type="F")
x3.featnorm.frobnorm <- x3.featnorm / norm(as.matrix(x3.featnorm), type="F")

# Split the features of the CNV matrix into two non-negative features each
x3.featnorm.frobnorm.nonneg <- t(split_neg_columns(t(x3.featnorm.frobnorm)))

# run the nmf function from the NMF package
require(NMF) 

r.nmf <- nmf(t(rbind(x1.featnorm.frobnorm,
                     x2.featnorm.frobnorm,
                     x3.featnorm.frobnorm.nonneg)),
             2,
             method='Frobenius')

# exctract the H and W matrices from the nmf run result
nmf.h <- NMF::basis(r.nmf)
nmf.w <- NMF::coef(r.nmf)
nmfw <- t(nmf.w)

# create a dataframe with the H matrix and the CMS label (subtype)
nmf_df <- as.data.frame(nmf.h)
colnames(nmf_df) <- c("dim1", "dim2")
nmf_df$subtype <- factor(covariates[rownames(nmf_df),]$cms_label)

# create the scatter plot
ggplot2::ggplot(nmf_df, ggplot2::aes(x=dim1, y=dim2, color=subtype)) +
  ggplot2::geom_point() +
  ggplot2::ggtitle("Scatter plot of 2-component NMF for multi-omics integration")

pheatmap::pheatmap(t(nmf_df[,1:2]),
                   annotation_col = anno_col,
                   show_colnames=FALSE,
                   main="Heatmap of 2-component NMF")



data_anno <- data.frame(
  omics=c(rep('expression',dim(x1)[1]),
          rep('mut',dim(x2)[1]),
          rep('cnv',dim(x3.featnorm.frobnorm.nonneg)[1])))
rownames(data_anno) <- c(rownames(x1),
                         paste0("mut:", rownames(x2)),
                         rownames(x3.featnorm.frobnorm.nonneg))
rownames(nmfw) <- rownames(data_anno)

# generate the heat map
pheatmap::pheatmap(nmfw,
                   cluster_cols = FALSE,
                   annotation_row = data_anno,
                   main="NMF coefficients",
                   clustering_distance_rows = "manhattan",
                   fontsize_row = 1)

require(enrichR)

# select genes associated preferentially with each factor
# by their relative loading in the W matrix
genes.factor.1 <- names(which(nmfw[1:dim(x1)[1],1] > nmfw[1:dim(x1)[1],2]))
genes.factor.2 <- names(which(nmfw[1:dim(x1)[1],1] < nmfw[1:dim(x1)[1],2]))

# call the enrichr function to find gene sets enriched
# in each latent factor in the GO Biological Processes 2018 library
go.factor.1 <- enrichR::enrichr(genes.factor.1,
                                databases = c("GO_Biological_Process_2018")
)$GO_Biological_Process_2018
go.factor.2 <- enrichR::enrichr(genes.factor.2,
                                databases = c("GO_Biological_Process_2018")
)$GO_Biological_Process_2018

# create a data frame holding covariates (age, gender, MSI status)
a <- data.frame(age=covariates$age,
                gender=as.numeric(covariates$gender),
                msi=covariates$msi)

b <- nmf.h
colnames(b) <- c('factor1', 'factor2')

# concatenate the covariate dataframe with the H matrix
cov_factor <- cbind(a,b)

# generate the figure
ggplot2::ggplot(cov_factor, ggplot2::aes(x=msi, y=factor1, group=msi)) +
  ggplot2::geom_boxplot() +
  ggplot2::ggtitle("NMF factor 1 microsatellite instability")

ggplot2::ggplot(cov_factor, ggplot2::aes(x=msi, y=factor2, group=msi)) +
  ggplot2::geom_boxplot() +
  ggplot2::ggtitle("NMF factor 2 and microsatellite instability")

