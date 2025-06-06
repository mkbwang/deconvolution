

gene_expressions <- read.csv("experiment/train_data.csv", row.names=1)
library(stringr)

# extract ages
ages <- sapply(colnames(gene_expressions),
               function(cname) str_extract_all(cname, "\\d+")[[1]]) |>
    as.integer()

gene_expressions <- t(gene_expressions) |> as.data.frame()

# retain prevalent genes
prevalences <- colMeans(gene_expressions > 0)
subset_genes <- names(prevalences[prevalences == 1])
gene_expressions <- gene_expressions[, subset_genes]
# log_gene_expressions <- log(gene_expressions) # log expression
# gene_expressions_ranks <- apply(log_gene_expressions, 2, rank) |>
#     as.data.frame() # calculate ranks


outcome <- f_test_filter(geneexp_df=gene_expressions,
                         labels=ages,
                         min_label=2,
                         max_label=23)


subset_genes <- outcome$subset_genes
variance_df <- outcome$variance_subset
dendrogram <- hcluster_genes(geneexp_df=gene_expressions,
                             labels=ages,
                             subset_genes=subset_genes,
                             min_label=2,
                             max_label=23)


plot(dendrogram, labels = FALSE, xlab="Genes")

clusters <- cutree(dendrogram, k = 15)
variance_df$Cluster <- clusters

library(dplyr)
variance_df <- variance_df %>% group_by(Cluster) %>%
    arrange(pval)

# select top 2.5% of features
variance_subset_df <- variance_df %>% group_by(Cluster) %>%
    slice_head(prop = 1/40)

marker_genes <- variance_subset_df$Gene


ages_subset <- ages[ages >=2 & ages <= 23]
gene_expressions_subset <- gene_expressions[ages >=2 & ages <= 23, marker_genes] |>
    as.matrix()
gene_expressions_train <- gene_expressions_subset[-12, ]
gene_expressions_test <- gene_expressions_subset[12, ]

ages_train <- ages_subset[-12]
ages_test <- ages_subset[12]


l2_result <- deconvolution(X=gene_expressions_train, Y=gene_expressions_test,
                           labels=ages_train, type="l2")


l1_result <- deconvolution(X=gene_expressions_train, Y=gene_expressions_test,
                           labels=ages_train, type="l1")

