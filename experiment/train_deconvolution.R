

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

View(outcome$variance_subset)

subset_genes <- outcome$subset_genes
dendrogram <- hcluster_genes(geneexp_df=gene_expressions,
                             labels=ages,
                             subset_genes=subset_genes,
                             min_label=2,
                             max_label=23)


plot(dendrogram, labels = FALSE, xlab="Genes")

clusters <- cutree(dendrogram, k = 15)



