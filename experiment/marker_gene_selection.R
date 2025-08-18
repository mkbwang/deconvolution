
rm(list=ls())
gene_expressions <- read.csv("experiment/train_data.csv", row.names=1)
metadata <- read.csv("experiment/metadata.csv", row.names=1)
library(stringr)
library(dplyr)

# extract ages
ages <- sapply(colnames(gene_expressions),
               function(cname) str_extract_all(cname, "\\d+")[[1]]) |>
    as.integer()

gene_expressions <- t(gene_expressions) |> as.data.frame()
age_filter <- ages >= 2 & ages <= 23
gene_expressions <- gene_expressions[age_filter, ]
ages <- ages[age_filter]


# retain prevalent genes
prevalences <- colMeans(gene_expressions > 0)
subset_genes <- names(prevalences[prevalences == 1])
gene_expressions <- gene_expressions[, subset_genes]

gene_expressions$Age <- ages
gene_expressions_meanbyage <- gene_expressions %>% group_by(Age) %>%
    summarise_all(median)

outcome <- f_test_filter(geneexp_df=gene_expressions,
                         labels=ages,
                         min_label=6,
                         max_label=23)



subset_genes <- outcome$subset_genes
variance_df <- outcome$variance
variance_df <- variance_df %>% arrange(pval)
markergene_anova <- variance_df[subset_genes, ]




# cluster genes
dendrogram <- hcluster_genes(geneexp_df=gene_expressions,
                             labels=ages,
                             subset_genes=subset_genes,
                             min_label=6,
                             max_label=23)




plot(dendrogram, labels = FALSE, xlab="Genes")
clusters <- cutree(dendrogram, k = 15)
markergene_anova$Cluster <- clusters
markergene_anova <- markergene_anova %>% arrange(Cluster, pval)



# select top 2.5% of features
# variance_subset_df <- variance_df %>% group_by(Cluster) %>%
#     slice_head(prop = 1/40)
write.csv(markergene_anova, "experiment/markergene_anova.csv", row.names=FALSE,
          quote=FALSE)


gene_expressions_subset <- gene_expressions[, subset_genes]

# check linear association
assoc_2_18 <- twas(ages[ages >= 2 & ages <= 18],
                   gene_expressions_subset[ages >= 2 & ages <= 18, ],
                   takelog=T)
colnames(assoc_2_18)[c(2,3)] <- c("effsize_2_18", "pval_2_18")



assoc_6_15 <- twas(ages[ages >= 6 & ages <= 15],
                   gene_expressions_subset[ages >= 6 & ages <= 15, ],
                   takelog=T)
colnames(assoc_6_15)[c(2,3)] <- c("effsize_6_15", "pval_6_15")


assoc_6_18 <- twas(ages[ages >= 6 & ages <= 18],
                   gene_expressions_subset[ages >= 6 & ages <= 18, ],
                   takelog=T)
colnames(assoc_6_18)[c(2,3)] <- c("effsize_6_18", "pval_6_18")


assoc_6_23 <- twas(ages[ages >= 6 & ages <= 23],
                   gene_expressions_subset[ages >= 6 & ages <= 23, ],
                   takelog=T)
colnames(assoc_6_23)[c(2,3)] <- c("effsize_6_23", "pval_6_23")

markergene_association <- markergene_anova %>% select(Gene, Cluster) %>% left_join(assoc_2_18, by="Gene") %>%
    left_join(assoc_6_15, by="Gene") %>%
    left_join(assoc_6_18, by="Gene") %>%
    left_join(assoc_6_23, by="Gene")

write.csv(markergene_association, "experiment/markergene_association.csv", row.names=FALSE,
          quote=FALSE)


# check monotonicity
monotonicity <- monotone_check(ages=ages,
                               expression_df=gene_expressions_subset)

monotonicity <- (monotonicity-0.5)*2
monotonicity$Gene <- rownames(monotonicity)

markergene_monotonicity <- markergene_anova %>% select(Gene, Cluster) %>%
    left_join(monotonicity, by="Gene")

write.csv(markergene_monotonicity, "experiment/markergene_monotonicity.csv", row.names=FALSE,
          quote=FALSE)

# visualize some genes

example_geneexp <- gene_expressions[, "SSMa011080"]
example_plot <- viz_gene(ages=ages, values=example_geneexp,
                         title="SSMa011080")



