
rm(list=ls())
gene_expressions <- read.csv("experiment/train_data.csv", row.names=1)
gene_expressions_new <- read.csv("experiment/train_data_new.csv", row.names=1)
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
gene_expressions_medianbyage <- gene_expressions %>% group_by(Age) %>%
    summarise_all(median)
gene_expressions$Age <- NULL

anova_log <- f_test_filter(geneexp_df=gene_expressions,
                         labels=ages,
                         transformation = "Log",
                         min_label=6,
                         max_label=23)

anova_rank <- f_test_filter(geneexp_df=gene_expressions,
                            labels=ages,
                            transformation = "Rank",
                            min_label=6,
                            max_label=23)


log_variance_df <- anova_log$variance %>% arrange(pval) %>% filter(adjusted_pvals < 0.05)
rank_variance_df <- anova_rank$variance %>% arrange(pval) %>% filter(adjusted_pvals < 0.05)

markergenes <- intersect(log_variance_df$Gene, rank_variance_df$Gene)

# length(intersect(anova_log$subset_genes, anova_rank$subset_genes))
# View(log_variance_df[anova_rank$subset_genes, ])
# View(rank_variance_df[anova_log$subset_genes, ])
#
# markergene_anova <- variance_df[subset_genes, ]




# cluster genes
# dendrogram <- hcluster_genes(geneexp_df=gene_expressions,
#                              labels=ages,
#                              subset_genes=markergenes,
#                              min_label=6,
#                              max_label=23)



# check the variance of ranks at each age for each gene
gene_expressions_subset <- gene_expressions[, markergenes]
gene_expressions_new_subset <- gene_expressions_new[markergenes, ] |> t() |> as.data.frame()

gene_expressions_subset_rank <- apply(gene_expressions_subset, 2, rank) |>
    as.data.frame() # calculate ranks
gene_expressions_subset_rank$Age <- ages
rank_var_byage <- gene_expressions_subset_rank %>% group_by(Age) %>% summarise_all(var) |> as.data.frame()
rownames(rank_var_byage) <- sprintf("Age_%d", rank_var_byage$Age)
rank_var_byage$Age <- NULL
rank_var_byage <- t(rank_var_byage) |> as.data.frame()
unique_ages <- unique(ages)


rank_var_byage <- rank_var_byage %>% arrange(Age_2)
rank_var_byage_sort <- apply(rank_var_byage, 2, function(x) rank(x, ties.method="min")) |>
    as.data.frame()


useful_genes_by_age <- list()
for (age in unique(ages)){
    cname <- sprintf("Age_%d", age)
    useful_genes_by_age[[cname]] <- rownames(rank_var_byage)[which(rank_var_byage_sort[, cname] <=10)]
}

all_genes <- unlist(useful_genes_by_age)
gene_counts <- table(all_genes)

log_gene_expressions <- log(gene_expressions)


test_data <- read.csv("experiment/test_data_new.csv", row.names=1) |> t() |> as.data.frame()
test_data <- test_data[, markergenes]
metadata_test <- read.csv("experiment/metadata.csv", row.names=1)
test_ages <- tail(metadata_test$chronological.age, 8) |> as.integer()

library(ggplot2)
# visualize these genes
for (age in c(6,8,10,12,15,17)){

    genes <- useful_genes_by_age[[sprintf("Age_%d", age)]]
    output_folder_old <- sprintf("experiment/plots/markergene_%d_old", age)
    output_folder_new <- sprintf("experiment/plots/markergene_%d_new", age)
    if(!dir.exists(output_folder_old)){
        dir.create(output_folder_old)
        dir.create(output_folder_new)
    }

    for (gene in genes){
        train_geneexp_old <- gene_expressions_subset[, gene]
        train_exp_df_old <- data.frame(Expression=train_geneexp_old,
                                   Age=ages,
                                   Label="Train")

        train_geneexp_new <- gene_expressions_new_subset[, gene]
        train_exp_df_new <- data.frame(Expression=train_geneexp_new,
                                       Age=ages,
                                       Label="Train")

        test_geneexp <- test_data[, gene]
        test_exp_df <- data.frame(Expression=test_geneexp,
                                  Age=test_ages,
                                  Label="Test")

        combined_df_old <- rbind(train_exp_df_old, test_exp_df)
        geneexp_plot_old <- ggplot(combined_df_old, aes(x=Age,
                                                y=Expression, color=Label)) +
            geom_point(alpha=0.7) + scale_x_continuous(breaks=unique(combined_df_old$Age)) +
            scale_color_manual(values=c("blue", "black")) + theme(legend.position = "none")+
            xlab("Age") +ylab("Expression") + ggtitle(gene)
        ggsave(filename=file.path(output_folder_old, sprintf("%s.pdf", gene)),
               plot=geneexp_plot_old, width=6, height=4)


        combined_df_new <- rbind(train_exp_df_new, test_exp_df)
        geneexp_plot_new <- ggplot(combined_df_new, aes(x=Age,
                                                        y=Expression, color=Label)) +
            geom_point(alpha=0.7) + scale_x_continuous(breaks=unique(combined_df_new$Age)) +
            scale_color_manual(values=c("blue", "black")) + theme(legend.position = "none")+
            xlab("Age") +ylab("Expression") + ggtitle(gene)
        ggsave(filename=file.path(output_folder_new, sprintf("%s.pdf", gene)),
               plot=geneexp_plot_new, width=6, height=4)

    }
}

# plot(dendrogram, labels = FALSE, xlab="Genes")
# clusters <- cutree(dendrogram, k = 15)
# markergene_anova$Cluster <- clusters
# markergene_anova <- markergene_anova %>% arrange(Cluster, pval)



# select top 2.5% of features
# variance_subset_df <- variance_df %>% group_by(Cluster) %>%
#     slice_head(prop = 1/40)
write.csv(markergene_anova, "experiment/markergene_anova.csv", row.names=FALSE,
          quote=FALSE)


# gene_expressions_subset <- gene_expressions[, subset_genes]
#
# # check linear association
# assoc_2_18 <- twas(ages[ages >= 2 & ages <= 18],
#                    gene_expressions_subset[ages >= 2 & ages <= 18, ],
#                    takelog=T)
# colnames(assoc_2_18)[c(2,3)] <- c("effsize_2_18", "pval_2_18")
#
#
#
# assoc_6_15 <- twas(ages[ages >= 6 & ages <= 15],
#                    gene_expressions_subset[ages >= 6 & ages <= 15, ],
#                    takelog=T)
# colnames(assoc_6_15)[c(2,3)] <- c("effsize_6_15", "pval_6_15")
#
#
# assoc_6_18 <- twas(ages[ages >= 6 & ages <= 18],
#                    gene_expressions_subset[ages >= 6 & ages <= 18, ],
#                    takelog=T)
# colnames(assoc_6_18)[c(2,3)] <- c("effsize_6_18", "pval_6_18")
#
#
# assoc_6_23 <- twas(ages[ages >= 6 & ages <= 23],
#                    gene_expressions_subset[ages >= 6 & ages <= 23, ],
#                    takelog=T)
# colnames(assoc_6_23)[c(2,3)] <- c("effsize_6_23", "pval_6_23")
#
# markergene_association <- markergene_anova %>% select(Gene, Cluster) %>% left_join(assoc_2_18, by="Gene") %>%
#     left_join(assoc_6_15, by="Gene") %>%
#     left_join(assoc_6_18, by="Gene") %>%
#     left_join(assoc_6_23, by="Gene")
#
# write.csv(markergene_association, "experiment/markergene_association.csv", row.names=FALSE,
#           quote=FALSE)
#
#
# # check monotonicity
# monotonicity <- monotone_check(ages=ages,
#                                expression_df=gene_expressions_subset)
#
# monotonicity <- (monotonicity-0.5)*2
# monotonicity$Gene <- rownames(monotonicity)
#
# markergene_monotonicity <- markergene_anova %>% select(Gene, Cluster) %>%
#     left_join(monotonicity, by="Gene")
#
# write.csv(markergene_monotonicity, "experiment/markergene_monotonicity.csv", row.names=FALSE,
#           quote=FALSE)
#
# # visualize some genes
#




