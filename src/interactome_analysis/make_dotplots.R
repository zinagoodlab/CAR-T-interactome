# file: make_dotplots.R
# author: Christine Yiwen Yeh
# ---------------------------
# script to make dotplots from cellphonedb outputs

library(reshape2)
library(dplyr)
library(ggplot2)

getwd()
setwd("~/Desktop/CellPhoneDB_Results")

#############
# FUNCTIONS #
#############

'%!in%' <- function(x,y)!('%in%'(x,y))

#' Magic Shortcut for "not in"
#'
#' This function is a negation of the R default magic
#' `%in%` it return a boolean of left hand side list of
#' whether each element is in the right hand side list.
plot_dotplot_broad <- function(tissue, plot_title, dir) {
  meansfilepath <- Sys.glob(paste0(dir, tissue,  "/means.txt"))
  means_m <- melt_cpdb_data(meansfilepath, "mean_expression")
  means_m$log2_mean_expression <- log2(means_m$mean_expression)

  pvalsfilepath <- Sys.glob(paste0(dir, tissue, "/pvalues.txt"))
  pvals_m <- melt_cpdb_data(pvalsfilepath, "pvalues")

  # bin pvalues
  pvals_m$neg_log10_pval <- unlist(lapply(pvals_m$pvalues, function(x){
    x_prime <- -log10(x)
    if (x_prime >= 3) {
      return (3)
    } else if (x_prime>=2) {
      return (2)
    } else if (x_prime >= 1) {
      return (1)
    } else {
      return (0)
    }
  }))

  # find order of genes
  gene_lists <- find_gene_levels(pvals_m)

  # merge
  dotplot_df <- merge(means_m,
                      pvals_m,
                      by = c("interacting_gene_pair", "interacting_cells"))
  #dotplot_df <- dotplot_df[grepl("CCL18", dotplot_df$interacting_gene_pair), ]
  dotplot_df$interacting_gene_pair <- factor(
    dotplot_df$interacting_gene_pair,
    levels = gene_lists[[1]])
  dotplot_df <- dotplot_df %>%
    mutate(interacting_cells = gsub("\\.", "_", interacting_cells))
  dotplot_df <- dotplot_df %>%
    mutate(interacting_cells = gsub("\\|", " > ", interacting_cells))

  plot_height <- 0.3*length(unique(filter(dotplot_df,
                                          interacting_gene_pair %!in%
                                            gene_lists[[2]])$interacting_gene_pair))

  # plot pdf
  pdf(file = paste0("CCL18_", tissue, "_sig_dotplot.pdf"), width = 10, height = 3.5) # 8, 3.5 
  p = ggplot(filter(dotplot_df, interacting_gene_pair %!in% gene_lists[[2]]) ,
             aes(x = interacting_cells, y = interacting_gene_pair,
                 col = log2_mean_expression, size = neg_log10_pval)) +
    geom_point() +
    scale_color_gradient2(low = "#E4E4E4",
                          mid ="#FFDD1D",
                          high = "#D10000", midpoint = -5) +
    scale_size(range = c(1, 4)) +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 90),
          legend.position = "top") +
    ggtitle(paste0("Significant Ligand-Receptor Pairs between ", plot_title)) +
    xlab("Interacting Cell [celltype1 > celltype2]") +
    ylab("Interacting Gene Pair [gene1_gene2]") +
    guides(color=guide_colorbar(title="Log2 Mean Expression", title.position="top"),
           size =guide_legend(title = "-Log10 P-Value", title.position = "top"))
  print(p)
  dev.off()
  ggsave(paste0(tissue,".png"), p, dpi=600, width = 8, height = 8) #8, 3.5

  dotplot_df$sample <- plot_title
  return(dotplot_df)
}

#' Get Gene Expression Levels
#'
#' This function takes the pvalue long format dataframe
#' and summarizes the data such that we obtain the
#' summary statistics
#'
#' @param pvals_m long format pvalues dataframe
#' @returns gene levels
#' @export
find_gene_levels <- function(pvals_m) {
  pvalue_summary <- pvals_m %>%
    group_by(interacting_gene_pair) %>%
    summarise(sum_pvals = sum(neg_log10_pval),
              max_pval = max(neg_log10_pval))

  gene_levels <- pvalue_summary[
    order(
      pvalue_summary$sum_pvals, decreasing = F
      ),
    ]$interacting_gene_pair

  zeros <- filter(pvalue_summary,
                  sum_pvals == 0 | max_pval <=1)$interacting_gene_pair

  return(list(gene_levels, zeros))
}

#' Melt CellPhoneDB Data
#'
#' This function takes the raw outputs from cellphoneDB and converts it to
#' long format for downstream analyses and plotting
#'
#' @param filepath filepath to the  outputs from cellphonedb
#' @param name name of the final column
#' @returns melted dataframe
#' @export
melt_cpdb_data <- function(filepath, name){
  df <- read.csv(filepath, sep = "\t", header = T, check.names = FALSE)
  interacting_gene_pair <- df$interacting_pair
  df <- df[,c(12:dim(df)[2])]
  df <- data.frame(interacting_gene_pair, df, check.names = FALSE)
  df_m <- melt(df, variable.name = "interacting_cells", value.name = name)
  return(df_m)
}

##############
# PRODUCTION #
##############

# BASE V1.0
dotploti <- plot_dotplot_broad("CR",
                               "30d_CR",
                               "Base_v1.0/30d/")
dotplotii <- plot_dotplot_broad("PD",
                               "30d_PD",
                               "Base_v1.0/30d/")
dotplot1 <- plot_dotplot_broad("CR",
                               "3m_CR",
                               "Base_v1.0/3m/")
dotplot2 <- plot_dotplot_broad("PD",
                               "3m_PD",
                               "Base_v1.0/3m/")
dotplot3 <- plot_dotplot_broad("CR",
                               "6m_CR",
                               "Base_v1.0/6m/")
dotplot4 <- plot_dotplot_broad("PD",
                               "6m_PD",
                               "Base_v1.0/6m/")

# TREG V2.0
dotploti <- plot_dotplot_broad("CR",
                               "30d_CR",
                               "Treg_v2.0/30d/")
dotplotii <- plot_dotplot_broad("PD",
                                "30d_PD",
                                "Treg_v2.0/30d/")
dotplot1 <- plot_dotplot_broad("CR",
                               "3m_CR",
                               "Treg_v2.0/3m/")
dotplot2 <- plot_dotplot_broad("PD",
                               "3m_PD",
                               "Treg_v2.0/3m/")
dotplot3 <- plot_dotplot_broad("CR",
                               "6m_CR",
                               "Treg_v2.0/6m/")
dotplot4 <- plot_dotplot_broad("PD",
                               "6m_PD",
                               "Treg_v2.0/6m/")

dotploti + dotplotii


# CAR V3.0
dotploti <- plot_dotplot_broad("CR",
                               "30d_CR",
                               "CAR_v3.0/30d/")
dotplotii <- plot_dotplot_broad("PD",
                                "30d_PD",
                                "CAR_v3.0/30d/")
dotplot1 <- plot_dotplot_broad("CR",
                               "3m_CR",
                               "CAR_v3.0/3m/")
dotplot2 <- plot_dotplot_broad("PD",
                               "3m_PD",
                               "CAR_v3.0/3m/")
dotplot3 <- plot_dotplot_broad("CR",
                               "6m_CR",
                               "CAR_v3.0/6m/")
dotplot4 <- plot_dotplot_broad("PD",
                               "6m_PD",
                               "CAR_v3.0/6m/")
