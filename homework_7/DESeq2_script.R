if (!require("BiocManager", quietly = TRUE)) install.packages("BiocManager")
BiocManager::install("DESeq2")
BiocManager::install("dplyr")
BiocManager::install("ggplot2")
BiocManager::install("ggrepel")
install.packages("dplyr")
install.packages("ggrepel")
library(DESeq2)
library(dplyr)
library(ggplot2)
library(ggrepel)
library(dplyr)


expr_raw <- read.csv("C:/Users/Александр/Downloads/Telegram Desktop/rna_seq_diff_exp/rna_seq_diff_exp/raw_counts_ici_samples.tsv",
                     sep = '\t', row.names = 1)

meta = data.frame(RNA = colnames(expr_raw))

dds <- DESeqDataSetFromMatrix(countData = expr_raw,
                              colData = meta,
                              design = ~ 1)

dds <- estimateSizeFactors(dds)
normalized_counts <- counts(dds, normalized = TRUE)

write.table(normalized_counts,
            file = "C:/Users/Александр/Downloads/normalized_counts_ici_samples.tsv",
            sep = "\t", quote = FALSE, col.names = NA)

sample_metadata <- read.csv("C:/Users/Александр/Downloads/Telegram Desktop/rna_seq_diff_exp/rna_seq_diff_exp/meta_responses.tsv", row.names = 1, sep='\t')
sample_metadata <- sample_metadata %>%
  filter(X0 %in% c('R', 'NR'))
rownames(sample_metadata) <- gsub("-", ".", rownames(sample_metadata))

normalized_counts <- read.csv('C:/Users/Александр/Downloads/normalized_counts_ici_samples.tsv', sep='\t', row.names = 1)

normalized_counts <- normalized_counts[row.names(sample_metadata)]

normalized_counts <- round(normalized_counts)


dds <- DESeqDataSetFromMatrix(countData = normalized_counts,
                              colData = sample_metadata,
                              design = ~ X0)

dds$X0 <- relevel(dds$X0, ref = "NR")
dds <- DESeq(dds)

res <- results(dds)
summary(res)

# --- Преобразуем результаты в data.frame ---
res_df <- as.data.frame(res)

# Убираем NA, чтобы не ломать графики
res_df <- res_df[!is.na(res_df$padj) & !is.na(res_df$log2FoldChange), ]

# --- Добавляем колонку Gene (HGNC символы) ---
hgnc_table <- read.csv(
  "C:/Users/Александр/Downloads/Telegram Desktop/rna_seq_diff_exp/rna_seq_diff_exp/hgnc_complete_set.txt",
  sep = "\t", row.names = 1
)

symbol_map <- setNames(hgnc_table$symbol, hgnc_table$ensembl_gene_id)

res_df$Gene <- sapply(rownames(res_df), function(id) {
  if (id %in% names(symbol_map) && !is.na(symbol_map[id])) {
    symbol_map[id]
  } else {
    id
  }
})

# --- Порог значимости ---
pvalue_threshold <- 0.01
log2fc_threshold <- 3

res_df$significance <- ifelse(
  res_df$padj < pvalue_threshold & abs(res_df$log2FoldChange) > log2fc_threshold,
  "Significant",
  "Not Significant"
)

# --- Volcano plot ---
library(ggplot2)
library(ggrepel)

pvalue_threshold <- 0.01
log2fc_threshold <- 3

# Метка значимости
res_df$significance <- ifelse(
  res_df$padj < pvalue_threshold & abs(res_df$log2FoldChange) > log2fc_threshold,
  "Significant", "Not Significant"
)

# Выбираем топ-15 значимых генов для подписей
top_labels <- res_df %>%
  dplyr::filter(significance == "Significant") %>%
  dplyr::arrange(padj) %>%
  head(15)


library(ggplot2)
library(ggrepel)

volcano_plot <- ggplot(res_df,
                       aes(x = log2FoldChange,
                           y = -log10(padj),
                           color = significance)) +
  geom_point(alpha = 0.8, size = 1.5) +
  scale_color_manual(values = c("grey", "red")) +
  labs(title = "Volcano plot: R vs NR",
       x = "Log2 Fold Change (R vs NR)",
       y = "-Log10 adjusted p-value") +
  theme_minimal(base_size = 12) +
  theme(
    legend.position = "none",
    plot.title = element_text(hjust = 0.5)
  ) +
  geom_hline(yintercept = -log10(pvalue_threshold),
             linetype = "dashed", color = "black", linewidth = 0.4) +
  geom_vline(xintercept = c(-log2fc_threshold, log2fc_threshold),
             linetype = "dashed", color = "black", linewidth = 0.4) +
  geom_text_repel(
    data = top_labels,
    aes(label = Gene),
    size = 3,
    max.overlaps = 50,
    box.padding = 0.4,
    point.padding = 0.2,
    segment.size = 0.3
  ) +
  xlim(-10, 10)

volcano_plot

# --- Сохраняем таблицу ---
write.csv(
  res_df[, c("Gene", "log2FoldChange", "pvalue", "padj")],
  "C:/Users/Александр/Downloads/differential_expression_results_ICI.csv",
  row.names = FALSE
)
