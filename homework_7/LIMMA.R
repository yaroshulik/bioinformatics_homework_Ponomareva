BiocManager::install("limma")


library(limma)

BiocManager::install("Biobase")
library(Biobase)

geo_id <- 'GSE63885'
exp = read.csv(file = 'C:/Users/Александр/Downloads/Telegram Desktop/rna_seq_diff_exp/rna_seq_diff_exp/GSE63885/expression_for_limma2.csv', header = TRUE, row.names = 'Gene.Symbol')
ann = read.csv(file = 'C:/Users/Александр/Downloads/Telegram Desktop/rna_seq_diff_exp/rna_seq_diff_exp/GSE63885/annotation_for_limma.csv', header = TRUE, row.names = 'X')


exp_set <- ExpressionSet(assayData = as.matrix(exp), phenoData = AnnotatedDataFrame(ann))

slope <- factor(ann$clinical.status.post.1st.line.chemotherapy..cr...complete.response..pr...partial.response..sd...stable.disease..p...progression..ch1, levels = c("pCR","pNC"), labels = c(1, 0))

pCR=as.integer(as.vector(slope))

iterse <- rep(1, length(slope))

design <- cbind(npCR=iterse,pCR=pCR)

fit <- lmFit(exp_set, design)
fit <- eBayes(fit)
top <-topTable(fit, coef="pCR", adjust="BH", n = Inf)

thresholds <- c(1, 2, 3)

for (t in thresholds) {
  top[[paste0("signif_logFC", t)]] <- ifelse(
    abs(top$logFC) > t & top$adj.P.Val < 0.05,
    TRUE,
    FALSE
  )
}

sapply(thresholds, function(t) sum(top[[paste0("signif_logFC", t)]]))
library(ggplot2)
library(ggrepel)

# Добавляем имена генов как отдельный столбец
top$Gene.Symbol <- rownames(top)

thresholds <- c(1, 2, 3)

results <- data.frame(
  logFC_threshold = thresholds,
  significant_genes = NA
)

for (i in seq_along(thresholds)) {
  t <- thresholds[i]
  
  # Значимые гены по сырым p-values
  top$signif <- abs(top$logFC) > t & top$P.Value < 0.05
  
  # Сохраняем количество генов
  results$significant_genes[i] <- sum(top$signif)
  
  # Строим volcano plot
  p <- ggplot(top, aes(x = logFC, y = -log10(P.Value))) +
    geom_point(aes(color = signif), alpha = 0.7, size = 2.2) +
    scale_color_manual(values = c("grey70", "red")) +
    geom_vline(xintercept = c(-t, t), linetype = "dashed", color = "black") +
    geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "black") +
    
    # 🔥 Подписи генов
    geom_text_repel(
      data = subset(top, signif),
      aes(label = Gene.Symbol),
      size = 3,
      max.overlaps = 20,
      color = "red"
    ) +
    
    labs(
      title = paste("Volcano Plot — |logFC| >", t),
      subtitle = "Significance based on raw p-values (P.Value < 0.05)",
      x = "Log2 Fold Change (pCR vs pNC)",
      y = expression(-log[10](p-value))
    )
    theme_classic(base_size = 14)
  
  # Сохраняем график
  ggsave(
    filename = paste0("C:/Users/Александр/Downloads/volcano_logFC", t, "_rawP_final.png"),
    plot = p,
    width = 6,
    height = 5
  )
}

results

write.csv(top,
          file = "C:/Users/Александр/Downloads/limma_diffexp_results_with_significance.csv",
          row.names = TRUE)

