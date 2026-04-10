genes <- c("BRCA1", "TP53", "EGFR")
expression <- c(12.5, 45.2, 30.1)
condition <- c("Control", "Treatment", "Treatment")

exp_data <- data.frame(genes, expression, condition)

str(exp_data)

png("expression_plot.png", width = 800, height = 600)
barplot(expression,
names.arg = genes,
col = c("skyblue", "lightcoral", "lightgreen"),
main = "Gene Expression Levels",
xlab = "Genes",
ylab = "Expression Value",
ylim = c(0,50))

text(x = 1:3,
y = expression + 2,
labels = expression,
cex = 0.8)
dev.off()

cat("График сохранен в файл expression_plot.png\n")
 
