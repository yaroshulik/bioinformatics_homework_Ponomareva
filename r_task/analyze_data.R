data <- read.csv("sample_data.csv")

cat("Структура данных:\n")
str(data)
cat("\n")

mean_score <- mean(data$Score)
cat("Среднее значение Score для всех образцов:", mean_score, "\n")

max_score_treatment <- max(data[data$Group == "Treatment", ]$Score)
cat("Максимальное значение Score в группе Treatment:", max_score_treatment, "\n")
cat("\n")

png("score_boxplot.png", width = 800, height = 600)
boxplot(Score ~ Group, 
data = data,
main = "SCore Distribution By Group",
xlab = "Group",
ylab = "Score",
col = c("skyblue", "lightcoral"),
border = "darkblue",
notch = FALSE,
vsrwidth = TRUE)

stripchart(Score ~ Group,
data = data,
method = "jitter",
vertical = TRUE,
add = TRUE,
pch = 19,
col = c("darkblue", "darkred"))

grid(nx = NA, ny = NULL, col = "lightgray", lty = "dotted")

dev.off()

cat("График сохранен в файл: score_boxplot.png\n")
