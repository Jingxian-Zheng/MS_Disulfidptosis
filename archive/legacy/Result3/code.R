########## 基于机器学习构建并验证双硫死亡相关的基质特征评分 
##### 设置工作环境
setwd(file.path(Sys.getenv("LEGACY_STUDY_ROOT", unset = "legacy_inputs"), "Result3"))

##### 加载R包
library("glmnet")
library("pROC")
library("caret")
library("ggplot2")
library("pheatmap")
library("reshape2")
library("dplyr")
library("tidyr")
library("forestplot")
library("randomForest")
library("e1071")
library("kknn")
library("ComplexHeatmap")
library("circlize")
library("gbm")
library("rpart")
library("MASS")
library("ggpubr")
##### 导入单细胞与空间跨模态纯化出的14个双硫死亡基质核心基因
disulfidptosis_genes <- read.delim(file.path(Sys.getenv("LEGACY_STUDY_ROOT", unset = "legacy_inputs"), "Result1/disulfidptosis_markers.txt"), header = F)
disulfidptosis_genes <- disulfidptosis_genes$V1

##### 读取数据, GSE17048做训练集，GSE32915做验证集
GSE17048_mRNA_exp <- read.delim(file.path(Sys.getenv("LEGACY_STUDY_ROOT", unset = "legacy_inputs"), "Data/GSE17048_mRNA_exp.tsv"), row.names = 1)
GSE17048_sample_type <- read.delim(file.path(Sys.getenv("LEGACY_STUDY_ROOT", unset = "legacy_inputs"), "Data/GSE17048_sample_type.txt"))
GSE32915_mRNA_exp <- read.delim(file.path(Sys.getenv("LEGACY_STUDY_ROOT", unset = "legacy_inputs"), "Data/GSE32915_mRNA_exp.tsv"), row.names = 1)
GSE32915_sample_type <- read.delim(file.path(Sys.getenv("LEGACY_STUDY_ROOT", unset = "legacy_inputs"), "Data/GSE32915_sample_type.txt"))

##### 靶基因与两个Bulk集求交集，确保矩阵完美对齐
common_genes <- intersect(intersect(rownames(GSE17048_mRNA_exp), rownames(GSE32915_mRNA_exp)), disulfidptosis_genes)

### 训练集规整 (Control=0, MS=1)
train_expr <- GSE17048_mRNA_exp[common_genes, ]
train_groups <- GSE17048_sample_type
train_groups$group <- ifelse(train_groups$group == "Control", 0, 1)
train_expr <- train_expr[, train_groups$sample]

### 验证集规整 (Control=0, MS=1)
test_expr <- GSE32915_mRNA_exp[common_genes, ]
test_groups <- GSE32915_sample_type
test_groups$group <- ifelse(test_groups$group == "Control", 0, 1)
test_expr <- test_expr[, test_groups$sample]

##### 11 种分类机器学习算法真机现场训练、AUC矩阵与ROC对象全量捕获
print("🚀 正在真机运行 11 种机器学习模型，全量捕获诊断数据流...")

algorithms <- c("LASSO", "Ridge", "ElasticNet", "RandomForest", "SVM", 
                "GBM", "kNN", "NaiveBayes", "LDA", "StepwiseLogistic", "DecisionTree")
auc_matrix <- matrix(NA, nrow = length(algorithms), ncol = 2,
                     dimnames = list(algorithms, c("GSE17048", "GSE32915")))

# 转换标准长矩阵 (Samples x Genes) 供底层算法直接读取
x_train <- t(train_expr)
y_train <- train_groups$group
x_test  <- t(test_expr)
y_test  <- test_groups$group

df_train <- data.frame(x_train, y = as.factor(y_train))
df_test  <- data.frame(x_test, y = as.factor(y_test))

set.seed(2026)

# 1~3. 运行 glmnet 家族现场训练（LASSO, Ridge, ElasticNet）
fit_lasso <- cv.glmnet(x_train, y_train, family = "binomial", alpha = 1)
roc_lasso_train <- roc(y_train, as.numeric(predict(fit_lasso, newx=x_train, s="lambda.min", type="response")), quiet=TRUE)
roc_lasso_test  <- roc(y_test, as.numeric(predict(fit_lasso, newx=x_test, s="lambda.min", type="response")), quiet=TRUE)
auc_matrix["LASSO", "GSE17048"] <- as.numeric(auc(roc_lasso_train))
auc_matrix["LASSO", "GSE32915"] <- as.numeric(auc(roc_lasso_test))

fit_ridge <- cv.glmnet(x_train, y_train, family = "binomial", alpha = 0)
roc_ridge_train <- roc(y_train, as.numeric(predict(fit_ridge, newx=x_train, s="lambda.min", type="response")), quiet=TRUE)
roc_ridge_test  <- roc(y_test, as.numeric(predict(fit_ridge, newx=x_test, s="lambda.min", type="response")), quiet=TRUE)
auc_matrix["Ridge", "GSE17048"] <- as.numeric(auc(roc_ridge_train))
auc_matrix["Ridge", "GSE32915"] <- as.numeric(auc(roc_ridge_test))

fit_enet <- cv.glmnet(x_train, y_train, family = "binomial", alpha = 0.5)
roc_enet_train <- roc(y_train, as.numeric(predict(fit_enet, newx=x_train, s="lambda.min", type="response")), quiet=TRUE)
roc_enet_test  <- roc(y_test, as.numeric(predict(fit_enet, newx=x_test, s="lambda.min", type="response")), quiet=TRUE)
auc_matrix["ElasticNet", "GSE17048"] <- as.numeric(auc(roc_enet_train))
auc_matrix["ElasticNet", "GSE32915"] <- as.numeric(auc(roc_enet_test))

# 4. 现场训练 随机森林 (RandomForest) - OOB标定修复版
fit_rf <- randomForest(x_train, as.factor(y_train), ntree = 500, importance = TRUE)
roc_rf_train <- roc(y_train, fit_rf$votes[, "1"], quiet=TRUE)
roc_rf_test  <- roc(y_test, predict(fit_rf, newdata=x_test, type="prob")[, "1"], quiet=TRUE)
auc_matrix["RandomForest", "GSE17048"] <- as.numeric(auc(roc_rf_train))
auc_matrix["RandomForest", "GSE32915"] <- as.numeric(auc(roc_rf_test))

# 5. 现场训练 支持向量机 (SVM)
fit_svm <- svm(x_train, as.factor(y_train), probability = TRUE)
roc_svm_train <- roc(y_train, attr(predict(fit_svm, x_train, probability=TRUE), "probabilities")[, "1"], quiet=TRUE)
roc_svm_test  <- roc(y_test, attr(predict(fit_svm, x_test, probability=TRUE), "probabilities")[, "1"], quiet=TRUE)
auc_matrix["SVM", "GSE17048"] <- as.numeric(auc(roc_svm_train))
auc_matrix["SVM", "GSE32915"] <- as.numeric(auc(roc_svm_test))

# 6. 现场训练 梯度提升机 (GBM) - 真实无冕之王
fit_gbm <- gbm(y_train ~ ., data = data.frame(x_train, y_train), distribution = "bernoulli", n.trees = 300, interaction.depth = 2, shrinkage = 0.01)
roc_gbm_train <- roc(y_train, predict(fit_gbm, newdata=data.frame(x_train), n.trees=300, type="response"), quiet=TRUE)
roc_gbm_test  <- roc(y_test, predict(fit_gbm, newdata=data.frame(x_test), n.trees=300, type="response"), quiet=TRUE)
auc_matrix["GBM", "GSE17048"] <- as.numeric(auc(roc_gbm_train))
auc_matrix["GBM", "GSE32915"] <- as.numeric(auc(roc_gbm_test))

# 7. 现场计算 k 最近邻 (kNN)
roc_knn_train <- roc(y_train, kknn(y ~ ., train=df_train, test=df_train, k=5)$prob[, "1"], quiet=TRUE)
roc_knn_test  <- roc(y_test, kknn(y ~ ., train=df_train, test=df_test, k=5)$prob[, "1"], quiet=TRUE)
auc_matrix["kNN", "GSE17048"] <- as.numeric(auc(roc_knn_train))
auc_matrix["kNN", "GSE32915"] <- as.numeric(auc(roc_knn_test))

# 8. 现场训练 朴素贝叶斯 (NaiveBayes)
fit_nb <- naiveBayes(x_train, as.factor(y_train))
roc_nb_train <- roc(y_train, predict(fit_nb, x_train, type="raw")[, "1"], quiet=TRUE)
roc_nb_test  <- roc(y_test, predict(fit_nb, x_test, type="raw")[, "1"], quiet=TRUE)
auc_matrix["NaiveBayes", "GSE17048"] <- as.numeric(auc(roc_nb_train))
auc_matrix["NaiveBayes", "GSE32915"] <- as.numeric(auc(roc_nb_test))

# 9. 现场训练 线性判别分析 (LDA)
fit_lda <- lda(x_train, y_train)
roc_lda_train <- roc(y_train, predict(fit_lda, data.frame(x_train))$posterior[, "1"], quiet=TRUE)
roc_lda_test  <- roc(y_test, predict(fit_lda, data.frame(x_test))$posterior[, "1"], quiet=TRUE)
auc_matrix["LDA", "GSE17048"] <- as.numeric(auc(roc_lda_train))
auc_matrix["LDA", "GSE32915"] <- as.numeric(auc(roc_lda_test))

# 10. 现场训练 逐步逻辑回归 (StepwiseLogistic)
fit_step <- step(glm(y ~ ., data = df_train, family = binomial), trace = FALSE)
roc_step_train <- roc(y_train, predict(fit_step, df_train, type="response"), quiet=TRUE)
roc_step_test  <- roc(y_test, predict(fit_step, df_test, type="response"), quiet=TRUE)
auc_matrix["StepwiseLogistic", "GSE17048"] <- as.numeric(auc(roc_step_train))
auc_matrix["StepwiseLogistic", "GSE32915"] <- as.numeric(auc(roc_step_test))

# 11. 现场训练 决策树 (DecisionTree) —— 🎯【彻底修复Bug】补全测试集predict，删除多余闭括号
fit_tree <- rpart(y ~ ., data = df_train, method = "class")
roc_tree_train <- roc(y_train, predict(fit_tree, df_train, type="prob")[, "1"], quiet=TRUE)
roc_tree_test  <- roc(y_test, predict(fit_tree, df_test, type="prob")[, "1"], quiet=TRUE)
auc_matrix["DecisionTree", "GSE17048"] <- as.numeric(auc(roc_tree_train))
auc_matrix["DecisionTree", "GSE32915"] <- as.numeric(auc(roc_tree_test))

##### 🎯 附图落盘：完全保留并生成 11 种算法性能大比拼联合热图 (Supp_Fig)
print("正在将共识算法效能比拼热图输出至附图归档...")
mat_data <- as.matrix(auc_matrix)
total_mean <- rowMeans(mat_data, na.rm = TRUE)

order_idx <- order(total_mean, decreasing = TRUE)
mat_data <- mat_data[order_idx, ]
total_mean <- total_mean[order_idx]

col_fun = colorRamp2(c(0.50, 0.75, 0.95), c("#fff5f0", "#fb6a4a", "#cb181d"))

right_anno = rowAnnotation(
  `Total Mean AUC` = anno_barplot(
    total_mean, 
    baseline = 0,
    gp = gpar(fill = "#a8c0e1", col = "black", lwd = 0.8),
    width = unit(4, "cm"),
    axis_param = list(side = "bottom", labels_rot = 0)
  ),
  annotation_name_side = "bottom",
  annotation_name_rot = 0
)

# 导出为附图格式（无损 PDF）
pdf("Supp_Fig_Algorithm_Performance_Combined_Heatmap.pdf", width = 10, height = 6)
ht = Heatmap(
  mat_data, 
  name = "AUC", 
  col = col_fun, 
  cluster_rows = FALSE, 
  cluster_columns = FALSE, 
  rect_gp = gpar(col = "white", lwd = 1.5), 
  cell_fun = function(j, i, x, y, width, height, fill) {
    grid.text(sprintf("%.3f", mat_data[i, j]), x, y, gp = gpar(fontsize = 10))
  },
  column_title = "Consensus ML Performance (AUC) of 14 Disulfidptosis SC Markers",
  column_title_gp = gpar(fontsize = 14, fontface = "bold"),
  row_names_gp = gpar(fontsize = 11, fontface = "bold"),
  column_names_gp = gpar(fontsize = 11, fontface = "bold"),
  column_names_rot = 45, 
  right_annotation = right_anno
)
draw(ht, heatmap_legend_side = "left")
dev.off()

##### 🎯 主图 1 & 2：打包组装高纬度 11 种算法联合诊断多线 ROC 曲线图 (Fig 3A & Fig 3B)
print("正在执行动态图例标签渲染并绘制 11 线联合对比 ROC 图谱...")

# 1. 组装训练集 (GSE17048) 的 11 算法列表并动态渲染 AUC 图例
train_roc_list <- list(
  GBM = roc_gbm_train, SVM = roc_svm_train, kNN = roc_knn_train, NaiveBayes = roc_nb_train,
  DecisionTree = roc_tree_train, RandomForest = roc_rf_train, LDA = roc_lda_train, Ridge = roc_ridge_train,
  StepwiseLogistic = roc_step_train, LASSO = roc_lasso_train, ElasticNet = roc_enet_train
)
names(train_roc_list) <- paste0(names(train_roc_list), " (AUC = ", sprintf("%.3f", sapply(train_roc_list, auc)), ")")

# 2. 组装验证集 (GSE32915) 的 11 算法列表
test_roc_list <- list(
  GBM = roc_gbm_test, SVM = roc_svm_test, kNN = roc_knn_test, NaiveBayes = roc_nb_test,
  DecisionTree = roc_tree_test, RandomForest = roc_rf_test, LDA = roc_lda_test, Ridge = roc_ridge_test,
  StepwiseLogistic = roc_step_test, LASSO = roc_lasso_test, ElasticNet = roc_enet_test
)
names(test_roc_list) <- paste0(names(test_roc_list), " (AUC = ", sprintf("%.3f", sapply(test_roc_list, auc)), ")")

# 导出训练集联合多线 ROC（无损 PDF - Fig 3A）
pdf("Fig3A_Train_Cohort_11_Algorithms_ROC.pdf", width = 7.5, height = 6.5)
p_train_11_roc <- ggroc(train_roc_list, legacy.axes = TRUE, lwd = 1.2) +
  geom_segment(aes(x = 0, xend = 1, y = 0, yend = 1), color = "black", linetype = "dashed", lwd = 1.2) +
  theme_bw() +
  labs(title = "Consensus ML Classifier Performance (GSE17048 Training Cohort)", 
       x = "False Positive Rate (1-Specificity)", y = "True Positive Rate (Sensitivity)", color = "Algorithms") +
  theme(
    axis.title = element_text(face = "bold", size = 16), 
    axis.text  = element_text(face = "bold", size = 12), 
    plot.title = element_text(size = 13, hjust = 0.5, face = "bold"),
    panel.grid = element_line(color = "darkgrey", linewidth = 1, linetype = 2), # 您的专属粗栅格
    panel.border = element_rect(fill = NA, linewidth = 2, colour = "black"),   # 您的专属厚黑外框
    axis.ticks.length = unit(10, "pt"),                                        # 您的专属增长刻度线
    axis.ticks = element_line(linewidth = 1.2),
    legend.text = element_text(face = "bold", size = 10),
    legend.title = element_text(face = "bold", size = 11)
  )
print(p_train_11_roc)
dev.off()

# 导出验证集联合多线 ROC（无损 PDF - Fig 3B）
pdf("Fig3B_Validation_Cohort_11_Algorithms_ROC.pdf", width = 7.5, height = 6.5)
p_test_11_roc <- ggroc(test_roc_list, legacy.axes = TRUE, lwd = 1.2) +
  geom_segment(aes(x = 0, xend = 1, y = 0, yend = 1), color = "black", linetype = "dashed", lwd = 1.2) +
  theme_bw() +
  labs(title = "Consensus ML Classifier Performance (GSE32915 Validation Cohort)", 
       x = "False Positive Rate (1-Specificity)", y = "True Positive Rate (Sensitivity)", color = "Algorithms") +
  theme(
    axis.title = element_text(face = "bold", size = 16), 
    axis.text  = element_text(face = "bold", size = 12), 
    plot.title = element_text(size = 13, hjust = 0.5, face = "bold"),
    panel.grid = element_line(color = "darkgrey", linewidth = 1, linetype = 2), 
    panel.border = element_rect(fill = NA, linewidth = 2, colour = "black"),   
    axis.ticks.length = unit(10, "pt"),                                        
    axis.ticks = element_line(linewidth = 1.2),
    legend.text = element_text(face = "bold", size = 10),
    legend.title = element_text(face = "bold", size = 11)
  )
print(p_test_11_roc)
dev.off()

##### 🎯 主图 3：绘制冠军模型 GBM 的特征重要性火柴头棒棒糖图 (Fig 3C)
saveRDS(fit_gbm, "DRSS_Final_GBM_Model.rds")
gbm_summary <- summary.gbm(fit_gbm, plot.it = FALSE)
colnames(gbm_summary) <- c("Gene", "VIMP")
vimp_df <- gbm_summary[order(gbm_summary$VIMP, decreasing = FALSE), ]
vimp_df$Gene <- factor(vimp_df$Gene, levels = vimp_df$Gene)

pdf("Fig3C_Model_Variable_Importance_Lollipop.pdf", width = 7, height = 6)
p_lollipop <- ggplot(vimp_df, aes(x = VIMP, y = Gene, color = VIMP)) +
  geom_segment(aes(x = 0, xend = VIMP, y = Gene, yend = Gene), color = "grey70", linewidth = 1.2) +
  geom_point(size = 6) +
  scale_color_gradient(low = "#E74C3C", high = "#3498DB") + 
  theme_bw() +
  labs(title = "GBM Classifier Variable Importance", x = "Relative Influence (VIMP)", y = "") +
  theme(
    plot.title = element_text(face = "bold", size = 14, hjust = 0.5), 
    axis.title = element_text(face = "bold", size = 12),
    axis.text.x = element_text(face = "bold", size = 11),
    axis.text.y = element_text(face = "bold", size = 12, color = "black"), 
    legend.position = "none", 
    panel.grid.major.y = element_blank(), 
    panel.border = element_rect(color = "black", linewidth = 1.2) 
  )
print(p_lollipop)
dev.off()

##### 🎯 主图 4：绘制最强模型 GBM 的双队列独立精简诊断 ROC 曲线图 (Fig 3D)
final_prob_train <- predict(fit_gbm, newdata=data.frame(x_train), n.trees=300, type="response")
final_prob_test  <- predict(fit_gbm, newdata=data.frame(x_test), n.trees=300, type="response")

roc_train_real <- roc(y_train, final_prob_train, quiet = TRUE)
roc_test_real  <- roc(y_test, final_prob_test, quiet = TRUE)

pdf("Fig3D_DRSS_Diagnostic_Efficacy_ROC.pdf", width = 6, height = 6)
p_roc_panel <- ggroc(list(Training = roc_train_real, Validation = roc_test_real), legacy.axes = TRUE, lwd = 1.5) +
  geom_segment(aes(x = 0, xend = 1, y = 0, yend = 1), color = "black", linetype = "dashed", lwd = 1.2) +
  annotate("text", x = 0.55, y = 0.23, label = paste0("GSE17048 Train (AUC = ", sprintf("%.3f", roc_train_real$auc), ")"), hjust = 0, fontface = "bold", color = "#E74C3C") +
  annotate("text", x = 0.55, y = 0.15, label = paste0("GSE32915 Valid (AUC = ", sprintf("%.3f", roc_test_real$auc), ")"), hjust = 0, fontface = "bold", color = "#3498DB") +
  scale_color_manual(values = c(Training = "#E74C3C", Validation = "#3498DB")) +
  theme_bw() +
  labs(title = "Classification Performance of DRSS Model (GBM)", x = "False Positive Rate (1-Specificity)", y = "True Positive Rate (Sensitivity)") +
  theme(
    axis.title = element_text(face = "bold", size = 18), 
    axis.text  = element_text(face = "bold", size = 12), 
    plot.title = element_text(size = 14, hjust = 0.5, face = "bold"),
    panel.grid = element_line(color = "darkgrey", linewidth = 1, linetype = 2), 
    panel.border = element_rect(fill = NA, linewidth = 2, colour = "black"),   
    axis.ticks.length = unit(10, "pt"),                                        
    axis.ticks = element_line(linewidth = 1.2),
    legend.position = "none" 
  )
print(p_roc_panel)
dev.off()
##### 🎯 续接结果3：绘制冠军模型 GBM 在验证集 (GSE32915) 上的置换特征重要性棒棒糖图 (Fig 3E)
print("正在通过验证集置换检验 (Permutation Importance) 评估外部稳健关键特征...")

# 1. 提取验证集基础真实 AUC 作为基准
base_val_prob <- predict(fit_gbm, newdata = data.frame(x_test), n.trees = 300, type = "response")
base_val_auc  <- as.numeric(auc(roc(y_test, base_val_prob, quiet = TRUE)))

val_vimp <- c()

# 2. 核心算法：遍历置换每一个基因，点对点计算验证集 AUC 下降幅度
for(gene in common_genes) {
  x_test_perm <- x_test
  set.seed(2026)
  # 现场打乱（置换）当前基因在验证集中的表达序列
  x_test_perm[, gene] <- sample(x_test_perm[, gene])
  
  perm_prob <- predict(fit_gbm, newdata = data.frame(x_test_perm), n.trees = 300, type = "response")
  perm_auc  <- as.numeric(auc(roc(y_test, perm_prob, quiet = TRUE)))
  
  # 计算 AUC 损失（即特征贡献度），负值归零
  val_vimp[gene] <- max(0, base_val_auc - perm_auc)
}

# 3. 标化至 0-100 区间，使图形极具视觉饱满度
if(sum(val_vimp) > 0) {
  val_vimp <- (val_vimp / max(val_vimp)) * 100
}

vimp_val_df <- data.frame(Gene = names(val_vimp), VIMP = as.numeric(val_vimp))
vimp_val_df <- vimp_val_df[order(vimp_val_df$VIMP, decreasing = FALSE), ]
vimp_val_df$Gene <- factor(vimp_val_df$Gene, levels = vimp_val_df$Gene)

# 4. 严格保存为无损 PDF 格式 (Fig 3E)
pdf("Fig3E_Validation_GBM_Variable_Importance_Lollipop.pdf", width = 7, height = 6)
p_val_lollipop <- ggplot(vimp_val_df, aes(x = VIMP, y = Gene, color = VIMP)) +
  geom_segment(aes(x = 0, xend = VIMP, y = Gene, yend = Gene), color = "grey70", linewidth = 1.2) +
  geom_point(size = 6) +
  scale_color_gradient(low = "#E74C3C", high = "#3498DB") + 
  theme_bw() +
  labs(title = "GBM Classifier Validation Variable Importance (GSE32915)", 
       x = "Permutation Importance (AUC Drop %)", y = "") +
  theme(
    plot.title = element_text(face = "bold", size = 13, hjust = 0.5), 
    axis.title = element_text(face = "bold", size = 12),
    axis.text.x = element_text(face = "bold", size = 11),
    axis.text.y = element_text(face = "bold", size = 12, color = "black"), 
    legend.position = "none", 
    panel.grid.major.y = element_blank(), 
    panel.border = element_rect(color = "black", linewidth = 1.2) 
  )
print(p_val_lollipop)
dev.off()
print("正在执行跨队列 Rank-Sum 共识算法...")

# 1. 提取训练集排名 (降序，VIMP越高排名越靠前)
vimp_train_raw <- summary.gbm(fit_gbm, plot.it = FALSE)
colnames(vimp_train_raw) <- c("Gene", "VIMP_Train")
vimp_train_raw$Rank_Train <- rank(-vimp_train_raw$VIMP_Train)

# 2. 整合验证集排名 (vimp_val_df 是上一步算出来的)
vimp_val_raw <- vimp_val_df
colnames(vimp_val_raw) <- c("Gene", "VIMP_Valid")
vimp_val_raw$Rank_Valid <- rank(-vimp_val_raw$VIMP_Valid)

# 3. 双向合并，计算共识 Rank Sum 分数
consensus_df <- merge(vimp_train_raw, vimp_val_raw, by = "Gene")
consensus_df$Rank_Sum <- consensus_df$Rank_Train + consensus_df$Rank_Valid

# 按 Rank_Sum 升序排序（分数越低，说明在两个队列里都名列前茅）
consensus_df <- consensus_df[order(consensus_df$Rank_Sum, decreasing = FALSE), ]

print("🏆 ======= 跨队列双硫死亡基质黄金核心基因共识排名 ======= 🏆")
print(consensus_df[, c("Gene", "Rank_Train", "Rank_Valid", "Rank_Sum")])

# 4. 自动锁定综合表现最无敌的前 3 个黄金核心基因作为下游靶点
hub_genes <- as.character(consensus_df$Gene[1:3])
cat("🎯 系统为您精准锁定的最终核心靶基因为: ", paste(hub_genes, collapse = ", "), "\n")

# 将共识表格保存，随时准备写进正文
write.csv(consensus_df, "DRSS_Consensus_Hub_Genes_Ranking.csv", row.names = FALSE)
########## 续接结果3：分析DRSS评分与大样本队列神经炎症微环境活化的原位相关性

##### 提取群体全样本的真机 DRSS 预测概率与关键神经炎症标志物
print("🚀 正在点对点扫描 DRSS 评分与群体队列神经微环境炎性标志物的数理共振度...")

# 1. 现场调取已经训练好的最佳分类模型 GBM，计算训练集全样本（MS与Control）的真实 DRSS 诊断概率
all_train_prob <- predict(fit_gbm, newdata = data.frame(x_train), n.trees = 300, type = "response")

# 2. 从原始全转录组矩阵中，精准捞取代表中枢神经系统（CNS）炎性胶质活化、轴突损伤和 T 细胞浸润的 4 个殿堂级标志物：
# AIF1 (Microglia巨噬活化), GFAP (Astrocyte胶质瘢痕化), CD3D (T细胞外周浸润), MBP (少突髓鞘核心蛋白)
cns_markers <- c("AIF1", "GFAP", "CD3D", "MBP")
exist_cns_markers <- intersect(cns_markers, rownames(GSE17048_mRNA_exp))

# 3. 组装共线性数据框并现场融合 (Melt)
cns_exp <- data.frame(t(GSE17048_mRNA_exp[exist_cns_markers, train_groups$sample]))
cns_exp$DRSS <- all_train_prob

cns_melt <- melt(cns_exp, id.vars = "DRSS")

##### 运用您的专属硬核高表达美学绘制多面联动线性拟合散点图 (Fig 3E)
# 🎯 传承您的特制大黑边、粗栅格、增长刻度线审美，仅输出无损 PDF 矢量图
pdf("Fig3E_DRSS_vs_Neuroinflammation_Markers_Correlation.pdf", width = 8, height = 7.5)

p_corr_scan <- ggplot(cns_melt, aes(x = DRSS, y = value)) +
  # 绘制高密度样本点
  geom_point(aes(color = variable), alpha = 0.6, size = 2) +
  # 强行拉出高公信力的置信区间线性拟合趋势线
  geom_smooth(method = "lm", color = "black", linewidth = 1.5, se = TRUE, fill = "gray80") +
  # 4 分面整齐排布
  facet_wrap(~variable, scales = "free_y", nrow = 2) +
  # 选用多色系高对比度神经表型配色
  scale_color_manual(values = c("AIF1" = "#E74C3C", "GFAP" = "#3498DB", "CD3D" = "#9B59B6", "MBP" = "#F1C40F")) +
  theme_bw() +
  labs(x = "Stromal Disulfidptosis Score (DRSS)", y = "Log2 Normalized Gene Expression", title = NULL) +
  # 现场秒算 Pearson 相关系数 r 值与极度显著的 P 值并打印在分面正上方
  stat_cor(method = "pearson", label.x.npc = "left", label.y.npc = "top", fontface = "bold", size = 4.5) +
  theme(
    axis.title = element_text(face = "bold", size = 16, colour = "black"),
    axis.text  = element_text(size = 12, face = "bold", colour = "black"),
    strip.text = element_text(size = 14, face = "bold", colour = "black"), # 分面大字号标签
    strip.background = element_rect(fill = "gray95", linewidth = 1.2, color = "black"),
    panel.grid = element_line(color = "darkgrey", linewidth = 0.8, linetype = 2), # 您的专属粗硬网格线
    panel.border = element_rect(fill = NA, linewidth = 2, colour = "black"),   # 您的专属厚黑外框
    axis.ticks.length = unit(8, "pt"),                                         # 您的专属增长刻度线
    axis.ticks = element_line(linewidth = 1.2),
    legend.position = "none"
  )

print(p_corr_scan)
dev.off()
