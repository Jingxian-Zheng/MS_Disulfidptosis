########## 进行Step3-下游机制与生物信息学分析 (14个核心双硫基质基因纯化版)
##### 设置工作环境
setwd(file.path(Sys.getenv("LEGACY_STUDY_ROOT", unset = "legacy_inputs"), "Result1"))

##### 加载R包
library(Seurat)
library(ggplot2)
library(dplyr)
library(corrplot)
library(pheatmap)
library(viridis)
library(ggpubr)

######## 1. 读取数据 (完美咬合结果1的处理后单细胞资产)
GSE279180_sn_atlas <- readRDS(Sys.getenv("ATLAS_RDS", unset = ""))
######## 🎯【核心替换】锁定最终确定的 14 个基质特异性双硫死亡核心特征基因
key_genes <- c("FMO2", "SLC6A20", "PTGIS", "BHMT2", "HSD11B1", "DSE", "OSBPL5",
               "ATP10D", "IRS1", "CALU", "SEC61A1", "RRBP1", "PRKAA1", "MTAP")

# 安全检查：确保所有挑选的基因都老老实实存在于当前的表达矩阵中
key_genes <- intersect(key_genes, rownames(GSE279180_sn_atlas))

######## 2. 核心基因在单细胞中的表达水平 (FeaturePlot)
print("🚀 正在绘制 14 个黄金靶点在单细胞流形空间的空间分布图 (FeaturePlot)...")
# 使用高容错率循环批量生成 PDF，拒绝复制粘贴14遍，排版更优雅
pdf("Fig1_Key_Genes_FeaturePlot_Summary.pdf", width = 10, height = 8)
for (gene in key_genes) {
  p_feat <- FeaturePlot(GSE279180_sn_atlas, features = gene, reduction = "X_umap",
                        pt.size = 1.5, order = TRUE, cols = c("grey90", "#E31812")) +
            theme(plot.title = element_text(face = "bold", size = 14))
  print(p_feat)
}
dev.off()
print("🚀 正在启动流形色彩拉伸算法，强行熔断 1.00 对角线并放大真实相关性信号...")

# 提取 14 个核心基因的表达矩阵
expr_matrix <- FetchData(GSE279180_sn_atlas, vars = key_genes, slot = "data")

# 计算基因之间的 Spearman 相关性
cor_matrix <- cor(expr_matrix, method = "spearman")

# 🎯【核心修改一】抹去对角线自我对比的 1，释放色彩天花板
diag(cor_matrix) <- NA

# 🎯【核心修改二】动态捕捉真实数据的边界值，实施像素级色彩拉伸
min_true <- min(cor_matrix, na.rm = TRUE)
max_true <- max(cor_matrix, na.rm = TRUE)

# 自动生成 101 个高精细色阶断点，强行把全谱红色压进 [min_true, max_true] 真实区间内
breaks_list <- seq(min_true, max_true, length.out = 101)

# 导出无损 PDF 发表级热图
pdf("Fig1_Key_Genes_Correlation_Heatmap.pdf", width = 7, height = 6.5)

pheatmap(cor_matrix,
         # 沿用高雅红系渐变
         col = colorRampPalette(c("#fff5f5", "#fca5a5", "#b2182b"))(100),
         breaks = breaks_list,         # 🎯 注入拉伸断点，彻底激活弱信号
         display_numbers = TRUE,
         number_format = "%.2f",
         fontsize_number = 9,
         number_color = "black",
         na_col = "white",             # 🎯 对角线留白！配合下面的白边框，视觉效果极其震撼
         border_color = "white",
         main = "Correlation Heatmap of 14 Core Genes (Stretched)")

dev.off()

cat("📊 绝美单色系相关性热图已成功洗涤并落盘为 Fig1_Key_Genes_Correlation_Heatmap.pdf！\n")
######## 5. 🎯【图层置顶+点缩小完美版】关键基因表达与细胞因子评分的差异
print("🚀 正在启动图层强行置顶与精细化像素渲染，全面复活 High 组红色信号...")

library(patchwork)
library(Seurat)
library(ggplot2)
library(dplyr)
library(ggpubr)

# 载入并整理细胞因子元数据基因集
cytokine_genes <- read.csv(file.path(Sys.getenv("LEGACY_STUDY_ROOT", unset = "legacy_inputs"), "Data/Cytokine_inf.csv"), sep = ',')
cytokine_genes <- unique(as.character(cytokine_genes$Human_gene_symbol))
cytokine_genes <- intersect(cytokine_genes, rownames(GSE279180_sn_atlas))

# 计算整个单细胞库里每个细胞反应性的 Cytokine 评分
GSE279180_sn_atlas <- AddModuleScore(GSE279180_sn_atlas,
                                     features = list(cytokine_genes),
                                     name = "Cytokine_Score")

# 创建独立的文件夹存放 14 个基因的分组精细对比 PDF
dir.create("./Cytokine_Diff_Plots", showWarnings = FALSE)

for (gene in key_genes) {
  cat("正在执行图层重排、计算显著性并拼图导出基因:", gene, "...\n")

  gene_vals <- FetchData(GSE279180_sn_atlas, vars = gene)[, 1]
  med_val <- median(gene_vals)

  # 如果某个基因表达0太多，改用平均值分层
  if(med_val == 0) { med_val <- mean(gene_vals) }

  GSE279180_sn_atlas$gene_group <- ifelse(gene_vals > med_val, "High", "Low")
  GSE279180_sn_atlas$gene_group <- factor(GSE279180_sn_atlas$gene_group, levels = c("Low", "High"))

  # 🎯【核心修复一】现场强行对十万个细胞重新排队：Low 在前（先画），High 在后（后画盖在最上面）
  cells_low <- colnames(GSE279180_sn_atlas)[GSE279180_sn_atlas$gene_group == "Low"]
  cells_high <- colnames(GSE279180_sn_atlas)[GSE279180_sn_atlas$gene_group == "High"]
  ordered_cells <- c(cells_low, cells_high)

  # UMAP 高低分层映射景观
  p1 <- DimPlot(GSE279180_sn_atlas,
                reduction = "X_umap",
                cells = ordered_cells,       # 🎯 强制注入重排后的条形码，实现红色绝对置顶
                pt.size =1,               # 🎯 核心修复二：点大小从 1 缩减到 0.4，消除大面积颜色污染
                alpha = 0.8,
                group.by = "gene_group",
                cols = c("Low" = "#1CB4B8", "High" = "#EB7369")) +
        labs(title = paste(gene, "Expression Stratification")) +
        theme_classic() + # 换成最干净的出版级经典主题
        theme(
          plot.title = element_text(face = "bold", size = 14, hjust = 0.5, color = "black"),
          axis.title = element_text(face = "bold", size = 11, color = "black"),
          axis.text = element_text(face = "bold", size = 9, color = "black"),
          axis.line = element_line(linewidth = 0.8, color = "black"),
          legend.text = element_text(face = "bold", size = 10),
          legend.position = "top" # 将图例移到顶部方便排版
        )

  # 细胞因子得分对比箱线图
  p2 <- ggplot(GSE279180_sn_atlas@meta.data, aes(x = gene_group, y = Cytokine_Score1, fill = gene_group, color = gene_group)) +
    stat_boxplot(geom = "errorbar", width = 0.2) +
    geom_boxplot(alpha = 0.6, lwd = 1.2, width = 0.5, outlier.shape = NA) +

    stat_compare_means(method = "wilcox.test",
                       label = "p.format",
                       label.x = 1.5,
                       label.y = 0.32,
                       size = 4.5,
                       fontface = "bold",
                       color = "black",
                       hjust = 0.5) +

    scale_y_continuous(breaks = seq(-0.5, 0.5, 0.1)) +
    coord_cartesian(ylim = c(-0.4, 0.4)) +
    scale_fill_manual(values = c("Low" = "#1CB4B8", "High" = "#EB7369")) +
    scale_color_manual(values = c("Low" = "#1CB4B8", "High" = "#EB7369")) +
    labs(title = paste("Cytokine Activity Response by", gene), x = paste(gene, "Group"), y = "Cytokine Score") +
    theme_bw() +
    theme(
      axis.title = element_text(colour = "black", face = "bold", size = 14),
      axis.text  = element_text(colour = "black", face = "bold", size = 11),
      plot.title = element_text(colour = "black", face = "bold", size = 13, hjust = 0.5),
      panel.background = element_blank(),
      panel.border = element_rect(fill = NA, linewidth = 1.8, colour = "black"),
      axis.title.x = element_text(margin = margin(0.3, 0, 0, 0, 'cm')),
      axis.title.y = element_text(margin = margin(0, 0.3, 0, 0, 'cm')),
      axis.ticks.length = unit(6, "pt"),
      axis.ticks = element_line(linewidth = 1.2),
      legend.position = "none"
    )

  # 左右拼图完美落盘
  combined_p <- p1 + p2 + plot_layout(ncol = 2)
  ggsave(filename = paste0("./Cytokine_Diff_Plots/Fig1_Cytokine_Response_via_", gene, ".pdf"),
         plot = combined_p, width = 11, height = 5, dpi = 300)
}

print("🏆 大功告成！红色信号已全部苏醒，组合矢量 PDF 已经完美刷新落盘！")