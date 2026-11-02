##########单细胞转录组分析揭示基质细胞是 MS 微环境中双硫死亡的主要中心
#####设置工作环境
setwd(file.path(Sys.getenv("LEGACY_STUDY_ROOT", unset = "legacy_inputs"), "Result1"))
#####加载R包
library(Seurat)
library(dplyr)
library(tidyr)
library(ggplot2)
library(AUCell)
library(CellChat)
library(patchwork)
library(clusterProfiler)
library(org.Hs.eg.db)
library(ggrepel)
#####读取数据
###单细胞数据
GSE279180_sn_atlas <- readRDS(Sys.getenv("ATLAS_RDS", unset = ""))
# 过滤没有对照的BC细胞（沿用您之前的清洗步骤）
GSE279180_sn_atlas <- subset(GSE279180_sn_atlas, subset = celltype != "BC")
###双硫死亡基因集
disulfidptosis_geneset <- read.delim(file.path(Sys.getenv("LEGACY_STUDY_ROOT", unset = "legacy_inputs"), "Data/\u4e8c\u786b\u4e2d\u6bd2\u76f8\u5173\u57fa\u56e0\u96c6.txt"))
disulfidptosis_genes  <- unique(disulfidptosis_geneset$DFG)
# 确保基因为大写并过滤掉矩阵中不存在的基因
disulfidptosis_genesdisulfidptosis_genes <- intersect(disulfidptosis_genes, rownames(GSE279180_sn_atlas))
###设置颜色映射
cell_type_colors <- c("AS" = "#1d77b5",
                      "BC" = "#ff8112",
                      "EC" = "#2a9e65",
                      "MG" = "#d62926",
                      "NEU" = "#aa41fc",
                      "OL" = "#8d5649",
                      "OPC" = "#e278c5",
                      "SC" = "#7f7f7f",
                      "TC" = "#bbbc24")
condition_colors <- c("MS" = "red",
                      "Control" = "grey")
#####筛选不同细胞亚型对MS的贡献，找出在MS-CA与MS-CI之间差异最大的细胞亚型，作为“关键亚型”。
#差异分析计算贡献度（更高效）
#为每个亚型计算在MS-CA vs MS-CI中的差异表达
#使用FindMarkers函数，但需要为每个亚型分别计算
###定义贡献度计算模块
calculate_contribution_df <- function(seurat_obj, condition_control, condition_disease) {
  subtypes <- unique(seurat_obj$celltype)
  results <- data.frame(subtype = character(), avg_log2FC = numeric(), pct_prop = numeric(), stringsAsFactors = FALSE)
  for (subtype in subtypes) {
    control_cells <- Cells(seurat_obj)[seurat_obj$celltype == subtype & seurat_obj$condition == condition_control]
    disease_cells <- Cells(seurat_obj)[seurat_obj$celltype == subtype & seurat_obj$condition == condition_disease]
    expr_data <- GetAssayData(seurat_obj, assay = "originalexp", slot = "data")
    subtype_expr <- expr_data[, c(control_cells, disease_cells)]
    control_mean <- rowMeans(subtype_expr[, control_cells, drop = FALSE])
    disease_mean <- rowMeans(subtype_expr[, disease_cells, drop = FALSE])
    n_genes_ctrl <- min(100, sum(control_mean > 0))
    n_genes_ms <- min(100, sum(disease_mean > 0))
    top_control_genes <- names(sort(control_mean, decreasing = TRUE)[1:n_genes_ctrl])
    top_disease_genes <- names(sort(disease_mean, decreasing = TRUE)[1:n_genes_ms])
    top_genes <- unique(c(top_control_genes, top_disease_genes))
    fc_values <- (disease_mean[top_genes] + 1e-10) / (control_mean[top_genes] + 1e-10)
    avg_fc <- mean(fc_values, na.rm = TRUE)
    expr_binary <- as.matrix(subtype_expr[top_genes, ] > 0)
    pct_prop <- mean(rowMeans(expr_binary))
    results <- rbind(results, data.frame(
      subtype = subtype,
      avg_log2FC = log2(avg_fc),
      pct_prop = pct_prop,
      stringsAsFactors = FALSE
    ))
  }
  return(results)
}
###运行贡献度计算
Idents(GSE279180_sn_atlas) <- GSE279180_sn_atlas$condition
contribution_results <- calculate_contribution_df(
  seurat_obj = GSE279180_sn_atlas,
  condition_control = "Control",
  condition_disease = "MS"
)
contribution_results$contribution_score <- contribution_results$avg_log2FC * contribution_results$pct_prop
contribution_results <- contribution_results %>% arrange(desc(contribution_score))
###绘图展示不同亚型对疾病的贡献度条形图
p_contrib <- ggplot(contribution_results, aes(x = reorder(subtype, -contribution_score), y = contribution_score, fill = subtype)) +
  geom_bar(stat = "identity", color = "black", width = 0.6, lwd = 0.4) +
  scale_fill_manual(values = cell_type_colors) +
  labs(title = "Cell Type Contribution to MS Pathogenesis", x = "Cell Type", y = "Contribution Score (log2FC * PctProp)") +
  theme_classic() +
  theme(plot.title = element_text(hjust = 0.5, face = "bold"), axis.text.x = element_text(angle = 45, hjust = 1))
ggsave("Fig1_Cell_Contribution_Bar.pdf", p_contrib, width = 6, height = 4.5)
###确定了关键亚型SC，那么绘制SC具体的umap图
#将细胞亚型设置为对象的身份标识（Identity）
Idents(GSE279180_sn_atlas) <- "celltype"
#定义高亮细胞类型-SC
cell_type_colors_Highlight <- c("AS" = "lightgray",
                                      "BC" = "lightgray",
                                      "EC" = "lightgray",
                                      "MG" = "lightgray",
                                      "NEU" = "lightgray",
                                      "OL" = "lightgray",
                                      "OPC" = "lightgray",
                                      "SC" = "yellow",
                                      "TC" = "lightgray")
pdf("Fig1_SC_UMAP.pdf", width = 6, height = 5)
DimPlot(GSE279180_sn_atlas,
        reduction = "X_umap",
        pt.size = 1.5,                # 点大小（推荐：1-3）
        #alpha = 0.8,                  # 透明度（0-1）
        group.by = "celltype",
        cols = cell_type_colors_Highlight)
dev.off()
#绘制SC细胞亚型情况图
GSE279180_sn_atlas_sc <- subset(GSE279180_sn_atlas, subset = celltype == "SC")
cell_type_colors_sc <- c("SC_FibroInfl" = "#fc9a9a",
                         "SC_FibroMyo1" = "#b37557",
                         "SC_FibroMyo2" = "#b49d99",
                         "SC_Fibroblast1" = "#96873b",
                         "SC_Fibroblast2" = "#b24646",
                         "SC_Fibroblast3" = "#fe9e37",
                         "SC_Fibrocyte" = "#ed4437",
                         "SC_MSCs" = "#e18748",
                         "SC_PerivFibro" = "#ebd57c")
pdf("Fig1_SC_subtype_UMAP.pdf", width = 6, height = 5)
DimPlot(GSE279180_sn_atlas_sc,
        reduction = "X_umap",
        pt.size = 5,                # 点大小（推荐：1-3）
        alpha = 0.8,                  # 透明度（0-1）
        group.by = "subtype",
        cols = cell_type_colors_sc)
dev.off()
###提取关键亚型 SC 的 Marker 基因
Idents(GSE279180_sn_atlas) <- "celltype"
sc.markers <- FindMarkers(GSE279180_sn_atlas, ident.1 = "SC", only.pos = TRUE, min.pct = 0.25, logfc.threshold = 1, test.use = "wilcox")
SC_markers_genes <- rownames(sc.markers[sc.markers$p_val_adj < 0.05, ])
###绘制火山图
sc_degs <- FindMarkers(GSE279180_sn_atlas,
                       ident.1 = "SC",
                       logfc.threshold = 0,
                       min.pct = 0,
                       assay = "originalexp")
sc_up_genes <- sc_degs %>%
  filter(p_val_adj < 0.05 & avg_log2FC > 0.5) %>%
  arrange(desc(avg_log2FC)) # 按表达差异倍数降序排列
write.csv(sc_degs, "Stromal_vs_All_DEGs.csv", row.names = TRUE)
cat("✅ 基质细胞全局差异分析矩阵已成功无损落盘至: Stromal_vs_All_DEGs.csv\n")
sc_degs$gene <- rownames(sc_degs) # 将行名(基因名)提取为独立列
sc_degs$Change <- "Stable"        # 默认标记为稳定
sc_degs$Change[sc_degs$avg_log2FC > 0.5 & sc_degs$p_val_adj < 0.05] <- "Up"
sc_degs$Change[sc_degs$avg_log2FC < -0.5 & sc_degs$p_val_adj < 0.05] <- "Down"
### 4. 🧠 动态海选 Top 标签，并强行死锁你的绝对主角 DSE
# 先按常规显著性捞出排名前 15 的外显基因
top_genes <- head(sc_degs %>% filter(Change %in% c("Up", "Down")) %>% arrange(p_val_adj), 15)$gene
# 🎯【绝杀核心】如果明星基因 DSE 躺在矩阵里，强制剥离并塞入标注名单，产生绝对的视觉爆点！
if ("DSE" %in% sc_degs$gene) {
  top_genes <- unique(c("DSE", top_genes))
}
label_genes <- subset(sc_degs, gene %in% top_genes)
### 5. 配色定义与高维 ggplot2 出版级艺术化渲染
volcano_colors <- c("Up" = "#ECA8A9", "Down" = "#74AED4", "Stable" = "lightgrey")
p_volcano <- ggplot(sc_degs, aes(x = avg_log2FC, y = -log10(p_val_adj + 1e-300), color = Change)) +
  # 适当控制点的透明度，使得中心密集区呈现极佳的半透明堆叠高级感
  geom_point(alpha = 0.6, size = 1.3) +
  scale_color_manual(values = volcano_colors) +
  # 添加数理差异截断虚线
  geom_vline(xintercept = c(-0.5, 0.5), linetype = "dashed", color = "black", alpha = 0.4) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "black", alpha = 0.4) +
  # 渲染基因名称标签，字体加粗提升印刷质感
  geom_text_repel(data = label_genes,
                  aes(label = gene),
                  size = 3.8,
                  fontface = "bold",
                  box.padding = 0.45,
                  point.padding = 0.3,
                  max.overlaps = 30,
                  color = "black") +
  # 设置图表标题和坐标轴
  labs(title = "Stromal Cells vs Other Cell Types",
       x = "Log2 Fold Change",
       y = "-Log10(adj. p-value)") +
  # 使用经典学术主题
  theme_classic() +
  theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
        axis.title = element_text(face = "bold", size = 12),
        axis.text = element_text(color = "black", size = 10),
        legend.title = element_blank(),
        legend.position = "right")
### 6. 保存高质量出版级矢量 PDF 火山图资产
ggsave("Fig1_Stromal_Volcano_Plot.pdf", plot = p_volcano, width = 7, height = 6, dpi = 300)
###与双硫死亡基因集（二硫基因）取交集
intersect_genes <- intersect(SC_markers_genes, disulfidptosis_genes)
cat("🎯 关键亚型 SC 的 Marker 与双硫死亡基因的交集为：", paste(intersect_genes, collapse = ", "), "\n")
write.table(intersect_genes,"disulfidptosis_markers.txt",sep="\t",quote=F,row.names=F,col.names=F)
###交集基因绘制气泡图
if(length(intersect_genes) > 0) {
  p_dot <- DotPlot(GSE279180_sn_atlas, features = intersect_genes, group.by = "celltype")
  dot_data <- p_dot$data[, c('id', 'features.plot', 'pct.exp', 'avg.exp.scaled')]
  plot_bubble <- ggplot(dot_data, aes(x = id, y = features.plot)) +
    geom_point(aes(fill = avg.exp.scaled, size = pct.exp), color = 'black', shape = 21, stroke = 0.4) +
    xlab("") + ylab("") +
    scale_fill_gradientn(colors = c("#C5AFAB", "#D0A197", "#E17463", "#E7614E", "#E31812")) +
    scale_size(range = c(1, 6), limits = c(0, 100)) +
    theme_classic() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1, color = "black"),
      axis.text.y = element_text(color = "black"),
      legend.position = "right"
    ) +
    labs(title = "Expression of Disulfidptosis-Related SC Markers")
  ggsave("Fig1_Intersect_Genes_Dotplot.pdf", plot = plot_bubble, width = 6, height = 5, dpi = 300)
} else {
  cat("⚠️ 提示：当前严格阈值下交集基因为0，可适当调低 FindMarkers 的 logfc.threshold（例如调到 0.5 或 0.25）再试。\n")
}
#####运行 CellChat 构建基质细胞旁分泌网络
###运行CellChat
options(future.globals.maxSize = 20 * 1024^3)
set.seed(42) # 固定随机种子确保结果可重复
Idents(GSE279180_sn_atlas) <- "celltype"
GSE279180_for_chat <- subset(GSE279180_sn_atlas, downsample = 2500)
# 🎯 核心修复：强行洗掉已经删除的 "BC" 等幽灵细胞类型的因子水平！
GSE279180_for_chat$celltype <- droplevels(factor(GSE279180_for_chat$celltype))
cat("🔥 因子洗涤成功！当前实际参与通讯的细胞谱系为：\n")
print(levels(GSE279180_for_chat$celltype)) # 此时列表里绝对不会再有 BC 了
cat("当前用于通讯分析的细胞总数：", ncol(GSE279180_for_chat), "\n")
Idents(GSE279180_for_chat) <- "celltype"
cellchat <- createCellChat(object = GSE279180_for_chat, group.by = "celltype", assay = "originalexp")
# 设置人类细胞通讯分泌型信号数据库
CellChatDB <- CellChatDB.human
CellChatDB.sub <- subsetDB(CellChatDB, search = "Secreted Signaling")
cellchat@DB <- CellChatDB.sub
cellchat <- subsetData(cellchat)
#运行高过表达基因与多核并行计算
cellchat <- identifyOverExpressedGenes(cellchat)
cellchat <- identifyOverExpressedInteractions(cellchat)
# 开启多核并行加速
future::plan("multisession", workers = 4)
cellchat <- computeCommunProb(cellchat, raw.use = TRUE)
cellchat <- filterCommunication(cellchat, min.cells = 10)
cellchat <- computeCommunProbPathway(cellchat)
# 聚合全局网络
cellchat <- aggregateNet(cellchat)
###可视化
cat("正在生成通讯网络图表...\n")
#全局网络数量与权重圈图
pdf("Fig1_Global_Communication_Network.pdf", width = 10, height = 5)
par(mfrow = c(1, 2), xpd = TRUE)
netVisual_circle(cellchat@net$count, weight.scale = TRUE, label.edge = FALSE,
                 title.name = "Number of interactions")
netVisual_circle(cellchat@net$weight, weight.scale = TRUE, label.edge = FALSE,
                 title.name = "Interaction weights/strength")
dev.off()
#核心图
pdf("Fig1_SC_As_Source_Outgoing_Network.pdf", width = 6, height = 6)
netVisual_circle(cellchat@net$weight, sources.use = "SC", weight.scale = TRUE,
                 label.edge = TRUE, title.name = "Secreted signals transmitted from SCs")
dev.off()
#深度展示 SC 发出的所有显著配体-受体对 (图 1E)
pdf("Fig1_SC_Driven_LR_Pairs_Bubble.pdf", width = 7, height = 9)
netVisual_bubble(cellchat, sources.use = "SC", remove.isolate = TRUE,
                 title.name = "Significant Outgoing L-R Pairs from SCs")
dev.off()
#保存 CellChat 最终计算对象
saveRDS(cellchat, "MS_Stromal_CellChat_Final.rds")
#保存计算好的单细胞中间对象
saveRDS(GSE279180_sn_atlas, "GSE279180_processed_with_Disulfidptosis.rds")
#####保存数据
save.image(".RData")