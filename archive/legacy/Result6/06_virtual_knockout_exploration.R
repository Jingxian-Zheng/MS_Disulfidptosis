###################### Step3-Results 6: 黄金靶点 DSE 的原位虚拟扰动实验与中枢神经微环境精准靶向
########## 设置工作环境
setwd(file.path(Sys.getenv("LEGACY_STUDY_ROOT", unset = "legacy_inputs"), "Result6"))
########## 加载R包
library(Seurat)
library(scTenifoldKnk)
library(ggplot2)
library(dplyr)
library(ggrepel)
########## 1. 读取多发性硬化症单细胞资产
GSE279180_sn_atlas <- readRDS(Sys.getenv("ATLAS_RDS", unset = ""))
# 严格执行结果 1 的清洗步骤：过滤没有对照的 BC 细胞
GSE279180_sn_atlas <- subset(GSE279180_sn_atlas, subset = celltype != "BC")
# 🎯【核心对齐】将身份标识锁定为结果 1 的 celltype，并切片提取绝对主角：基质细胞 (SC)
Idents(GSE279180_sn_atlas) <- "celltype"
target_cell_cluster <- "SC"
stromal_cells <- subset(GSE279180_sn_atlas, subset = celltype == target_cell_cluster)
##### 寻找基质高变特征基因，并执行目标靶点防漏安全盾
# 🎯【核心对齐】显式指定使用结果 1 的 originalexp 默认 assay 槽位，防止读取报错
DefaultAssay(stromal_cells) <- "originalexp"
stromal_cells <- FindVariableFeatures(stromal_cells, selection.method = "vst", nfeatures = 5000)
hvg_genes <- VariableFeatures(stromal_cells)
if (!"DSE" %in% hvg_genes) {
  hvg_genes <- c(hvg_genes, "DSE")
  cat("🎯 安全拦截：已将核心黄金基质靶点 DSE 强制逆向补回特征基因列表。\n")
}
# 提取未经 Scale 的原始 Counts 稠密矩阵准备网络流形对齐
counts_matrix <- as.matrix(GetAssayData(stromal_cells, assay = "originalexp", slot = "counts")[hvg_genes, ])
########## 2. 发挥服务器多核算力，执行基质细胞 DSE 虚拟基因敲除 (Virtual Knockout)
print("🚀 正在运行 scTenifoldKnk 核心对抗网络算法，执行基质细胞 DSE 虚拟流形对齐敲除...")
knk_results <- scTenifoldKnk(countMatrix = counts_matrix,
                             gKO = "DSE",
                             qc = TRUE)
##### 提取受 DSE 调控级联影响最深重的下游基质驱动基因 (DRGs)
drg_df <- knk_results$diffRegulation
drg_df <- drg_df %>% mutate(p.adj = ifelse(p.adj == 0, 1e-300, p.adj))
drg_df$Significance <- ifelse(drg_df$p.adj < 0.05, "Significant", "Not Significant")
# 掐尖挑出网络扰动强度最剧烈的前 20 个基质靶向基因打上标签
top_genes <- drg_df %>% arrange(desc(abs(Z))) %>% slice(1:20) %>% pull(gene)
drg_df$Label <- ifelse(drg_df$gene %in% top_genes, drg_df$gene, "")
##### 绘制基质细胞中 DSE 虚拟敲除网络火山图 (Fig 6A)
# 🎯 格式安全盾：提前加载 ggrepel，绝对不写在 ggplot 的加号链条内部！
library(ggrepel)
p_volcano <- ggplot(drg_df, aes(x = Z, y = -log10(p.adj), color = Significance)) +
  geom_point(alpha = 0.7, size = 2.5) +
  scale_color_manual(values = c("Not Significant" = "grey80", "Significant" = "#E74C3C")) +
  geom_text_repel(aes(label = Label), size = 4, max.overlaps = 50, box.padding = 0.4, fontface = "bold", segment.size = 0.5) +
  theme_bw() +
  labs(
    title = paste("Virtual KO of DSE in Stromal Cells (SC)"),
    subtitle = "Perturbation of scGRN via scTenifoldKnk",
    x = "Z-score (Network Perturbation Strength)",
    y = "-Log10(FDR)"
  ) +
  theme(
    plot.title = element_text(face = "bold", size = 15, hjust = 0.5),
    plot.subtitle = element_text(face = "italic", size = 11, hjust = 0.5),
    axis.title = element_text(face = "bold", size = 13),
    axis.text = element_text(face = "bold", size = 11),
    legend.position = "top",
    panel.grid = element_blank(),
    panel.border = element_rect(color = "black", linewidth = 1.5)
  )
ggsave("Fig6A_Virtual_KO_DSE_Stromal_Volcano.pdf", p_volcano, width = 8, height = 7, dpi = 300)
print("✅ Fig 6A 基质虚拟敲除网络火山图已完美落盘！")
# =====================================================================
# 🎯【安全标化重构】Fig 6B：解构 DSE 虚拟敲除所瓦解的核心基质功能通路
# =====================================================================
cat("\n🔄 STAGE 1/2: 正在对 DSE 调控的核心基因集进行 GO (BP) 通路解构...\n")

# 筛选网络扰动极其显著的上游核心响应基因进行功能富集
perturbed_genes <- drg_df %>%
  filter(p.adj < 0.05) %>%
  arrange(desc(abs(Z))) %>%
  head(300) %>% # 掐尖网络流形发生剧烈位移的前 300 个核心靶标
  pull(gene)

# 启动 clusterProfiler 引擎
library(clusterProfiler)
library(org.Hs.eg.db)

go_knk_enrich <- enrichGO(
  gene          = perturbed_genes,
  OrgDb         = org.Hs.eg.db,
  keyType       = "SYMBOL",
  ont           = "BP",
  pAdjustMethod = "BH",
  pvalueCutoff  = 0.05
)

# 转化为标准数据框并清洗
go_knk_res <- as.data.frame(go_knk_enrich)

# 🧠 过滤并海选出与基质微环境、炎症重塑、细胞死亡相关的 MS 病理高相关条目
pathway_keywords <- c("death", "apoptosis", "inflammatory", "cytokine", "matrix",
                      "adhesion", "stress", "immune", "wound", "tissue", "response to")
ms_go_res <- go_knk_res %>%
  filter(stringr::str_detect(Description, paste(pathway_keywords, collapse = "|")))

# 降级防御垫：若过滤条目不足，自动切回前 12 个最显著条目
top_go_knk <- if(nrow(ms_go_res) >= 6) head(ms_go_res, 12) else head(go_knk_res, 12)

# 转换 GeneRatio 用于高级可视化排版
top_go_knk <- top_go_knk %>%
  mutate(GeneRatio_num = sapply(GeneRatio, function(x) {
    num <- as.numeric(strsplit(x, "/")[[1]])
    return(num[1] / num[2])
  })) %>%
  arrange(desc(p.adjust)) %>%
  mutate(Description = factor(Description, levels = Description))

# 绘制顶级学术品相的渐变条形图
p_bar <- ggplot(top_go_knk, aes(x = GeneRatio_num, y = Description, fill = -log10(p.adjust))) +
  geom_bar(stat = "identity", color = "black", width = 0.65, linewidth = 0.5) +
  scale_fill_gradientn(colors = c("#FADBD8", "#E74C3C", "#943126"), name = "-log10(p.adjust)") +
  labs(
    title = "Biological Processes Disrupted by DSE KO",
    subtitle = "Top Targeted Matrix & Neuroinflammatory Pathways",
    x = "Gene Ratio (Enrichment Fraction)",
    y = NULL
  ) +
  theme_classic() +
  theme(
    plot.title = element_text(face = "bold", size = 13, hjust = 0),
    plot.subtitle = element_text(face = "italic", size = 10, color = "gray30"),
    axis.text.y = element_text(face = "bold", size = 10, color = "black"),
    axis.text.x = element_text(color = "black", size = 10),
    panel.border = element_rect(color = "black", linewidth = 1.2, fill = NA),
    legend.position = "right"
  )

ggsave("Fig6B_DSE_KO_Disrupted_Pathways.pdf", plot = p_bar, width = 8, height = 5.2, dpi = 300)
cat("✅ Fig 6B 扰动通路条形图已完美落盘！\n")


# =====================================================================
# 🎯【顺序修正版】Fig 6C：具象化绘制 DSE 核心调控的【网络拓扑扰动球】
# =====================================================================
cat("\n🔄 STAGE 2/2: 正在利用 igraph 渲染 DSE 核心拓扑重塑网络图...\n")

if(!require(igraph)) install.packages("igraph")
library(igraph)

# 🎯【核心修复：函数定义强行置顶】先声明弹平缩放工具，彻底消除 could not find 报错
rescale <- function(x, to = c(0, 1)) {
  if(max(x) == min(x)) return(rep(mean(to), length(x))) # 防御单值常数
  (x - min(x)) / (max(x) - min(x)) * (to[2] - to[1]) + to[1]
}

# 掐尖挑出受到 DSE 震荡波波及最惨烈的 Top 20 下游靶基因
network_data <- drg_df %>%
  arrange(desc(abs(Z))) %>%
  filter(gene != "DSE") %>%
  head(20)

# 构建边线矩阵：DSE 居于中央放射状连接所有受害者基因
edges <- data.frame(
  from = "DSE",
  to = network_data$gene,
  weight = abs(network_data$Z) # 扰动强度越强，边线越粗
)

# 构建节点属性矩阵 (此时 rescale 函数已在上方注册，可完美畅通调用)
nodes <- data.frame(
  name = c("DSE", network_data$gene),
  size = c(35, rescale(-log10(network_data$p.adj), to = c(12, 28))),
  z_score = c(max(abs(network_data$Z)) + 2, abs(network_data$Z))
)

# 颜色映射平滑渐变
color_palette <- colorRampPalette(c("#FADBD8", "#E74C3C", "#7B241C"))(100)
node_colors <- c("#2C3E50", color_palette[as.numeric(cut(nodes$z_score[-1], 100))])

# 实例化图对象
g <- graph_from_data_frame(d = edges, vertices = nodes, directed = FALSE)

# 导出高质量无损 PDF 拓扑网络结构
pdf("Fig6C_DSE_Knockout_EgoNetwork.pdf", width = 6.8, height = 6.8)
par(mar = c(1, 1, 2, 1))

set.seed(42)
plot(
  g,
  layout = layout_with_fr(g),
  vertex.size = nodes$size,
  vertex.color = node_colors,
  vertex.frame.color = "black",
  vertex.frame.width = 0.8,
  vertex.label = nodes$name,
  vertex.label.color = "black",
  vertex.label.font = c(2, rep(2, 20)),
  vertex.label.cex = c(1.1, rep(0.85, 20)),
  vertex.label.dist = c(0, rep(0.4, 20)),
  edge.width = rescale(edges$weight, to = c(1.5, 5)),
  edge.color = "gray75",
  main = "DSE Targeted Regulatory Ego-Network"
)

legend("bottomleft", legend = c("Target Gene (DSE)", "High Perturbation", "Low Perturbation"),
       fill = c("#2C3E50", "#7B241C", "#FADBD8"), bty = "n", cex = 0.8)
dev.off()

cat("✅ Fig 6C 拓扑网络结构图已成功落盘！\n")
print("🏁 终极通盘！全套虚拟敲除火山图+富集图+拓扑球资产全部宣告大圆满！")
#######保存数据
save.image("scTenifoldKnk.RData")
