######### 🎯【成果六·Input篇】DSE 治疗特征谱之 GO & KEGG 双通道功能富集全图链
# =====================================================================
# 1. 设置工作环境与加载核心包
# =====================================================================
setwd(file.path(Sys.getenv("LEGACY_STUDY_ROOT", unset = "legacy_inputs"), "Result6"))

library(clusterProfiler)
library(org.Hs.eg.db)
library(ggplot2)
library(dplyr)
library(stringr)

# =====================================================================
# 2. 读取 DSE 转录特征谱并筛选核心响应基因
# =====================================================================
print("🔓 正在解冻 DSE 基质响应转录特征谱...")
sig_data <- read.csv("DSE_therapeutic_signature.csv", row.names = 1)
sig_data$gene <- rownames(sig_data)

# 筛选核心基因 (|log2FC| > 0.25 且调整后 p 值显著)
sig_genes <- sig_data %>%
  filter(p_val_adj < 0.05 & abs(avg_log2FC) > 0.25) %>%
  pull(gene)

# =====================================================================
# 3. 【通道一】运行 GO (Biological Process) 富集分析与可视化
# =====================================================================
print("🔄 正在轰鸣 GO (Biological Process) 富集分析引擎...")
go_enrich <- enrichGO(
  gene          = sig_genes,
  OrgDb         = org.Hs.eg.db,
  keyType       = "SYMBOL",
  ont           = "BP",
  pAdjustMethod = "BH",
  pvalueCutoff  = 0.05,
  qvalueCutoff  = 0.05
)

go_res <- as.data.frame(go_enrich)

# 提取显著性最高的前 12 个 GO 条目
top_go <- go_res %>%
  arrange(p.adjust) %>%
  head(12) %>%
  mutate(GeneRatio_num = sapply(GeneRatio, function(x) {
    num <- as.numeric(strsplit(x, "/")[[1]])
    return(num[1] / num[2])
  })) %>%
  mutate(Description = factor(Description, levels = rev(Description)))

# 绘制 GO 发表级暖色系气泡图
p_go <- ggplot(top_go, aes(x = GeneRatio_num, y = reorder(Description, GeneRatio_num))) +
  geom_point(aes(size = Count, fill = -log10(p.adjust)), shape = 21, color = "black", stroke = 0.6) +
  scale_fill_gradientn(
    colors = c("#fef0d9", "#fdcc8a", "#fc8d59", "#d7301f"), 
    name = "-log10(p.adjust)",
    limits = c(0, max(-log10(top_go$p.adjust)))
  ) +
  scale_size_continuous(name = "Gene Count", range = c(3.8, 8.5)) +
  labs(
    title = "Column 1: Input Signature (GO)",
    subtitle = "DSE-altered Biological Processes",
    x = "Gene Ratio",
    y = NULL
  ) +
  theme_bw() + 
  theme(
    plot.title = element_text(size = 14, face = "bold", color = "black"),
    plot.subtitle = element_text(size = 11, face = "italic", color = "gray30", margin = margin(b = 10)),
    axis.text.y = element_text(size = 10.5, face = "bold", color = "black", hjust = 1), 
    axis.text.x = element_text(size = 10, color = "black"),
    axis.title.x = element_text(size = 12, face = "bold", margin = margin(t = 10)),
    panel.grid.major.y = element_line(color = "gray90", linetype = "dashed"),
    panel.grid.major.x = element_line(color = "gray95"),
    panel.grid.minor = element_blank(),
    legend.title = element_text(size = 10, face = "bold"),
    legend.text = element_text(size = 9),
    legend.position = "right"
  )

ggsave("Fig6_Column1_DSE_GO_Pathways.pdf", plot = p_go, width = 7.8, height = 5.2, dpi = 300)
print("✅ GO 气泡图已无损落盘。")

# =====================================================================
# 4. 【通道二】运行 KEGG 通路富集分析（含 MS 病理靶向纯化）与可视化
# =====================================================================
print("🔄 正在执行 SYMBOL 到 ENTREZID 的矩阵级编码转换...")
gene_ids <- bitr(sig_genes, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)

print("🚀 正在全面轰鸣 KEGG 富集分析引擎...")
kegg_enrich <- enrichKEGG(
  gene          = gene_ids$ENTREZID,
  organism      = "hsa",
  keyType       = "kegg",
  pAdjustMethod = "BH",
  pvalueCutoff  = 0.05,
  qvalueCutoff  = 0.05
)

kegg_res <- as.data.frame(kegg_enrich)

# 智能化海选并过滤出与多发性硬化症 (MS) 免疫炎症高度相关的核心条目
ms_keywords <- c("sclerosis", "Th17", "Th1", "NF-kappa B", "TNF", "Chemokine", 
                 "Cytokine", "JAK-STAT", "Antigen", "Apoptosis", "T cell", 
                 "B cell", "Toll-like", "NOD-like", "Necroptosis", "Autoimmune")

ms_kegg_res <- kegg_res %>%
  filter(stringr::str_detect(Description, paste(ms_keywords, collapse = "|")))

# 容错降级机制
if(nrow(ms_kegg_res) < 5) {
  top_kegg <- kegg_res %>% arrange(p.adjust) %>% head(12)
} else {
  top_kegg = ms_kegg_res %>% arrange(p.adjust) %>% head(12)
}

top_kegg <- top_kegg %>%
  mutate(GeneRatio_num = sapply(GeneRatio, function(x) {
    num <- as.numeric(strsplit(x, "/")[[1]])
    return(num[1] / num[2])
  })) %>%
  mutate(Description = factor(Description, levels = rev(Description)))

# 绘制 KEGG 发表级暖色系气泡图
p_kegg <- ggplot(top_kegg, aes(x = GeneRatio_num, y = reorder(Description, GeneRatio_num))) +
  geom_point(aes(size = Count, fill = -log10(p.adjust)), shape = 21, color = "black", stroke = 0.6) +
  scale_fill_gradientn(
    colors = c("#fef0d9", "#fdcc8a", "#fc8d59", "#d7301f"), 
    name = "-log10(p.adjust)",
    limits = c(0, max(-log10(top_kegg$p.adjust)))
  ) +
  scale_size_continuous(name = "Gene Count", range = c(3.8, 8.5)) +
  labs(
    title = "Column 1: Input Signature (KEGG)",
    subtitle = "DSE-associated Pathogenic KEGG Pathways in MS",
    x = "Gene Ratio",
    y = NULL
  ) +
  theme_bw() + 
  theme(
    plot.title = element_text(size = 14, face = "bold", color = "black"),
    plot.subtitle = element_text(size = 11, face = "italic", color = "gray30", margin = margin(b = 10)),
    axis.text.y = element_text(size = 11, face = "bold", color = "black", hjust = 1), 
    axis.text.x = element_text(size = 10, color = "black"),
    axis.title.x = element_text(size = 12, face = "bold", margin = margin(t = 10)),
    panel.grid.major.y = element_line(color = "gray90", linetype = "dashed"),
    panel.grid.major.x = element_line(color = "gray95"),
    panel.grid.minor = element_blank(),
    legend.title = element_text(size = 10, face = "bold"),
    legend.text = element_text(size = 9),
    legend.position = "right"
  )

ggsave("Fig6_Column1_DSE_KEGG_Pathways.pdf", plot = p_kegg, width = 7.8, height = 5.2, dpi = 300)

print("🏁 普天同庆！GO 与 KEGG 双通道富集分析全部圆满落盘！")
