############## 多发性硬化症单细胞全量细胞图谱丰度统计级级联排版
setwd(file.path(Sys.getenv("LEGACY_STUDY_ROOT", unset = "legacy_inputs"), "Result1"))

############## 加载核心可视化包
library(ggplot2)
library(dplyr)

############## 1. 基于本地 table 实际产出构建高精统计矩阵
# 严格对齐您的实际细胞数，BC 因细胞数为 0 已被物理过滤
cell_counts_summary <- data.frame(
  cell_type = c("OL", "AS", "MG", "NEU", "OPC", "EC", "TC", "SC"),
  Num_Cells = c(64834, 13987, 9239, 8171, 3945, 2230, 835, 277)
)

####### 锁定官方发表级细胞类型颜色表
cell_type_colors <- c("AS" = "#1d77b5",
                      "BC" = "#ff8112",
                      "EC" = "#2a9e65",
                      "MG" = "#d62926",
                      "NEU" = "#aa41fc",
                      "OL" = "#8d5649",
                      "OPC" = "#e278c5",
                      "SC" = "#7f7f7f",
                      "TC" = "#bbbc24")

####### 按照细胞丰度从大到小进行平滑阶梯重排，并锁定因子水平
cell_counts_summary <- cell_counts_summary %>%
  arrange(desc(Num_Cells)) %>%
  mutate(cell_type = factor(cell_type, levels = cell_type))

####### 2. 绘制全量细胞数统计图 (主图)
# 🎯 针对最大值 OL (约 65k) 进行数理平齐：以 10^3 (千) 为单位，Y轴上限锁定 70
p_cells_main <- ggplot(cell_counts_summary, aes(x = cell_type, y = Num_Cells/1000, fill = cell_type)) +
  geom_bar(stat = "identity", color = "black", width = 0.65, linewidth = 0.4) + 
  scale_fill_manual(values = cell_type_colors) +                  
  scale_y_continuous(expand = c(0,0), breaks = seq(0, 70, 10), limits = c(0, 70)) +  
  labs(y = "Number of Cells (*10^3)", x = NULL, title = "All Cell Types") +
  theme_classic() +
  theme(
    axis.text.x = element_text(size = 12, face = "bold", color = "black", angle = 45, hjust = 1),        
    axis.text.y = element_text(size = 12, color = "black"),          
    panel.background = element_blank(),
    panel.border = element_rect(color = "black", linewidth = 1, fill = NA),
    axis.ticks.length = unit(1, "mm"),        
    axis.ticks = element_line(color = "black", linewidth = 1),
    plot.title = element_text(face = "bold", size = 14, hjust = 0.5),
    legend.position = "none"
  )

# 预览并落盘主图
print(p_cells_main)
ggsave("MS_Cell_Counts_Main.pdf", p_cells_main, width = 6.5, height = 5, dpi = 300)


####### 3. 绘制去除高丰度(OL, AS)后的低丰度亚群局部放大图 (Fig1 子图)
# 🎯 筛选标准：锁定细胞数小于 10,000 的稀有/核心响应亚群 (包含 MG, NEU, OPC, EC, TC, SC)
cell_counts_sub <- cell_counts_summary %>% filter(Num_Cells < 10000)

# 🎯 针对局部最大值 MG (9239个) 进行数理标化：以 10^2 (百) 为单位，Y轴上限锁定 100 完美充满画布
p_cells_sub <- ggplot(cell_counts_sub, aes(x = cell_type, y = Num_Cells/100, fill = cell_type)) +
  geom_bar(stat = "identity", color = "black", width = 0.65, linewidth = 0.4) + 
  scale_fill_manual(values = cell_type_colors) +                  
  scale_y_continuous(expand = c(0,0), breaks = seq(0, 100, 20), limits = c(0, 100)) +  
  labs(y = "Number of Cells (*10^2)", x = NULL, title = "Low-abundance Cell Types") +
  theme_classic() +
  theme(
    axis.text.x = element_text(size = 12, face = "bold", color = "black", angle = 45, hjust = 1),        
    axis.text.y = element_text(size = 12, color = "black"),          
    panel.background = element_blank(),
    panel.border = element_rect(color = "black", linewidth = 1, fill = NA),
    axis.ticks.length = unit(1, "mm"),        
    axis.ticks = element_line(color = "black", linewidth = 1),
    plot.title = element_text(face = "bold", size = 12, hjust = 0.5),
    legend.position = "none"
  )

# 预览并落盘子图
print(p_cells_sub)
ggsave("MS_Cell_Counts_Sub.pdf", p_cells_sub, width = 5.5, height = 5, dpi = 300)

cat("🏁 成果一：全量大盘与稀有亚群局部放大图已全部成功无损落盘！\n")