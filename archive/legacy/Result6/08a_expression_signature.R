###################### 多发性硬化症单细胞进行 DSE 虚拟药物筛选
############ R语言提取基质细胞中 DSE-High vs DSE-Low 的全局差异表达谱
###### 设置工作环境
setwd(file.path(Sys.getenv("LEGACY_STUDY_ROOT", unset = "legacy_inputs"), "Result6"))

###### 加载R包
library(Seurat)

###### 读取多发性硬化症单细胞全景图谱数据
print("🚀 正在载入中枢神经系统单细胞表达谱大数据集...")
GSE279180_sn_atlas <- readRDS(Sys.getenv("ATLAS_RDS", unset = ""))

### 🎯【核心增加】兼容性对齐身份：洗涤并纯化出基质细胞 (SC) 亚群进行靶向计算
if ("celltype" %in% colnames(GSE279180_sn_atlas@meta.data)) {
  Idents(GSE279180_sn_atlas) <- "celltype"
} else {
  Idents(GSE279180_sn_atlas) <- "cell_type"
}
endo_cells <- subset(GSE279180_sn_atlas, idents = "SC")

# 自动锁定结果 1 的关键 Assay 槽位 (优先 originalexp)
assay_name <- ifelse("originalexp" %in% Assays(endo_cells), "originalexp", "RNA")
DefaultAssay(endo_cells) <- assay_name

### 以 DSE 的中位表达量，将基质细胞分为 High 和 Low 组
med_val <- median(endo_cells@assays[[assay_name]]@data["DSE", ])

# 容错机制：如果中位数为 0（单细胞零膨胀常见），则改用均值切分，防止高低两组细胞数极端失衡
if(med_val == 0) { med_val <- mean(endo_cells@assays[[assay_name]]@data["DSE", ]) }

endo_cells$DSE_status <- ifelse(endo_cells@assays[[assay_name]]@data["DSE", ] > med_val, "High", "Low")

### 设置聚类标识为 DSE_status
Idents(endo_cells) <- "DSE_status"

### 计算 High 组对比 Low 组的差异表达
### 💡 提示：阈值设为 0，因为神经网络需要看全景基因的变化趋势，而不仅是显著基因
print("🚀 正在计算基质细胞中 DSE 暴走表型的全景神经网络特征签名...")
therapeutic_sig <- FindMarkers(endo_cells, ident.1 = "High", ident.2 = "Low",
                               logfc.threshold = 0, min.pct = 0, assay = assay_name)

### 导出为 CSV 文件供 Python 读取
write.csv(therapeutic_sig, "DSE_therapeutic_signature.csv", row.names = TRUE)
cat("✅ 逆向基质特征签名已成功导出为: DSE_therapeutic_signature.csv\n")
############ Python 中运行 DrugReflector 虚拟筛选
######加载库
