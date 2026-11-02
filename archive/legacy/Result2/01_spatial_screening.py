import os
########## 查看Cell type mapping reveals tissue niches and interactions in subcortical multiple sclerosis lesions提供的空转切片那个最适合用在我们论文里面，然后作图
##### 加载库
import os
import scanpy as sc
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.stats import pearsonr

# 🎯【绝杀乱码】配置全局绘图引擎，锁死 TrueType 字体嵌入模式，彻底扫平 Adobe Illustrator 导入文本乱码与方块问题
plt.rcParams['pdf.fonttype'] = 42
plt.rcParams['ps.fonttype'] = 42
plt.rcParams['font.sans-serif'] = ['Arial', 'Helvetica', 'DejaVu Sans']  # 优先采用 AI 百分之百原生支持的标准无衬线字体
plt.rcParams['font.family'] = 'sans-serif'

##### 环境配置与全局变量
data_dir = os.path.join(os.environ.get('LEGACY_STUDY_ROOT', 'legacy_inputs'), 'Data/ST_Data')
output_dir = os.path.join(os.environ.get('LEGACY_STUDY_ROOT', 'legacy_inputs'), 'Result2/')
os.makedirs(output_dir, exist_ok=True)

##### 14 个基质特异性双硫死亡核心基因
disulfidptosis_genes = [
    "FMO2", "SLC6A20", "PTGIS", "BHMT2", "HSD11B1", "DSE", "OSBPL5",
    "ATP10D", "IRS1", "CALU", "SEC61A1", "RRBP1", "PRKAA1", "MTAP"
]

##### 2 个核心炎症神经重构转录因子
tfs = ["CEBPB", "RELA"]
h5ad_files = [f for f in os.listdir(data_dir) if f.endswith('.h5ad')]
summary_stats = []

print(f"🚀 开始对 {len(h5ad_files)} 张多发性硬化症空转切片进行双硫死亡评分与金牌筛选...")

##### 遍历全切片：计算双硫死亡评分与炎症关联度并筛选最佳切片
for file_name in h5ad_files:
    file_path = os.path.join(data_dir, file_name)
    try:
        adata = sc.read_h5ad(file_path)
        # 强制基因名大写对齐
        adata.var_names = adata.var_names.str.upper()
        # 筛选当前切片中存在的双硫死亡基因
        valid_drgs = [g for g in disulfidptosis_genes if g in adata.var_names]
        if len(valid_drgs) < 5:
            continue
        # 核心算法：利用 Scanpy 官方标准的 score_genes 计算空间双硫死亡评分
        sc.tl.score_genes(adata, gene_list=valid_drgs, score_name='disulfidptosis_score')

        # 🎯【安全标化】使用 obs_vector 提取分子空间丰度，彻底斩断遍历扫描时的维度崩溃隐患
        if "CEBPB" in adata.var_names and "RELA" in adata.var_names:
            cebpb_val = adata.obs_vector("CEBPB")
            rela_val = adata.obs_vector("RELA")
            dis_score = adata.obs['disulfidptosis_score'].values
            r_cebpb, _ = pearsonr(dis_score, cebpb_val)
            r_rela, _ = pearsonr(dis_score, rela_val)
        else:
            r_cebpb, r_rela = 0, 0

        # 获取病灶分型标签参考
        lesion_types = "/".join(list(adata.obs['lesion_type'].dropna().unique())[:2]) if 'lesion_type' in adata.obs.columns else "Unknown"
        summary_stats.append({
            "File_Name": file_name,
            "Total_Spots": adata.n_obs,
            "Mean_Disulfidptosis": adata.obs['disulfidptosis_score'].mean(),
            "R_with_CEBPB": r_cebpb,
            "R_with_RELA": r_rela,
            "Score_Sum": r_cebpb + r_rela, # 用于综合评价炎症共定位强度的指标
            "Lesion_Types": lesion_types
        })
    except Exception as e:
        print(f"❌ 预扫描切片 {file_name} 时发生微小跳过: {e}")

# 输出并展示筛选排名表
df_summary = pd.DataFrame(summary_stats).sort_values(by="Score_Sum", ascending=False)
print("\n🏆 ======= 空间转录组切片病理特征与炎症关联度最终金牌排名 ======= 🏆")
print(df_summary.to_string(index=False))

##### 锁定最佳切片并配置 Scanpy 绘图导出环境
best_file = df_summary.iloc[0]["File_Name"]
print(f"\n🎯 系统已自动为您锁定综合质量最佳的冠军切片: [{best_file}]，正在全自动加载绘图...")

# 直接从原始路径读取金牌切片，确保内存和数据逻辑绝对自愈合
best_adata = sc.read_h5ad(os.path.join(data_dir, best_file))
best_adata.var_names = best_adata.var_names.str.upper()

# 重新在内存金牌切片中精准算分
valid_drgs_best = [g for g in disulfidptosis_genes if g in best_adata.var_names]
sc.tl.score_genes(best_adata, gene_list=valid_drgs_best, score_name='disulfidptosis_score')

# 全局出图配置归档至 Result2 路径下
sc.settings.figdir = output_dir
sc.settings.set_figure_params(dpi=300, facecolor='white', frameon=False)

try:
    library_id = list(best_adata.uns['spatial'].keys())[0]
    available_img_keys = list(best_adata.uns['spatial'][library_id]['images'].keys())
    actual_img_key = available_img_keys[0] if len(available_img_keys) > 0 else None
except Exception as e:
    actual_img_key = None

##### 绘制空间原位生态位图 (Fig 2A/B)
# 🎯 受惠于顶部的 TrueType 配置，这张由 Scanpy 吐出的原位切片 PDF 里的文字在 AI 里也将能够无损二次编辑！
print("正在导出 Fig2_Niche.pdf ...")
sc.pl.spatial(
    best_adata,
    color=['disulfidptosis_score', 'CEBPB', 'RELA'],
    img_key=actual_img_key,
    size=1.5,
    cmap="jet",
    show=False,
    save="Fig2_Niche.pdf"
)

##### 绘制空间原位 Pearson 相关性热图 (Fig 2C)
# 🎯【核心重构】全部升级为最稳健的 obs_vector 圣剑机制，彻底粉碎矩阵维度退化陷阱，拒绝任何 NaN 引起的空白图
corr_df = pd.DataFrame({
    "Disulfidptosis": best_adata.obs['disulfidptosis_score'].values,
    "CEBPB (TF)": best_adata.obs_vector("CEBPB"),
    "RELA (TF)": best_adata.obs_vector("RELA"),
    "DSE (SC Marker)": best_adata.obs_vector("DSE"),
    "HSPG2 (SC Marker)": best_adata.obs_vector("HSPG2")
}).fillna(0) # 兜底填充，保全全样本矩阵完整性

corr_matrix = corr_df.corr(method="pearson")
corr_matrix.to_csv(os.path.join(output_dir, "Fig2C_Spatial_Corr_Data.csv"))

# 🎯 显式创建独立的子图轴，强制锁定渲染上下文
fig, ax = plt.subplots(figsize=(6.5, 5.5))

# 利用 seaborn 绘制精美出版热图
sns.heatmap(
    corr_matrix,
    annot=True,
    fmt=".3f",
    cmap="RdBu_r",
    vmin=-1, vmax=1,
    linewidths=0.5,
    cbar_kws={"label": "Spatial Pearson Correlation (r)"},
    annot_kws={"weight": "bold", "size": 10}, # 矩阵内部数值加粗显示
    ax=ax
)

# 对热图的文本轴实施无损加粗控制
ax.set_title("Spatial Correlation Matrix in Local Niche", fontsize=12, fontweight='bold', pad=12)
plt.xticks(rotation=45, ha="right", weight="bold")
plt.yticks(rotation=0, weight="bold")

plt.tight_layout()

# 严格执行：连续保存无损 PDF 格式，最后统一关闭画布防止数据流失
print("正在导出 Fig2C_Spatial_Correlation_Heatmap.pdf (已完美嵌入AI原生编辑文本)...")
plt.savefig(os.path.join(output_dir, "Fig2C_Spatial_Correlation_Heatmap.pdf"), dpi=300, transparent=True)
plt.close()

print("\n🏁 全自动高精筛选与 AI 友好型矢量图件导出全盘大功告成！")
##### 病灶微解剖区 + 双硫评分 + SC 标志物原位共定位四联图 (Fig 2D)
# 采用更细腻的空间 areas 标签完美勾勒解剖异质性，仅保存无损 PDF 版本
print("正在为您总攻绘制结果2的顶级微解剖空间共定位多联图...")
plot_features = ['areas', 'disulfidptosis_score', 'DSE', 'HSPG2']
sc.pl.spatial(
    best_adata,
    color=plot_features,
    img_key=actual_img_key,
    size=1.5,
    cmap="jet",
    show=False,
    save="Fig2_Areas_vs_Disulfidptosis_SC_Markers.pdf"
)
##### 跨组学对比：健康对照组 vs 疾病组各亚区双硫死亡评分横向大比拼 (Fig 2E)
# 自动化载入扫描出的第一张干净健康对照白质数据（visium_CO41.h5ad）进行横向绝杀
try:
    print("正在智能计算健康对照组的空间双硫死亡基线评分...")
    adata_ctrl = sc.read_h5ad(os.path.join(data_dir, "visium_CO41.h5ad"))
    adata_ctrl.var_names = adata_ctrl.var_names.str.upper()
    valid_drgs_ctrl = [g for g in disulfidptosis_genes if g in adata_ctrl.var_names]
    sc.tl.score_genes(adata_ctrl, gene_list=valid_drgs_ctrl, score_name='disulfidptosis_score')
    # 构建并洗涤融合数据矩阵
    df_ms = pd.DataFrame({
        "Area": best_adata.obs['areas'].astype(str),
        "Score": best_adata.obs['disulfidptosis_score'].values
    })
    df_ctrl = pd.DataFrame({
        "Area": adata_ctrl.obs['areas'].astype(str) + "_Ctrl",
        "Score": adata_ctrl.obs['disulfidptosis_score'].values
    })
    plot_df = pd.concat([df_ms, df_ctrl], axis=0)
    plot_df = plot_df[plot_df['Area'].isin(['LC', 'LR', 'PPWM', 'WM_Ctrl'])]
    # 映衬微环境专业英文出版名称映射
    name_map = {
        "WM_Ctrl": "Healthy WM\n(Control)",
        "PPWM": "Peri-Plaque WM\n(MS)",
        "LR": "Lesion Rim\n(MS)",
        "LC": "Lesion Core\n(MS)"
    }
    plot_df['Display_Area'] = plot_df['Area'].map(name_map)
    order = ["Healthy WM\n(Control)", "Peri-Plaque WM\n(MS)", "Lesion Rim\n(MS)", "Lesion Core\n(MS)"]
    # 绘制高雅的出版箱线小提琴组合图
    plt.figure(figsize=(7, 5))
    my_colors = ["#3B4992FF", "#6388B4FF", "#EE0000FF", "#B2182BFA"]
    sns.violinplot(data=plot_df, x="Display_Area", y="Score", order=order, palette=my_colors, inner=None, linewidth=0.6, alpha=0.4)
    sns.boxplot(data=plot_df, x="Display_Area", y="Score", order=order, width=0.25, palette=my_colors, boxprops=dict(alpha=0.9), showfliers=False, linewidth=0.8)
    plt.title("Disulfidptosis Score across Spatial Microenvironments", fontsize=12, fontweight='bold')
    plt.xlabel("")
    plt.ylabel("Spatial Disulfidptosis Score")
    plt.xticks(fontsize=9)
    plt.tight_layout()
    print("正在保存横向大比拼图 Fig2E_MS_vs_Control_Score_Comparison.pdf ...")
    plt.savefig(os.path.join(output_dir, "Fig2E_MS_vs_Control_Score_Comparison.pdf"), dpi=300)
    plt.close()
    print(f"\n🏆 功德圆满！结果 2 的全套空间重构无损 PDF 矢量景观图谱已完美完工！")
    print(f"📦 产出归档目录：{output_dir}")
except Exception as e:
    print(f"⚠️ 在进行健康组学横向比拼图绘制时发生微小异常: {e}")