import os
# 切换至结果6工作目录
os.chdir(file.path(Sys.getenv("LEGACY_STUDY_ROOT", unset = "legacy_inputs"), "Result6/"))
import os
import pandas as pd
import drugreflector as dr

######加载并整理 R 语言导出的治疗签名
df = pd.read_csv("DSE_therapeutic_signature.csv", index_col=0)

######DrugReflector 需要一个 Pandas Series 作为输入
#使用 avg_log2FC 作为基因的特征分 (V-score 替代物)
signature = df['avg_log2FC']
#给这个转换命名，输出结果时会用到这个名字
signature.name = 'Stromal:DSE_Low->DSE_High'

######初始化 DrugReflector 深度学习模型
# 🎯【绝对锁定】指向您 Linux 服务器上真实的深度学习权重核心路径
model_paths = [
    os.path.join(os.environ.get('LEGACY_STUDY_ROOT', 'legacy_inputs'), 'Result6/checkpoints/model_fold_0.pt'),
    os.path.join(os.environ.get('LEGACY_STUDY_ROOT', 'legacy_inputs'), 'Result6/checkpoints/model_fold_1.pt'),
    os.path.join(os.environ.get('LEGACY_STUDY_ROOT', 'legacy_inputs'), 'Result6/checkpoints/model_fold_2.pt')
]

######实例化模型 (它会自动检测你的服务器是否有 GPU/CUDA 并加速)
model = dr.DrugReflector(checkpoint_paths=model_paths)

######执行极速虚拟药物筛选
#传入我们的特征签名，让模型去寻找能诱导这种转变的药物
predictions = model.predict(signature, n_top=50)

######查看并保存绝杀结果！
#提取概率最高 (最有可能促使 DSE 表达上调) 的前 10 种化合物
prob_col = ('prob', signature.name)
print("\n🌟 极具潜力的前 10 种靶向/逆转药物 (按神经网络概率排序):")
print(predictions[prob_col].nlargest(10))
# 将预测结果完整保存为 CSV，用作文章附表
predictions.to_csv("DrugReflector_DSE_Candidates_Raw.csv")
######读取LINCS官方文件对药物进行匹配
# 1. 准备你的前 50 名预测结果，并强制提取 13 位核心库编号
prob_col = ('prob', signature.name)
top_predictions = predictions[prob_col].nlargest(50).reset_index()
top_predictions.columns = ['pert_id', 'Probability']
top_predictions['pert_id'] = top_predictions['pert_id'].astype(str).str[:13]

# 2. 读取你最新下载的 compoundinfo 字典
print("正在加载最新版 LINCS 药物字典并执行 13 位核心化清洗...")
drug_info = pd.read_csv(os.path.join(os.environ.get('LEGACY_STUDY_ROOT', 'legacy_inputs'), 'Result6/compoundinfo_beta.txt'), sep='\t')

# 🎯【终极修复盾】必须把官方大字典里的 pert_id 也强行截断为 13 位！熔断因批号后缀不同导致的错位！
drug_info['pert_id'] = drug_info['pert_id'].astype(str).str[:13]

# 提取关键列
cols_to_keep = ['pert_id', 'cmap_name', 'target', 'moa', 'canonical_smiles']
cols_to_keep = [c for c in cols_to_keep if c in drug_info.columns]

# 🎯 此时基于已经清洗干净的 13 位核心编号进行去重，确保合并时绝对一对一
drug_subset = drug_info[cols_to_keep].drop_duplicates(subset=['pert_id'])

# 3. 终极精准合并！
final_annotated_results = pd.merge(top_predictions, drug_subset, on='pert_id', how='left')

print("\n🌟 全新字典精细解码大获全胜！排名前 10 的潜力药物真实身份如下：")
print(final_annotated_results.head(10))

# 导出最终的绝杀附表
final_annotated_results.to_csv("Fig6E_Final_Decoded_DSE_Targeted_Drugs.csv", index=False)
print("\n✅ 带有确切靶点与机制的终极表格已保存为 Fig6E_Final_Decoded_DSE_Targeted_Drugs.csv")