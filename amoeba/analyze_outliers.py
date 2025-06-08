#%%
import pandas as pd
import numpy as np

# 读取数据
df = pd.read_csv('/pubhome/xtzhang/interaction/FF/amoeba/data/ACEM_ACET.debug.csv')

# 计算FF和QM的差值
df['FF_QM_diff'] = df['FF_Energy'] - df['QM_Energy']

# 提取差值<-10的行
outliers = df[df['FF_QM_diff'] < -10]

# 按差值排序
outliers_sorted = outliers.sort_values('FF_QM_diff')

# 打印结果
print("找到的异常值数量:", len(outliers))
print("\n差值从小到大排序的结果:")
print(outliers_sorted[['Index', 'FF_Energy', 'QM_Energy', 'FF_QM_diff']].to_string())

# 保存结果到文件
output_file = 'ACET_N1PA_outliers.csv'
outliers_sorted.to_csv(output_file, index=False)
print(f"\n结果已保存到 {output_file}")
# %%
