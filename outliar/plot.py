#%%
import pandas as pd
file = "/pubhome/xtzhang/interaction/FF/outliar/outliar_head/ETAM_MSM.csv"
df = pd.read_csv(file)
# filter that diff > 20
df = df[df['diff'] > 5]
# 取前100个
# df = df.head(50)
# plot data in ETAM_MSM,柱状图,vdw.静电一样粗细,qm能量最粗
import matplotlib.pyplot as plt
import numpy as np
x = np.arange(len(df['Index']))
width = 0.25  # 柱子的宽度

plt.bar(x - width, df['QM_Energy'], width, color='blue', label='QM_Energy')
plt.bar(x, df['vdw'], width, color='red', label='vdw')
plt.bar(x + width, df['chg'], width, color='green', label='chg')

plt.xlabel('Molecule')
plt.ylabel('Energy (kcal/mol)')
plt.title('ETAM_MSM head outliar')
plt.xticks([])  # 移除x轴刻度
plt.legend()
plt.show()
# %%
