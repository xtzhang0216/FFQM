#%%
import pandas as pd
file = "/pubhome/xtzhang/interaction/FF/outliar/tail/ACEM_ETAM/ACEM_ETAM.charmm.csv"
df = pd.read_csv(file)
# filter that diff > 20
# df = df[df['diff'] > 5]
print(len(df))
df = df[df['QM_Energy'] < 100]
# 取前100个
df = df.sort_values(by='diff', key=abs, ascending=False)
df = df.head(500)
# plot data in MIMM_N1PA,柱状图,Vdw_Energy.静电一样粗细,qm能量最粗
import matplotlib.pyplot as plt
import numpy as np
x = np.arange(len(df['Index']))
width = 0.25  # 柱子的宽度
plt.figure(figsize=(50, 20))
plt.bar(x - width, df['QM_Energy'], width, color='blue', label='QM_Energy')
plt.bar(x, df['Vdw_Energy'], width, color='red', label='Vdw_Energy')
plt.bar(x + width, df['Chg_Energy'], width, color='green', label='Ele_Energy')

plt.xlabel('Molecule', fontsize=50)
plt.ylabel('Energy (kcal/mol)', fontsize=50)
plt.title(f'{pair} outliar', fontsize=50)
plt.xticks([])  # /移除x轴刻度
plt.legend(fontsize=50)
plt.xticks(fontsize=40)  # 调大x轴刻度字体
plt.yticks(fontsize=40)  # 调大y轴刻度字体

plt.show()
# plt.savefig("/pubhome/xtzhang/interaction/FF/outliar/outliar_head/MIMM_N1PA/MIMM_N1PA_outliar.png")
# %%
