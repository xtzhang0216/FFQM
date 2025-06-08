#%%
import pandas as pd
import matplotlib.pyplot as plt
# 写脚本plot/pubhome/xtzhang/interaction/FF/amoeba/data/ACET_PRPA.csv中Tinker_Supermolecular_EnergyvsTinker_Intermolecular_Energy
ACET_PRPA = pd.read_csv("/pubhome/xtzhang/interaction/FF/amoeba/data/ACET_PRPA.csv")
plt.scatter(ACET_PRPA["Tinker_Supermolecular_Energy"], ACET_PRPA["Tinker_Intermolecular_Energy"])
plt.show()
# %%
# 计算差异
ACET_PRPA['diff'] = abs(ACET_PRPA["Tinker_Supermolecular_Energy"] - ACET_PRPA["Tinker_Intermolecular_Energy"])
# 排除能量>50的数据
# ACET_PRPA = ACET_PRPA[ACET_PRPA["Tinker_Supermolecular_Energy"] < 50]
# ACET_PRPA = ACET_PRPA[ACET_PRPA["Tinker_Intermolecular_Energy"] < 50]
# 分离差异大于1kcal的数据
outliers = ACET_PRPA[ACET_PRPA['diff'] > 5]
ACET_PRPA_filtered = ACET_PRPA[ACET_PRPA['diff'] <= 5]

# 重新绘图
plt.figure(figsize=(10,6))
plt.scatter(ACET_PRPA_filtered["Tinker_Supermolecular_Energy"], 
           ACET_PRPA_filtered["Tinker_Intermolecular_Energy"],
           label='Normal points')
# plt.scatter(outliers["Tinker_Supermolecular_Energy"],
#            outliers["Tinker_Intermolecular_Energy"],
#            color='red', label='Outliers (diff > 1 kcal)')
plt.xlabel('Tinker Supermolecular Energy')
plt.ylabel('Tinker Intermolecular Energy')
plt.plot([-10, 10], [-10, 10], 'k--')  # Add y=x line
plt.legend()
plt.show()

# %%
outliers.to_csv("/pubhome/xtzhang/interaction/FF/outliar/wrong/ACET_PRPA/ACET_PRPA_outliers.csv", index=False)
# %%