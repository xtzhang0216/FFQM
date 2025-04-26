# plop ff vs qm
#%%
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
file = "/pubhome/xtzhang/interaction/FF/scripts/results.csv"
df = pd.read_csv(file)
df = df.dropna()
df = df[df['QM_Energy'] <= 100]
df = df[df['FF_Energy'] <= 100]
plt.scatter(df['FF_Energy'], df['QM_Energy'])
plt.xlabel('FF Energy')
plt.ylabel('QM Energy')
x = np.linspace(min(df['FF_Energy']), max(df['FF_Energy']), 100)
y = x
plt.plot(x, y, color='red')
plt.show()

# %%
