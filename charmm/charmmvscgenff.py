#%%
import pandas as pd

name = "ETAM_ETSH"
charmm = "/pubhome/xtzhang/interaction/FF/charmm/data"
cgenff = "/pubhome/xtzhang/interaction/FF/cgenff/data"
charmm_file = f"{charmm}/{name}.csv"
cgenff_file = f"{cgenff}/{name}.csv"
charmm_df = pd.read_csv(charmm_file)
cgenff_df = pd.read_csv(cgenff_file)
cgenff_df = cgenff_df.round(2)

# %%

common_df = pd.merge(charmm_df, cgenff_df, on='QM_Energy', suffixes=('_charmm', '_cgenff'))
# filter lines that Index_charmm != Index_cgenff
common_df = common_df[common_df['Index_charmm'] != common_df['Index_cgenff']]
# filter enegy > 50
common_df = common_df[common_df['FF_Energy_charmm'] < 50]
#find lines that difference > 5
diff_df = common_df[abs(common_df['FF_Energy_charmm'] - common_df['FF_Energy_cgenff']) > 5]
# %%
# plot the FF_Energy_charmm vs FF_Energy_cgenff
import matplotlib.pyplot as plt
plt.scatter(common_df['FF_Energy_charmm'], common_df['FF_Energy_cgenff'])
plt.xlabel('FF_Energy_charmm')
plt.ylabel('FF_Energy_cgenff')
# plot y=x and corr
plt.plot([-20, 50], [-20, 50], color='red')
plt.text(0, 50, f"corr: {common_df['FF_Energy_charmm'].corr(common_df['FF_Energy_cgenff'])}")
plt.show()

# %%
