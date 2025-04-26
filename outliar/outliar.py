#%%
import pandas as pd
file = "/pubhome/xtzhang/interaction/FF/charmm/data/MIMM_N1PA.csv"
df = pd.read_csv(file)
# 找到df[ff]-df[qm] > -10的
df["diff"] = df["FF_Energy"] - df["QM_Energy"]
outliar = df[df["diff"]  > 5]
print(outliar)
with open("/pubhome/xtzhang/interaction/FF/outliar/MIMM_N1PA_outliar.charmm.csv", "w") as f:
    outliar.to_csv(f, index=False)
# %%
import json
# 读取所有MIMM_N1PA的内容
qmfile = "/pubhome/xtzhang/interaction/data/pdbpairs/total_inteng_comb.json"
with open(qmfile, 'r') as f:
    qmdata = json.load(f)
key = "MIMM_N1PA.xyz"
name = key[:-4]
qm_df = pd.DataFrame.from_dict(qmdata[key], orient='index', columns=['QM_Energy'])
# 把index新建一列为inedx
qm_df["Index"] = qm_df.index
#保留后两位
qm_df["QM_Energy"] = qm_df["QM_Energy"].round(2)

# %%
# 遍历outliar，对没个item在qm——df找到对应
match = pd.DataFrame()
for index, row in outliar.iterrows():
    idx = int(row["Index"])
    qmene = row["QM_Energy"]
    match = match.append(qm_df[((qm_df["QM_Energy"] == qmene)) & (qm_df["Index"] == str(idx))])
# %%
from rdkit import Chem
with open ("/pubhome/xtzhang/interaction/data/pdbpairs/files.txt", "r") as f:
    lines = f.readlines()
sdfpath = next((line for line in lines if name in line), None)
sdfpath = sdfpath.strip()
sdf = Chem.SDMolSupplier(sdfpath, removeHs=False)
# write into outliar.sdf
writer = Chem.SDWriter("/pubhome/xtzhang/interaction/FF/outliar/MIMM_N1PA_outliar.sdf")

for idx, mol in enumerate(sdf):
    if str(idx) not in match["Index"].values:
        continue
    ffene = outliar[outliar["Index"] == idx]["FF_Energy"].values[0]
    qmene = match[match["Index"] == str(idx)]["QM_Energy"].values[0]
    mol.SetProp("FF_Energy", str(ffene))    
    mol.SetProp("QM_Energy", str(qmene))

    writer.write(mol)
# %%
