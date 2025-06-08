#%%
import pandas as pd
import sys
import os
from rdkit import Chem
# model = sys.argv[1]
# pair = sys.argv[2]
model="amoeba"
pair="ACET_PRPA"
frag1 = "ACET"
frag2 = "PRPA"
# key = sys.argv[2]
# name = key[:-4]
# cutoff = float(sys.argv[3])
file = f"/pubhome/xtzhang/interaction/FF/{model}/data/{pair}.csv"
df = pd.read_csv(file)
df = df.drop_duplicates(subset=["Index", "QM_Energy", "Tinker_Supermolecular_Energy"])
df = df[~df['Index'].astype(str).str.startswith('rot-')]
# Filter out rows where Tinker_Supermolecular_Energy is less than 50
df = df[df["Tinker_Supermolecular_Energy"] < 50]
df = df[df["QM_Energy"] < 50]
df = df[df["Tinker_Supermolecular_Energy"] < 0]
# 找到df[ff]-df[qm] > -10的
df["diff"] = df["QM_Energy"] - df["Tinker_Supermolecular_Energy"]
outliar = df[abs(df["diff"]) > 10]
print(outliar)
if not os.path.exists(f"/pubhome/xtzhang/interaction/FF/outliar/wrong/{pair}"):
    os.makedirs(f"/pubhome/xtzhang/interaction/FF/outliar/wrong/{pair}")
# with open(f"/pubhome/xtzhang/interaction/FF/outliar/wrong/{pair}/{pair}_outliar.{model}.csv", "w") as f:
#     outliar.to_csv(f, index=False)
# %%
import json
# 读取所有MIMM_N1PA的内容
qmfile = "/pubhome/lzeng/data/pair25/NequipData/TOTAL/total_inteng_comb.json"
with open(qmfile, 'r') as f:
    qmdata = json.load(f)
writer = Chem.SDWriter(f"/pubhome/xtzhang/interaction/FF/outliar/wrong/{pair}/{pair}_outliar.{model}.sdf")

for key in qmdata.keys():
    if frag1 in key and frag2 in key:
        name = key[:-4]
        # break
        qm_df = pd.DataFrame.from_dict(qmdata[key], orient='index', columns=['QM_Energy'])
        # 把index新建一列为inedx
        qm_df["Index"] = qm_df.index
        #保留后两位
        qm_df["QM_Energy"] = qm_df["QM_Energy"].round(2)

        # 遍历outliar，对每个item在qm——df找到对应
        match = pd.DataFrame()
        for index, row in outliar.iterrows():
            idx = int(row["Index"])
            qmene = row["QM_Energy"]
            ffene = row["Tinker_Supermolecular_Energy"]
            # 要求index、qm_ene都匹配，并记录ff_ene用于后续验证
            qm_match = qm_df[((qm_df["QM_Energy"] == qmene)) & (qm_df["Index"] == str(idx))]
            if not qm_match.empty:
                # 添加Tinker_Supermolecular_Energy信息到匹配结果中
                qm_match = qm_match.copy()
                qm_match["Tinker_Supermolecular_Energy"] = ffene
                match = match.append(qm_match)
        if match.empty:
            print(f"Warning: No match found for {name} in {key}")
            continue
        with open ("/pubhome/xtzhang/interaction/data/pdbpairs/files.txt", "r") as f:
            lines = f.readlines()
        sdfpath = next((line for line in lines if name in line), None)
        sdfpath = sdfpath.strip()
        sdf = Chem.SDMolSupplier(sdfpath, removeHs=False)
        # write into outliar.sdf

        for idx, mol in enumerate(sdf):
            match["Index"] = match["Index"].astype(int)
            if idx not in match["Index"].values:
                continue
            # matching_rows = outliar[outliar["Index"] == str(idx)]
            # if matching_rows.empty:
            # 获取match中的能量值用于匹配
            match_row = match[match["Index"] == (idx)]
            if match_row.empty:
                print(f"Warning: No match found for Index {idx}")
                continue

            expected_qmene = match_row["QM_Energy"].values[0]
            expected_ffene = match_row["Tinker_Supermolecular_Energy"].values[0]

            # 在outliar中查找Index、QM_Energy和Tinker_Supermolecular_Energy都匹配的行
            matching_rows = outliar[
                (outliar["Index"] == idx) &
                (abs(outliar["QM_Energy"] - expected_qmene) < 0.01) &
                (abs(outliar["Tinker_Supermolecular_Energy"] - expected_ffene) < 0.01)
            ]

            if matching_rows.empty:
                print(f"Warning: No matching row found for Index {idx} with QM_Energy={expected_qmene} and Tinker_Supermolecular_Energy={expected_ffene}")
                continue

            # 获取匹配行中的Tinker_Supermolecular_Energy和QM_Energy
            outliar_ffene = matching_rows["Tinker_Supermolecular_Energy"].values[0]
            outliar_qmene = matching_rows["QM_Energy"].values[0]

            # 使用已经验证过的能量值
            ffene = outliar_ffene
            qmene = outliar_qmene
            mol.SetProp("Tinker_Supermolecular_Energy", str(ffene))
            mol.SetProp("QM_Energy", str(qmene))
            mol.SetProp("source", name)
            mol.SetProp("Index", str(idx))

            writer.write(mol)
    # %%
