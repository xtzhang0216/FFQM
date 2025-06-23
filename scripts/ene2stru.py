#%%
import pandas as pd
import sys
import os
from rdkit import Chem
model="charmm"
pair="ETOH_MIME"
frag1 = pair.split("_")[0]
frag2 = pair.split("_")[1]
file = f"/pubhome/xtzhang/interaction/FF/{model}/data/{pair}.csv"
writer = Chem.SDWriter(f"/pubhome/xtzhang/interaction/FF/test/{pair}_test_df.{model}.sdf")
df = pd.read_csv(file)
df = df.drop_duplicates(subset=["Index", "QM_Energy", "FF_Energy"])
df = df[~df['Index'].astype(str).str.startswith('rot-')]
# take the first 100 for test
test_df = df.head(100)
print(test_df)

# %%
import json
qmfile = "/pubhome/lzeng/data/pair25/NequipData/TOTAL/total_inteng_comb.json"
with open(qmfile, 'r') as f:
    qmdata = json.load(f)

for key in qmdata.keys():
    if frag1 in key and frag2 in key:
        name = key[:-4]
        qm_df = pd.DataFrame.from_dict(qmdata[key], orient='index', columns=['QM_Energy'])
        qm_df["Index"] = qm_df.index
        qm_df["QM_Energy"] = qm_df["QM_Energy"].round(2)

        match = pd.DataFrame()
        for index, row in test_df.iterrows():
            idx = int(row["Index"])
            qmene = row["QM_Energy"]
            ffene = row["FF_Energy"]
            qm_match = qm_df[((qm_df["QM_Energy"] == qmene)) & (qm_df["Index"] == str(idx))]
            if not qm_match.empty:
                qm_match = qm_match.copy()
                qm_match["FF_Energy"] = ffene
                match = match.append(qm_match)
        if match.empty:
            print(f"Warning: No match found for {name} in {key}")
            continue

        # save as sdf
        with open ("/pubhome/xtzhang/interaction/data/pdbpairs/files.txt", "r") as f:
            lines = f.readlines()
        sdfpath = next((line for line in lines if name in line), None)
        sdfpath = sdfpath.strip()
        sdf = Chem.SDMolSupplier(sdfpath, removeHs=False)
        # write into test_df.sdf

        for idx, mol in enumerate(sdf):
            match["Index"] = match["Index"].astype(int)
            if idx not in match["Index"].values:
                continue
            match_row = match[match["Index"] == (idx)]
            if match_row.empty:
                print(f"Warning: No match found for Index {idx}")
                continue

            expected_qmene = match_row["QM_Energy"].values[0]
            expected_ffene = match_row["FF_Energy"].values[0]

            # 在test_df中查找Index、QM_Energy和FF_Energy都匹配的行
            matching_rows = test_df[
                (test_df["Index"] == idx) &
                (abs(test_df["QM_Energy"] - expected_qmene) < 0.01) &
                (abs(test_df["FF_Energy"] - expected_ffene) < 0.01)
            ]

            if matching_rows.empty:
                print(f"Warning: No matching row found for Index {idx} with QM_Energy={expected_qmene} and FF_Energy={expected_ffene}")
                continue

            test_df_ffene = matching_rows["FF_Energy"].values[0]
            test_df_qmene = matching_rows["QM_Energy"].values[0]

            ffene = test_df_ffene
            qmene = test_df_qmene
            mol.SetProp("FF_Energy", str(ffene))
            mol.SetProp("QM_Energy", str(qmene))
            mol.SetProp("source", name)
            mol.SetProp("Index", str(idx))

            writer.write(mol)
    # %%
