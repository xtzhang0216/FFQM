#%%
import json
from pathlib import Path

# Read the energy data from the JSON file
energy_json_path = Path('/pubhome/lzeng/data/pair25/NequipData/TOTAL/total_inteng_comb.json')
with open(energy_json_path, 'r') as energy_file:
    qm_energy_data = json.load(energy_file)

# Now you can use qm_energy_data dictionary in your code
acet = qm_energy_data["ACEM_ACET_00.xyz"]
# find keys in acet that startswith rot
rot_keys = [key for key in acet.keys() if key.startswith("rot")]
# write acet 
with open("/pubhome/xtzhang/interaction/FF/scripts/acet.json", "w") as f:
    json.dump(acet, f)
# %%
frag_1 = "ACET"
frag_2 = "PRPA"
index = 5191
qm_ene = -2.77
energy_json_path = Path('/pubhome/lzeng/data/pair25/NequipData/TOTAL/total_inteng_comb.json')
with open(energy_json_path, 'r') as energy_file:
    qm_energy_data = json.load(energy_file)
with open ("/pubhome/xtzhang/interaction/data/pdbpairs/files.txt", "r") as f:
    lines = f.readlines()
for key in qm_energy_data.keys():
    # local = qm_energy_data[key]
    # print(key,len(local))
    if frag_1 in key and frag_2 in key:
        if str(index) in qm_energy_data[key].keys() and abs(qm_energy_data[key][str(index)] - qm_ene) < 0.1:
            print(key)
        


# %%
frag="ETOH"
for key in qm_energy_data.keys():
    if frag in key:
        print(key)
# %%