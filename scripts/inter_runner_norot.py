from pathlib import Path
import json
import os
import textwrap
energy_json_path = Path('/pubhome/lzeng/data/pair25/NequipData/TOTAL/total_inteng_comb.json')
with open(energy_json_path, 'r') as energy_file:
    qm_energy_data = json.load(energy_file)
with open ("/pubhome/xtzhang/interaction/data/pdbpairs/files.txt", "r") as f:
    lines = f.readlines()
for key in qm_energy_data.keys():
    # if "MIMM" not in key:
    #     continue
    # if not any(x in key for x in ["ACEM", "MIMM", "MIME", "MIMD"]):
    #     continue
    # if any(pair in key for pair in ["ETAM_MPHE", "MGDM_MPHE", "ACET_MPHE", "MIMM_MPHE", "ACET_ETOH"]):
    #     continue
    qm_index = list(qm_energy_data[key].keys())
    # Check if all keys can be converted to integers
    non_int_keys = [k for k in qm_index if not k.isdigit()]
    if non_int_keys:
        # with open('/pubhome/xtzhang/interaction/FF/amoeba/rot.log', 'a') as log_file:
        #     log_file.write(f"Warning: Non-integer keys found in {key}: {', '.join(non_int_keys)}\n")
        pass
    # Sort the remaining keys (which are integers) in ascending order
    int_keys = [k for k in qm_index if k.isdigit()]
    sorted_int_keys = sorted(int_keys, key=int)
    batch_size = 500
    for start in range(0, len(sorted_int_keys), batch_size):
        start_idx = sorted_int_keys[start]
        end_idx = sorted_int_keys[min(start + batch_size, len(sorted_int_keys) - 1)]
        # check if can be int, else print out and continue
        try:
            start_idx = int(start_idx)
            end_idx = int(end_idx)
        except ValueError:
            # with open('/pubhome/xtzhang/interaction/FF/amoeba/rot.log', 'a') as log_file:
            #     log_file.write(f"Warning: {key} not found in files.txt\n")
            continue
        cmd = f"python /pubhome/xtzhang/interaction/FF/charmm/openmm_charmm.py --xyz {key[:-4]} --start {start_idx} --end {end_idx} --check_nh True"
        sbatch = textwrap.dedent(f"""\
                #!/bin/bash
                #SBATCH --job-name={key[:-4]}
                #SBATCH --output=/pubhome/xtzhang/interaction/FF/log/{key[:-4]}_%j.out
                #SBATCH --error=/pubhome/xtzhang/interaction/FF/log/{key[:-4]}_%j.err
                #SBATCH --nodes 1
                #SBATCH --ntasks 1
                #SBATCH --partition mazda,benz
                #SBATCH --qos short_many
                ##SBATCH --begin 2024-08-05T03:00:00
                #SBATCH --time 200:00
                source ~/.bashrc
                conda activate htmd
                {cmd}
                """)
        # print(cmd)
        os.system(f"echo '{sbatch}' | sbatch")
    