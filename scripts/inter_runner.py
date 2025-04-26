from pathlib import Path
import json
import os
import textwrap
energy_json_path = Path('/pubhome/xtzhang/interaction/data/pdbpairs/total_inteng_comb.json')
with open(energy_json_path, 'r') as energy_file:
    qm_energy_data = json.load(energy_file)
with open ("/pubhome/xtzhang/interaction/data/pdbpairs/files.txt", "r") as f:
    lines = f.readlines()
for key in qm_energy_data.keys():
    # if key[:-4] not in ["ACET_nNMA_76", "ACET_nNMA_37", "ACET_nNMA_29","ACET_ETAM_07","ACET_ETAM_06","ACET_ETAM_03","ACET_ACET_01"]:
    #     continue

    cmd = f"python /pubhome/xtzhang/interaction/FF/opls/openmm_opls.py {key[:-4]}"
    sbatch = textwrap.dedent(f"""\
            #!/bin/bash
            #SBATCH --job-name={key[:-4]}
            #SBATCH --output=/pubhome/xtzhang/interaction/FF/log/{key[:-4]}_%j.out
            #SBATCH --error=/pubhome/xtzhang/interaction/FF/log/{key[:-4]}_%j.err
            #SBATCH --nodes 1
            #SBATCH --ntasks 1
            #SBATCH --partition mazda,benz
            #SBATCH --qos normal
            ##SBATCH --begin 2024-08-05T03:00:00
            #SBATCH --time 720:00
            source ~/.bashrc
            conda activate htmd
            {cmd}
            """)
    os.system(f"echo '{sbatch}' | sbatch")
    