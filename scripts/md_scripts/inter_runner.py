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
    mapped_line = next((line for line in lines if key[:-4] in line), None)
    if mapped_line is None:
        print(f"Warning: {key} not found in files.txt")
        continue
    mapped_line = mapped_line.strip()
    cmd = f"python /pubhome/xtzhang/interaction/scripts/md_runner.py --entry {key} --sdf {mapped_line}"
    sbatch = textwrap.dedent(f"""\
            #!/bin/bash
            #SBATCH --job-name={key[:-4]}
            #SBATCH --output=/pubhome/xtzhang/interaction/FF/log/{key[:-4]}.out
            #SBATCH --error=/pubhome/xtzhang/interaction/FF/log/{key[:-4]}.err
            #SBATCH --nodes 1
            #SBATCH --ntasks 1
            #SBATCH --partition mazda,benz
            #SBATCH --qos short_many
            ##SBATCH --begin 2024-08-05T03:00:00
            #SBATCH --time 200:00
            source ~/.bashrc
            conda activate htmd
            """)
    os.system(f"echo '{sbatch}' | sbatch")
    