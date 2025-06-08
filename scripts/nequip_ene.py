
from ase import io 
from ase.units import eV,mol,kcal,Hartree
from nequip.ase.nequip_calculator import NequIPCalculator
import numpy as np
from sklearn.metrics import median_absolute_error
from sklearn.metrics import r2_score
from scipy import stats
from pathlib import Path
import pandas as pd
import os
from ase.io import write, read
from fullspace import tools
from ase.calculators.singlepoint import SinglePointCalculator
from ase.io import write, read
from ase import Atoms

file="/pubhome/xtzhang/interaction/FF/tmp.xyz"
asef = read(file)
mol_name = Path(file).stem

# orca_file = f"/pubhome/xtzhang/interaction/data/CSD/{central}/{central}_{contact}_orca/{mol_name}.json.gz"
# data= tools.json_gzip_load(orca_file)
# true_energy = data[0]["energy"] - data[1]["energy"] - data[2]["energy"]
# 把au转为kcal/mol
# true_energy = true_energy * Hartree/(kcal/mol)
path_model_list = "/pubhome/xtzhang/interaction/scripts/para/ACET"
path_model_list = [os.path.join(path_model_list, f) for f in os.listdir(path_model_list)]
total_pred_list = []
# 从xyz读取真实能量
for path_model in path_model_list:
    predicted_energy_list = []
    calc_neq = NequIPCalculator.from_deployed_model(model_path = path_model, 
                                                    species_to_type_name = {"C" : "C",
                                                                            "H" : "H",
                                                                            "N" : "N",
                                                                            "O" : "O",
                                                                            "S" : "S" },
                                                    device='cpu',
                                                    energy_units_to_eV = (kcal/mol)/eV
                                                    )
    reader_xyz_iter = io.iread(file)
    for _, reader_xyz_item in enumerate(reader_xyz_iter):
        predicted_energy = calc_neq.get_potential_energy(reader_xyz_item)
        predicted_energy_kcal = predicted_energy * eV/(kcal/mol)
        predicted_energy_list.append(predicted_energy_kcal)
    total_pred_list.append(predicted_energy_list)