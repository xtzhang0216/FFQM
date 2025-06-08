#!/bin/env python
"""calculate non-bond interaction including Vdw & Coulomb force with charmm36 parameters
17 fragment str files , ffnonbonded.itp are needed
get charge, atomtype, sigma, epsilon, from str files and ffnonbonded.itp
get coordination from structure files

run under conda env py37
usage:
python xxx.py input_file_path
"""
#%%
from pathlib import Path
from rdkit import Chem
import numpy as np
import json
import sys
import csv
import pandas as pd
import os
FragAtomNum = {
    "ACEM": 9,
    "PCEM": 12,
    "ACET": 7,
    "ACEH": 8,
    "MBZ": 15,
    "MIMD": 12,
    "MIME": 12,
    "MIMM": 13,
    "MIND": 19,
    "ETAM": 11,
    "ETOH": 9,
    "ETSH": 9,
    "MGDM": 13,
    "MSM": 9,
    "NMA": 12,
    "MPHE": 16,
    "PRPA": 11,
    "N1PA": 19,
    # "HOH" : 3,
    "BUT": 14, 
    "IBUT": 14,
    "NMP": 15,
    "NEA": 15,
    "MTA": 5,
    "ETA": 8,
}

    # fg_file[key] = next((file for file in files if fga in file and fgb in file), None)

def CalDist(coord_i, coords_j):
    """Calculate distance between atom_i from fragment a and atoms from fragment b"""
    coord_i = np.asarray(coord_i)
    coords_j = np.asarray(coords_j)
    dists_array = np.sqrt(np.sum(((coords_j - coord_i) ** 2), axis=1))
    return dists_array

def CalVdwDist(fga_sig, fga_eps, fga_coords, fgb_sig, fgb_eps, fgb_coords):
    '''calculate vdw between fragment a and fragment b output values in array'''
    fga_sig = np.asarray(fga_sig)
    fgb_sig = np.asarray(fgb_sig)
    fga_eps = np.asarray(fga_eps)
    fgb_eps = np.asarray(fgb_eps)
    # len_fga = fga_sig.shape[0]
    # len_fgb = fgb_sig.shape[0]
    # sig_ijs_array = np.zeros((len_fga, len_fgb))
    # eps_ijs_array = np.zeros((len_fga, len_fgb))
    # dist_ijs_array = np.zeros((len_fga, len_fgb))
    # for i in range(len_fga):
    #     sig_ijs_array[i] = (fga_sig[i] + fgb_sig) / 2
    #     eps_ijs_array[i] = np.sqrt(fga_eps[i] * fgb_eps)
    #     dist_ijs_array[i] = CalDist(fga_coords[i], fgb_coords)
    sig_ijs_array = (fga_sig[:, np.newaxis] + fgb_sig) / 2
    eps_ijs_array = np.sqrt(fga_eps[:, np.newaxis] * fgb_eps)
    dist_ijs_array = np.linalg.norm(fga_coords[:, np.newaxis] - fgb_coords, axis=2)
    _ = (sig_ijs_array / dist_ijs_array) ** 6
    vdw_ijs_array = 4 * eps_ijs_array * (_ * _ - _)
    return vdw_ijs_array, dist_ijs_array

def CalChg(fga_chgs, fgb_chgs, dist_ijs_array):
    '''calculate coulomb between fragment a and fragment b output values in array'''
    f = 138.935458
    fga_chgs = np.asarray(fga_chgs)
    fgb_chgs = np.asarray(fgb_chgs)
    len_fga = fga_chgs.shape[0]
    len_fgb = fgb_chgs.shape[0]
    # q_ijs_array = np.zeros((len_fga, len_fgb))
    # for i in range(len_fga):
    #     q_ijs_array[i] = fga_chgs[i] * fgb_chgs
    q_ijs_array = fga_chgs[:, np.newaxis] * fgb_chgs
    chg_ijs_array = f * (q_ijs_array / dist_ijs_array)
    chg_ijs_array = f * (q_ijs_array / dist_ijs_array)
    return chg_ijs_array

def GetResName(pairname):
    fga_name, fgb_name = pairname.split(':')
    res_a = fga_name.split('-')[-2]
    res_b = fgb_name.split('-')[-2]
    return res_a, res_b

from openmm.app import *
from openmm import *
from openmm.unit import *

def get(pdb_file):
    # pdb_file = "/pubhome/xtzhang/interaction/FF/mol2/ACEM.pdb"
    psf = ForceField('/pubhome/xtzhang/interaction/FF/opls/param/opls.zxt.xml','tip3p.xml') 

    pdb = PDBFile(pdb_file)
    for res in pdb.topology.residues():
        print(f"PDB Residue: {res.name}")
    system = psf.createSystem(
        pdb.topology,
        nonbondedMethod=NoCutoff,
        constraints=None,
    )
    integrator = VerletIntegrator(0.001 * picoseconds)  # 步长设为极小值（不影响单点能量计算）
    context = Context(system, integrator)
    context.setPositions(pdb.positions)
    nonbonded = [f for f in system.getForces() if isinstance(f, NonbondedForce)][0]
    charges,sigmas,epsilons = [],[],[]
    for i in range(system.getNumParticles()):
        charge, sigma, epsilon = nonbonded.getParticleParameters(i)
        charges.append(charge.value_in_unit(elementary_charge))
        sigmas.append(sigma.value_in_unit(nanometers))
        epsilons.append(epsilon.value_in_unit(kilojoules_per_mole))
    return charges,sigmas,epsilons

xyz="ACEM_ACET_00"

energy_json_path = Path('/pubhome/lzeng/data/pair25/NequipData/TOTAL/total_inteng_comb.json')
with open(energy_json_path, 'r') as energy_file:
    qm_energy_data = json.load(energy_file)
# for key in qm_energy_data.keys():
ene = qm_energy_data[f"{xyz}.xyz"]
results = {}
with open ("/pubhome/xtzhang/interaction/data/pdbpairs/files.txt", "r") as f:
    lines = f.readlines()
sdfpath = next((line for line in lines if xyz in line), None)
sdfpath = sdfpath.strip()
sdf = Chem.SDMolSupplier(sdfpath, removeHs=False)
for idx, mol in enumerate(sdf):
    if str(idx) not in ene.keys():
        continue
    if idx > 500:
        break
    os.chdir("/pubhome/xtzhang/interaction/FF/test")
    fga_name, fgb_name = mol.GetProp("FRAG_NAME").split()
    fraga, fragb = Chem.GetMolFrags(mol, asMols=True, sanitizeFrags=True)
    if fga_name == "ACEH" or fgb_name == "ACEH":
        continue
    if fga_name == "ACET" and fraga.GetNumAtoms() != 7:
        continue
    if fgb_name == "ACET" and fragb.GetNumAtoms() != 7:
        continue
    _, fraga_n, _, fragb_n = mol.GetProp("FRAG_ATOM_NUM").split()
    fga_name, fgb_name = mol.GetProp("FRAG_NAME").split()
    Chem.MolToPDBFile(fraga, f"{fga_name}_a.pdb")
    Chem.MolToPDBFile(fragb, f"{fgb_name}_b.pdb")
    fga_chgs,fga_sig,fga_eps = get(f"{fga_name}_a.pdb")
    fgb_chgs,fgb_sig,fgb_eps = get(f"{fgb_name}_b.pdb")
    fraga_n, fragb_n = int(fraga_n), int(fragb_n)

    traj = [mol.GetConformer().GetAtomPosition(a) for a in range(len(mol.GetAtoms()))]
    coords = [[float(pos.x), float(pos.y), float(pos.z)] for pos in traj]
    coords_fa = np.asarray(coords[:fraga_n])
    coords_fb = np.asarray(coords[fraga_n:])
    # convert to nm
    fga_coords = coords_fa / 10
    fgb_coords = coords_fb / 10
    vdw_array, dist_array = CalVdwDist(fga_sig, fga_eps, fga_coords, fgb_sig, fgb_eps, fgb_coords)
    # print('calculation coulomnb ....')
    chg_array = CalChg(fga_chgs, fgb_chgs, dist_array)
    # kJ/mol to kcal/mol
    vdw_pair = np.sum(vdw_array) / 4.184
    chg_pair = np.sum(chg_array) / 4.184
    ff_ene = np.sum([vdw_pair, chg_pair])
    qm_ene = ene[str(idx)]
    results[str(idx)] = [round(vdw_pair,2),round(chg_pair,2), round(ff_ene,2), round(qm_ene,2)]
# Save results to a CSV file
results_df = pd.DataFrame.from_dict(results, orient='index', columns=['Vdw_Energy', 'Chg_Energy', 'FF_Energy', 'QM_Energy'])
results_df.index.name = 'Index'
if fgb_name[0] < fga_name[0]:
    pair = f"{fgb_name}_{fga_name}"
else:
    pair = f"{fga_name}_{fgb_name}"
output_file = Path(f'/pubhome/xtzhang/interaction/FF/opls/compose/{pair}.csv')
results_df.to_csv(output_file, mode='a', header=not output_file.exists())
# print(f'Results saved to {output_file}')
#%%
