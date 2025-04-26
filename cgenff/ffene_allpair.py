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
#%%
# get sigma, epsilon from ffnonbonded.itp
cwd = Path('/pubhome/xtzhang/interaction/FF/scripts')
# cwd = Path(__file__).absolute().parent
# vdw_params_f = cwd / 'ffnonbonded.itp'
# path = /pubhome/xtzhang/interaction/FF/scripts/charmm36-jul2022.ff
vdw_params_f = Path("/pubhome/xtzhang/interaction/FF/scripts/charmm36-jul2022.ff/ffnonbonded.itp")
siep = {}  #{ffname1:[sig,eps],...}
with open(vdw_params_f, 'r') as f:
    for line in f.readlines():
        if '; The following atom types are NOT part of the CHARMM distribution.' in line:
            break
        line = line.strip()
        # if not line.startswith(';') and not line.startswith('['):
        if len(line) == 0:
            continue
        if line[0].isupper():
            _ = line.split()
            # if len(_) > 1:
            siep[_[0]] = []
            siep[_[0]].append(float(_[5])) # sigma
            siep[_[0]].append(float(_[6])) # eps
# Read energy json
energy_json_path = Path('/pubhome/xtzhang/interaction/data/pdbpairs/total_inteng_comb.json')
# central = "PRPA"
with open(energy_json_path, 'r') as energy_file:
    qm_energy_data = json.load(energy_file)
qm_energy_key = qm_energy_data.keys()
#%%
# qm_energy_data
# # Filter energy data to keep only those starting with 'central'
# filtered_energy_data = {k: v for k, v in qm_energy_data.items() if k.startswith(central)}
# ene = filtered_energy_data
fg_file = {}
files_path = Path('/pubhome/xtzhang/interaction/data/pdbpairs/files.txt')
with open(files_path, 'r') as files_file:
    files = files_file.read().splitlines()
for key in qm_energy_data.keys():
    key = key.replace('.xyz', '')
    fga = key.split('_')[0]
    fgb = key.split('_')[1]
    pair = f'{fga}_{fgb}'
    sdfpath = next((file for file in files if key in file), None)

    if sdfpath and pair in sdfpath:
        if pair in fg_file.keys():
            fg_file[pair].append(sdfpath)
        else:
            fg_file[pair] = [sdfpath]
    else:
        print(f'No file found for {key}')
#%%
aa_resitypes = ['ARG', 'HID', 'HIE', 'HIP', 'LYS', 'ASP','ASH', 'GLU', 'GLH', 'SER', 'THR', 'ASN', 'GLN', 
'CYS', 'GLY', 'PRO', 'ALA', 'VAL', 'ILE', 'LEU', 'MET', 'PHE', 'TYR', 'TRP']
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
def GetChgAtype(fragname, cwd):
    """get charge and atomtype from str file according to fragment"""
    if 'NMA' in fragname:
        fragname = 'NMA'
    chg_fpath = cwd / 'ChgParam' / f'{fragname}.str'
    atomtypes_l = []
    chg_l = []
    with open(chg_fpath, 'r') as f:
        for line in f.readlines():
            if line.startswith('ATOM'):
                _ = line.split()
                atomtypes_l.append(_[2])
                chg_l.append(float(_[3]))
    return atomtypes_l, chg_l

# assign sigma, epsilon value to atoms
def AssignSiep(lista, siep):
    sig_list = []
    eps_list = []
    for k in lista:
        sig_value = float(siep[k][0])
        eps_value = float(siep[k][1])
        sig_list.append(sig_value)
        eps_list.append(eps_value)
    return sig_list, eps_list


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
# get charge, atomtype from str file
# sdffile = '/pubhome/lzeng/data/pair25/rot_split/rot_ACEH_ACEH.sdf' #!!to modify
# sdffile = sys.argv[1]
# sdffile = "/pubhome/lzeng/data/pair25/ProximityEffect/ACEM_ACEM_00_noproxim.sdf"
# len_f = len(Path(sdffile).stem.split('_'))
# if len_f == 3:
#     _, fga, fgb = Path(sdffile).stem.split('_')
# elif len_f == 4:
#     _, fga, fgb, idx_f = Path(sdffile).stem.split('_')
# fga = "ACEM"
# fgb = "ACEM"

# List all files in the directory
calculated_pairs = os.listdir('/pubhome/xtzhang/interaction/FF/data/charmm_ene')
calculated_pairs = [pair.split('.')[0] for pair in calculated_pairs]

# # Print the calculated pairs
# print("Calculated pairs:")
# for pair in calculated_pairs:
#     print(pair)
for pair, sdffiles in fg_file.items():
    # if "HOH" in pair:
    #     continue

    ref_fga, ref_fgb = pair.split('_')
    if 'NMA' in ref_fga:
        ref_fga = 'NMA'
    if 'NMA' in ref_fgb:
        ref_fgb = 'NMA'
    pair = f'{ref_fga}_{ref_fgb}'
    if pair in calculated_pairs:
        continue    
    ref_fga_atyps, ref_fga_chgs = GetChgAtype(ref_fga, cwd)
    ref_fgb_atyps, ref_fgb_chgs = GetChgAtype(ref_fgb, cwd)
    ref_fga_sig, ref_fga_eps = AssignSiep(ref_fga_atyps, siep)
    ref_fgb_sig, ref_fgb_eps = AssignSiep(ref_fgb_atyps, siep)

    for sdffile in sdffiles:
        print(f'loading {Path(sdffile).name}')
        ene = qm_energy_data[f'{Path(sdffile).stem.replace("_noproxim","").replace("delCA_","").replace("delNMA_","")}.xyz']
        # if "ACET_ETOH_07_noproxim.sdf" not in sdffile:
        #     continue


        geo_d = {} 
        # Write molecule to SDF
        # writer = Chem.SDWriter(f'/pubhome/xtzhang/interaction/FF/data//{Path(sdffile).stem}_{idx}.sdf')
        # writer.write(mol)
        # writer.close()
        # {pair_idx:[pairname, [nparray coords of frag1], [nparray coords of frag2]],...} 
        # pairs from ligand are not included
        results = {}
        pairs_lig = {} # pairs from ligands
        with Chem.SDMolSupplier(sdffile, removeHs=False) as suppl:
            for idx, mol in enumerate(suppl):
                _, fraga_n, _, fragb_n = mol.GetProp("FRAG_ATOM_NUM").split()
                fga, fgb = mol.GetProp("FRAG_NAME").split()
                # if fraga_n == '8' or fga == 'ACET':

                if fga == ref_fga and fgb == ref_fgb:
                    fga_atyps, fga_chgs, fga_sig, fga_eps, fgb_atyps, fgb_chgs, fgb_sig, fgb_eps = \
                        ref_fga_atyps, ref_fga_chgs, ref_fga_sig, ref_fga_eps, ref_fgb_atyps, ref_fgb_chgs, ref_fgb_sig, ref_fgb_eps
                elif fga == ref_fgb and fgb == ref_fga:
                    fga_atyps, fga_chgs, fga_sig, fga_eps, fgb_atyps, fgb_chgs, fgb_sig, fgb_eps = \
                        ref_fgb_atyps, ref_fgb_chgs, ref_fgb_sig, ref_fgb_eps, ref_fga_atyps, ref_fga_chgs, ref_fga_sig, ref_fga_eps
                else:
                    print(f'Fragment name not match !')
                    print(f'fga: {fga} fgb: {fgb}')
                    continue
                fraga_n = int(fraga_n)
                if idx == 0:
                    symbols = [i.GetSymbol() for i in mol.GetAtoms()]
                    symb_fa = symbols[:fraga_n]
                    symb_fb = symbols[fraga_n:]
                    symb_fa_str = [i[0] for i in fga_atyps]
                    symb_fb_str = [j[0] for j in fgb_atyps]
                    if symb_fa != symb_fa_str or symb_fb != symb_fb_str:
                        print("Atom order not match !")
                        print(f'symbol of {fga} in sdf: {symb_fa}\nin str:{symb_fa_str}')
                        print(f'symbol of {fgb} in sdf: {symb_fb}\nin str:{symb_fb_str}')
                        exit(0)
                if str(idx) not in ene.keys():
                    continue
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
                results[str(idx)] = [ff_ene, qm_ene]
        # Save results to a CSV file
        results_df = pd.DataFrame.from_dict(results, orient='index', columns=['FF_Energy', 'QM_Energy'])
        results_df.index.name = 'Index'
        output_file = Path(f'/pubhome/xtzhang/interaction/FF/data/charmm_ene/{pair}.csv')
        results_df.to_csv(output_file, mode='a', header=not output_file.exists())
        # print(f'Results saved to {output_file}')
#%%
