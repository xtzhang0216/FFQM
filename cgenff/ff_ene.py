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
#%%
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

# get charge, atomtype from str file
# sdffile = '/pubhome/lzeng/data/pair25/rot_split/rot_ACEH_ACEH.sdf' #!!to modify
# sdffile = sys.argv[1]
sdffile = "/pubhome/xtzhang/interaction/FF/data/delCA_ETOH_MBZ.sdf"
len_f = len(Path(sdffile).stem.split('_'))
if len_f == 3:
    _, fga, fgb = Path(sdffile).stem.split('_')
elif len_f == 4:
    _, fga, fgb, idx_f = Path(sdffile).stem.split('_')
fga = "ETOH"
fgb = "MBZ"
print(f'loading {Path(sdffile).name}')
fga_atyps, fga_chgs = GetChgAtype(fga, cwd)
fgb_atyps, fgb_chgs = GetChgAtype(fgb, cwd)

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

fga_sig, fga_eps = AssignSiep(fga_atyps, siep)
fgb_sig, fgb_eps = AssignSiep(fgb_atyps, siep)
#%%
def GetResName(pairname):
    fga_name, fgb_name = pairname.split(':')
    res_a = fga_name.split('-')[-2]
    res_b = fgb_name.split('-')[-2]
    return res_a, res_b
#%%

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
    len_fga = fga_sig.shape[0]
    len_fgb = fgb_sig.shape[0]
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

#%%
aa_resitypes = ['ARG', 'HID', 'HIE', 'HIP', 'LYS', 'ASP','ASH', 'GLU', 'GLH', 'SER', 'THR', 'ASN', 'GLN', 
'CYS', 'GLY', 'PRO', 'ALA', 'VAL', 'ILE', 'LEU', 'MET', 'PHE', 'TYR', 'TRP']
# extract structure data
geo_d = {} 
# {pair_idx:[pairname, [nparray coords of frag1], [nparray coords of frag2]],...} 
# pairs from ligand are not included
pairs_lig = {} # pairs from ligands
idx = 0 # avoid name repetition
with Chem.SDMolSupplier(sdffile, removeHs=False) as suppl:
    for mol in suppl:
        idx += 1
        traj = [mol.GetConformer().GetAtomPosition(a) for a in range(len(mol.GetAtoms()))]
        coords = [[float(pos.x), float(pos.y), float(pos.z)] for pos in traj]
        _, fraga_n, _, _ = mol.GetProp("FRAG_ATOM_NUM").split()
        fraga_n = int(fraga_n)
        pair_idx = f'pair{idx}'
        coords_fa = np.asarray(coords[:fraga_n])
        coords_fb = np.asarray(coords[fraga_n:])
        pairname = mol.GetProp('_Name')
        res_a, res_b = GetResName(pairname)
        # check atom order
        if idx == 1:
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
        # excluding fragments from ligands and record them
        if res_a not in aa_resitypes or res_b not in aa_resitypes:
            pairs_lig[pair_idx] = pairname
            continue
        geo_d[pair_idx] = [pairname, coords_fa, coords_fb]

# js = json.dumps(pairs_lig)
# if len_f == 3:
#     jf = open(f'{Path.home()}/data/pair25/FFData/{fga}_{fgb}_ligand.json', 'w')
# elif len_f == 4:
#     jf = open(f'{Path.home()}/data/pair25/FFData/{fga}_{fgb}_{idx_f}_ligand.json', 'w')
# jf.write(js)
# jf.close()
# %%
# Calculate energy
eng_d = {}
 #{pair_idx:{name:pairname, vdw:value,charge:value,nonbond:value},...}
for p in geo_d:
    eng_d[p] = {}
    pairname = geo_d[p][0]
    fga_coords = geo_d[p][1]
    fgb_coords = geo_d[p][2]
    # print('calculatin vdw ...')
    vdw_array, dist_array = CalVdwDist(fga_sig, fga_eps, fga_coords, fgb_sig, fgb_eps, fgb_coords)
    # print('calculation coulomnb ....')
    chg_array = CalChg(fga_chgs, fgb_chgs, dist_array)
    # kJ/mol to kcal/mol
    vdw_pair = np.sum(vdw_array) / 4.184
    chg_pair = np.sum(chg_array) / 4.184
    nonbond_eng = np.sum([vdw_pair, chg_pair])
    eng_d[p] = {'name': pairname, 'vdw':vdw_pair, 'chg':chg_pair, 'nonbond':nonbond_eng}

js_eng = json.dumps(eng_d)
if len_f == 3:
    jf_eng = open(f'{Path.home()}/interaction/FF/data/{fga}_{fgb}_ffeng.json2','w')
elif len_f == 4:
    jf_eng = open(f'{Path.home()}/interaction/FF/data/{fga}_{fgb}_{idx_f}_ffeng.json2','w')
jf_eng.write(js_eng)
jf_eng.close()
# if len_f == 3:
#     print(f'data are saved in {Path.home()}/data/pair25/FFData/{fga}_{fgb}_ffeng.json')
# elif len_f == 4:
#     print(f'data are saved in {Path.home()}/data/pair25/FFData/{fga}_{fgb}_{idx_f}_ffeng.json')
# %%
