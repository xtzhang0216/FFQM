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
import shutil
import subprocess
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
def write_cmx_pdb(mol1,mol2,pdb_file):
    length1 = mol1.GetNumAtoms()
    length2 = mol2.GetNumAtoms()
    combined = Chem.CombineMols(mol1, mol2)
    Chem.MolToPDBFile(combined, pdb_file)
    with open(pdb_file, 'r') as f:
        lines = f.readlines()

    with open(pdb_file, 'w') as f:
        for i,line in enumerate(lines):
            if line.startswith("ATOM") or line.startswith("HETATM"):
                if i < length1:
                    f.write(line[:25] + "1" + line[26:])
                elif length1 <= i < length1 + length2:
                    f.write(line[:25] + "2" + line[26:])
            else:
                f.write(line)

def mapping(frag, query,id):
    ref_mol_folder = "/pubhome/xtzhang/interaction/FF/amber/parameter"

    with open(os.path.join(ref_mol_folder, f"{frag}.mol2"), "r") as f:
        line1 = f.readlines()
    # Modify coordinates in line1 based on uery):
    crd_start=8
    for i in range(query.GetNumAtoms()):
        atom = query.GetAtomWithIdx(i)
        query_symbol = atom.GetSymbol()
        ref_symbol = line1[crd_start + i].split()[1][0]
        if query_symbol != ref_symbol:
            raise ValueError(f"Atom symbol mismatch: {query_symbol} != {ref_symbol}")
        crd = query.GetConformer().GetAtomPosition(i)
        line_idx = crd_start + i
        line1[line_idx] = line1[line_idx][:16] + f"{crd[0]:11.4f}" + line1[line_idx][27:]
        line1[line_idx] = line1[line_idx][:27] + f"{crd[1]:11.4f}" + line1[line_idx][38:]
        line1[line_idx] = line1[line_idx][:38] + f"{crd[2]:11.4f}" + line1[line_idx][49:]
    with open(f"{frag}_{id}.mol2", "w") as f:
        f.writelines(line1)

def calc_ene(crd_file):
    # Load the topology and coordinates
    # psf = ForceField("/pubhome/xtzhang/interaction/FF/charmm/charmm36.zxt.xml")
    psf = ForceField('/pubhome/xtzhang/interaction/FF/charmm/charmm36.zxt.xml','charmm36/water.xml')

    pdb = PDBFile(crd_file)
    # for res in pdb.topology.residues():
    #     print(f"PDB Residue: {res.name}")
    system = psf.createSystem(
        pdb.topology,
        nonbondedMethod=NoCutoff,
        constraints=None,
    )
    integrator = VerletIntegrator(0.001 * picoseconds)  # 步长设为极小值（不影响单点能量计算）
    context = Context(system, integrator)
    context.setPositions(pdb.positions)
    state = context.getState(getEnergy=True)
    potential_energy = state.getPotentialEnergy().value_in_unit(kilocalories_per_mole)
    return potential_energy

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

def get_opls(pdb_file):
    # pdb_file = "/pubhome/xtzhang/interaction/FF/mol2/ACEM.pdb"
    psf = ForceField('/pubhome/xtzhang/interaction/FF/opls/param/opls.zxt.xml','tip3p.xml')
    pdb = PDBFile(pdb_file)
    # for res in pdb.topology.residues():
    #     print(f"PDB Residue: {res.name}")
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

def get_charmm(pdb_file):
    # pdb_file = "/pubhome/xtzhang/interaction/FF/mol2/ACEM.pdb"
    psf = ForceField('/pubhome/xtzhang/interaction/FF/outliar/charmm/charmm36.zxt.fake.xml','charmm36/water.xml')
    pdb = PDBFile(pdb_file)
    # for res in pdb.topology.residues():
    #     print(f"PDB Residue: {res.name}")
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

def get_amber(top_file, crd_file):
    psf = AmberPrmtopFile(top_file)
    pdb = AmberInpcrdFile(crd_file)
    system = psf.createSystem(
        nonbondedMethod=NoCutoff,
        constraints=None,
        implicitSolvent=None
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

def write_leap(prefix):
    with open("leap.in", "w") as f:
        f.write("source leaprc.protein.ff14SB\n")
        f.write("source leaprc.water.tip3p\n")
        f.write(f"LIG = loadmol2 {prefix}.mol2\n")
        f.write(f"loadamberparams {prefix}.frcmod\n")
        f.write(f"saveamberparm LIG {prefix}.prmtop {prefix}.inpcrd\n")# generate topology and coordinate for specific conformation 
        f.write("quit\n")
        f.write("quit\n")
    subprocess.run(["/pubhome/xtzhang/mambaforge/envs/htmd/bin/tleap", "-f", "leap.in"])


# xyz="ACEM_ACET_00"
results = {}
model = sys.argv[1]
pair = sys.argv[2]
# model = "charmm"
# pair = "ACET_PRPA"
sdfpath = f"/pubhome/xtzhang/interaction/FF/outliar/tail/{pair}/{pair}_outliar.charmm.sdf"
sdf = Chem.SDMolSupplier(sdfpath, removeHs=False)
print(f"Total {len(sdf)} mols")
for idx, mol in enumerate(sdf):
    # if str(idx) not in ene.keys():
    #     continue

    mol_index = mol.GetProp("Index") if mol.HasProp("Index") else str(idx)
    # mol_index = str(idx)
    print(f"Processing {idx}")
    # if mol_index != "695":
    #     continue
    os.chdir("/pubhome/xtzhang/interaction/FF/tmp")
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
    # e_a = calc_ene(f"{fga_name}_a.pdb")
    # e_b = calc_ene(f"{fgb_name}_b.pdb")
    # write_cmx_pdb(fraga, fragb, f"{fga_name}_{fgb_name}.pdb")
    # e_complex = calc_ene(f"{fga_name}_{fgb_name}.pdb")
    # interaction_energy = e_complex - (e_a + e_b)
    if model == "opls":
        fga_chgs,fga_sig,fga_eps = get_opls(f"{fga_name}_a.pdb")
        fgb_chgs,fgb_sig,fgb_eps = get_opls(f"{fgb_name}_b.pdb")
    elif model == "charmm":
        fga_chgs,fga_sig,fga_eps = get_charmm(f"{fga_name}_a.pdb")
        fgb_chgs,fgb_sig,fgb_eps = get_charmm(f"{fgb_name}_b.pdb")
    elif model == "amber":
        shutil.copy(f"/pubhome/xtzhang/interaction/FF/amber/parameter/{fga_name}.frcmod", f"./{fga_name}_a.frcmod")
        shutil.copy(f"/pubhome/xtzhang/interaction/FF/amber/parameter/{fgb_name}.frcmod", f"./{fgb_name}_b.frcmod")
        mapping(fga_name, fraga,'a')
        mapping(fgb_name, fragb,'b')
        write_leap(f"{fga_name}_a")
        write_leap(f"{fgb_name}_b")
        fga_chgs,fga_sig,fga_eps = get_amber(f"{fga_name}_a.prmtop", f"{fga_name}_a.inpcrd")
        fgb_chgs,fgb_sig,fgb_eps = get_amber(f"{fgb_name}_b.prmtop", f"{fgb_name}_b.inpcrd")
        

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
    qm_ene = float(mol.GetProp("QM_Energy"))
    ff_supermolecule = float(mol.GetProp("FF_Energy"))
    # 使用分子的Index属性作为结果字典的键
    results[mol_index] = [round((qm_ene-ff_supermolecule), 2), round(qm_ene, 2), round(ff_supermolecule, 2), round(ff_ene, 2), round(vdw_pair, 2), round(chg_pair, 2)]
# Save results to a CSV file
results_df = pd.DataFrame.from_dict(results, orient='index', columns=['diff', 'QM_Energy', 'FF_Supermolecule_Energy', 'Vdw_Chg', 'Vdw_Energy', 'Chg_Energy'])
results_df.index.name = 'Index'
if fgb_name[0] < fga_name[0]:
    pair = f"{fgb_name}_{fga_name}"
else:
    pair = f"{fga_name}_{fgb_name}"
output_file = Path(f'/pubhome/xtzhang/interaction/FF/outliar/tail/{pair}/{pair}.{model}.csv')
results_df.to_csv(output_file, mode='a', header=not output_file.exists())
# print(f'Results saved to {output_file}')
#%%
