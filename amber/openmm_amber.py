#%%
import os
import shutil
from rdkit import Chem
from rdkit.Chem import rdmolfiles
import subprocess
import pandas as pd
import json
from pathlib import Path
import sys
from openmm.app import AmberPrmtopFile, AmberInpcrdFile
from openmm import System, Context, Platform
from openmm import unit
from openmm.app import NoCutoff,PDBFile
from openmm import VerletIntegrator
import argparse
import tempfile
import time
#%%

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

def write_leap_cmx(fga,fgb):
# source leaprc.protein.ff14SB
# loadamberparams etoh.frcmod
# loadamberparams tolu.frcmod
# ETOH = loadmol2 etoh.mol2
# TOLU = loadmol2 tolu.mol2
# COMPLEX = combine {ETOH TOLU}
# saveamberparm COMPLEX complex.prmtop complex.inpcrd
# quit
    with open("leap.in", "w") as f:
        f.write("source leaprc.protein.ff14SB\n")
        f.write("source leaprc.water.tip3p\n")
        f.write(f"loadamberparams {fga}_a.frcmod\n")
        f.write(f"loadamberparams {fgb}_b.frcmod\n")
        f.write(f"{fga}_a = loadmol2 {fga}_a.mol2\n")
        f.write(f"{fgb}_b = loadmol2 {fgb}_b.mol2\n")
        f.write(f"COMPLEX = combine {{{fga}_a {fgb}_b}}\n")
        f.write(f"saveamberparm COMPLEX {fga}_{fgb}.prmtop {fga}_{fgb}.inpcrd\n")
        f.write("quit\n")
    subprocess.run(["/pubhome/xtzhang/mambaforge/envs/htmd/bin/tleap", "-f", "leap.in"])
        
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

def calc_ene(top_file, crd_file):
    # Load the topology and coordinates
    psf = AmberPrmtopFile(top_file)
    pdb = AmberInpcrdFile(crd_file)
    system = psf.createSystem(
        nonbondedMethod=NoCutoff,
        constraints=None,
        implicitSolvent=None
    )
    integrator = VerletIntegrator(0.001 * unit.picoseconds)
    context = Context(system, integrator)
    context.setPositions(pdb.positions)
    state = context.getState(getEnergy=True)
    potential_energy = state.getPotentialEnergy().value_in_unit(unit.kilocalories_per_mole)
    return potential_energy



def run_sander(prefix):
    subprocess.run([
        "/pubhome/xtzhang/mambaforge/envs/htmd/bin/sander", "-O", "-i", "min.in", "-o", f"{prefix}.out",
        "-p", f"{prefix}.prmtop", "-c", f"{prefix}.inpcrd",
        "-r", f"{prefix}.rst", "-ref", f"{prefix}.inpcrd"
    ])
    with open(f"{prefix}.out") as f:
        for line in f:
            energy_line = None
            for line in f:
                if "ENERGY" in line:
                    energy_line = next(f, None)
            if energy_line:
                return float(energy_line.split()[1])
    return None


def check_nh_distance(fragA, fragB, fragA_name, idx):
    if fragA_name not in ["ACEM", "MIMM", "MIME", "MIMD"]:
        return True
    for atomA in fragA.GetAtoms():
        if atomA.GetSymbol() == 'N':
            for neighbor in atomA.GetNeighbors():
                if neighbor.GetSymbol() == 'H':
                    posH = fragA.GetConformer().GetAtomPosition(neighbor.GetIdx())
                    for atomB in fragB.GetAtoms():
                        posB = fragB.GetConformer().GetAtomPosition(atomB.GetIdx())
                        dist = ((posH.x-posB.x)**2 + (posH.y-posB.y)**2 + (posH.z-posB.z)**2)**0.5
                        if dist < 1.5:
                            print(f"N-H distance too short between atom {neighbor.GetIdx()} (NH) and atom {atomB.GetIdx()} ({atomB.GetSymbol()}): {dist:.2f} Å in conformer {idx}")
                            # Save fraga and fragb as SDF for inspection
                            # Chem.SDWriter(f"{fragA_name}_a_conf{idx}.sdf").write(fragA)
                            # Chem.SDWriter(f"{fragB_name}_b_conf{idx}.sdf").write(fragB)
                            return False
    return True


def main():
    parser = argparse.ArgumentParser(description='Amber OpenMM energy calculation')
    parser.add_argument('--xyz', type=str, default="ACEM_ACEM_00",help='xyz name, e.g. ACEM_HOH_00')
    parser.add_argument('--start', type=int, default=0, help='start index for sdf conformers')
    parser.add_argument('--end', type=int, default=100, help='end index for sdf conformers (exclusive)')
    parser.add_argument('--check_nh', type=bool, default=False, help='check N-H distance for specific molecules')
    args = parser.parse_args()

    xyz = args.xyz
    start = args.start
    end = args.end

    energy_json_path = Path('/pubhome/lzeng/data/pair25/NequipData/TOTAL/total_inteng_comb.json')
    with open(energy_json_path, 'r') as energy_file:
        qm_energy_data = json.load(energy_file)
    ene = qm_energy_data[f"{xyz}.xyz"]
    results = {}
    with open ("/pubhome/xtzhang/interaction/data/pdbpairs/files.txt", "r") as f:
        lines = f.readlines()
    sdfpath = next((line for line in lines if xyz in line), None)
    sdfpath = sdfpath.strip()
    sdf = Chem.SDMolSupplier(sdfpath, removeHs=False)

    WORKDIR = tempfile.mkdtemp(prefix=f"{xyz}_")
    total = len(sdf)
    if end is None or end > total:
        end = total
    for idx in range(start, end):
        mol = sdf[idx]
        if str(idx) not in ene.keys():
            continue
        if mol is None:
            continue
        fga_name, fgb_name = mol.GetProp("FRAG_NAME").split()
        fraga, fragb = Chem.GetMolFrags(mol, asMols=True, sanitizeFrags=True)
        print(f"Processing conformer {idx} / {total}")
        conf_dir = os.path.join(WORKDIR, f"conf_{idx}")
        if fga_name == "ACEH" or fgb_name == "ACEH":
            continue
        if fga_name == "ACET" and fraga.GetNumAtoms() != 7:
            continue
        if fgb_name == "ACET" and fragb.GetNumAtoms() != 7:
            continue

        if args.check_nh:
            if not check_nh_distance(fraga, fragb, fga_name, idx) or not check_nh_distance(fragb, fraga, fgb_name, idx):
                continue
        os.makedirs(conf_dir, exist_ok=True)
        os.chdir(conf_dir)
        mapping(fga_name, fraga,'a')
        mapping(fgb_name, fragb,'b')
        shutil.copy(f"/pubhome/xtzhang/interaction/FF/amber/parameter/{fga_name}.frcmod", f"{conf_dir}/{fga_name}_a.frcmod")
        shutil.copy(f"/pubhome/xtzhang/interaction/FF/amber/parameter/{fgb_name}.frcmod", f"{conf_dir}/{fgb_name}_b.frcmod")
        write_leap(f"{fga_name}_a")
        write_leap(f"{fgb_name}_b")
        write_leap_cmx(fga_name, fgb_name)
        e_a = calc_ene(f"{fga_name}_a.prmtop", f"{fga_name}_a.inpcrd")
        e_b = calc_ene(f"{fgb_name}_b.prmtop", f"{fgb_name}_b.inpcrd")
        e_complex = calc_ene(f"{fga_name}_{fgb_name}.prmtop", f"{fga_name}_{fgb_name}.inpcrd")
        interaction_energy = e_complex - (e_a + e_b)
        qm_ene = ene[str(idx)]
        results[str(idx)] = [round(interaction_energy, 2), round(qm_ene, 2)]
        os.chdir("../")
        shutil.rmtree(conf_dir)
    results_df = pd.DataFrame.from_dict(results, orient='index', columns=['FF_Energy', 'QM_Energy'])
    results_df.index.name = 'Index'
    if fgb_name[0] < fga_name[0]:
        pair = f"{fgb_name}_{fga_name}"
    else:
        pair = f"{fga_name}_{fgb_name}"
    output_file = Path(f'/pubhome/xtzhang/interaction/FF/amber/data/{pair}.csv')
    results_df.to_csv(output_file, mode='a', header=not output_file.exists())

if __name__ == "__main__":
    main()
# %%
# rot_sdf_path = next((path for path in Path('/pubhome/lzeng/data/pair25/rot_split/').glob(f'*{xyz}*')), None)
# if rot_sdf_path:
#     rot_sdf = Chem.SDMolSupplier(str(rot_sdf_path), removeHs=False)
# else:
#     print(f"No rotated SDF file found for {xyz}")
#     rot_sdf = None

# rot_ene = {k: v for k, v in ene.items() if k.startswith('rot')}
# total = len(rot_sdf)
# if rot_ene:
#     for idx, mol in enumerate(rot_sdf):
#         # 
#         # start_time = time.time()
#         # with open(f"{WORKDIR}/log.txt", "a") as log_file:
#         #     log_file.write(f"Processing conformer {idx + 1} / {total}\n")
#         # if idx > 500:
#         #     break
#         if f"rot-{idx}" not in rot_ene.keys():
#             continue
#         if mol is None:
#             continue
#         fga_name, fgb_name = mol.GetProp("FRAG_NAME").split()

#         print(f"Processing conformer {idx + 1} / {total}")

#         conf_dir = os.path.join(WORKDIR, f"conf_{idx}")
#         os.makedirs(conf_dir, exist_ok=True)
#         os.chdir(conf_dir)

#         fraga, fragb = Chem.GetMolFrags(mol, asMols=True, sanitizeFrags=True)
#         if fga_name == "ACEH" or fgb_name == "ACEH":
#             continue
#         if fga_name == "ACET" and fraga.GetNumAtoms() != 7:
#             continue
#         if fgb_name == "ACET" and fragb.GetNumAtoms() != 7:
#             continue
#         # get mol2
#         mapping(fga_name, fraga,'a')
#         mapping(fgb_name, fragb,'b')
#         # get frcmod
#         shutil.copy(f"/pubhome/xtzhang/interaction/FF/amber/parameter/{fga_name}.frcmod", f"{conf_dir}/{fga_name}_a.frcmod")
#         shutil.copy(f"/pubhome/xtzhang/interaction/FF/amber/parameter/{fgb_name}.frcmod", f"{conf_dir}/{fgb_name}_b.frcmod")
#         # with Chem.SDWriter(os.path.join(conf_dir, f"{fga_name}.sdf")) as writer:
#         #     writer.write(fragb)
#         # with Chem.SDWriter(os.path.join(conf_dir, f"{fgb_name}.sdf")) as writer:
#         #     writer.write(fraga)
#         # with Chem.SDWriter(os.path.join(conf_dir, f"{fga_name}_{fgb_name}.sdf")) as writer:
#         #     writer.write(mol)
#         # run_antechamber(f"{fga_name}.sdf", f"{fga_name}")
#         write_leap(f"{fga_name}_a")
#         # run_antechamber(f"{fgb_name}.sdf", f"{fgb_name}")
#         write_leap(f"{fgb_name}_b")
#         # complex
#         # run_antechamber(f"{fga_name}_{fragb_name}.sdf",   f"{fga_name}_{fgb_name}")
#         write_leap_cmx(fga_name, fgb_name)
        
#         # Energy calculation
#         # write_min_input()
#         e_a = calc_ene(f"{fga_name}_a.prmtop", f"{fga_name}_a.inpcrd")
#         e_b = calc_ene(f"{fgb_name}_b.prmtop", f"{fgb_name}_b.inpcrd")
#         e_complex = calc_ene(f"{fga_name}_{fgb_name}.prmtop", f"{fga_name}_{fgb_name}.inpcrd")
#         # e_a = run_sander(f"{fga_name}_a")
#         # e_b = run_sander(f"{fgb_name}_b")
#         # e_complex = run_sander(f"{fga_name}_{fgb_name}")

#         interaction_energy = e_complex - (e_a + e_b)
#         # interaction_energy_kcal = interaction_energy/4.184 
#         qm_ene = rot_ene[f"rot-{idx}"]
#         # save .2f
#         results[f"rot-{idx}"] = [round(interaction_energy, 2), round(qm_ene, 2)]

#         os.chdir("../")
#         shutil.rmtree(conf_dir)
# os.chdir("/")
# shutil.rmtree(WORKDIR)
# results_df = pd.DataFrame.from_dict(results, orient='index', columns=['FF_Energy', 'QM_Energy'])
# results_df.index.name = 'Index'
# if fgb_name[0] < fga_name[0]:
#     pair = f"{fgb_name}_{fga_name}"
# else:
#     pair = f"{fga_name}_{fgb_name}"
# output_file = Path(f'/pubhome/xtzhang/interaction/FF/amber/data/{pair}.csv')
# results_df.to_csv(output_file, mode='a', header=not output_file.exists())
