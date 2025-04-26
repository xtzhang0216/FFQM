#%%
from __future__ import print_function
import argparse
import sys
import os

from omm_readinputs import *
from omm_readparams import *
from omm_vfswitch import *

from simtk.unit import *
from simtk.openmm import *
from simtk.openmm.app import *

from rdkit import Chem
import pandas as pd
import json
from pathlib import Path

parser = argparse.ArgumentParser()
parser.add_argument("-i", dest="inpfile", help="Input parameter file", required=True)
args = parser.parse_args(
    "-i /pubhome/xtzhang/interaction/FF/Drude/run_drude/md_scripts/cmx.inp".split()
)


#%%


class fragmol:
    def __init__(self, mol, resname, segid, het):
        self.mol = mol
        self.resname = resname
        self.segid = segid
        self.n_atoms = mol.GetNumAtoms()
        self.het = het
        self.ref_mol = Chem.MolFromPDBFile(f"/pubhome/xtzhang/interaction/FF/drude/ref_mols/{resname}.pdb",removeHs=False)
        with open(f"/pubhome/xtzhang/interaction/FF/drude/ref_mols/{resname}.pdb") as f:
            self.drudename = f.readline().split()[3]
    def assign_crd(self):
        match = self.ref_mol.GetSubstructMatch(self.mol)
        for atom in match:
            pos = self.mol.GetConformer().GetAtomPosition(atom)
            self.ref_mol.GetConformer().SetAtomPosition(atom, pos)
    def write_crd(self, crd_file):
        self.assign_crd()
        atoms = [f"{self.n_atoms:>10}  EXT"]
        for i, atom in enumerate(self.ref_mol.GetAtoms(), start=1):
            pos = self.ref_mol.GetConformer().GetAtomPosition(atom.GetIdx())
            atom_name = atom.GetMonomerInfo().GetName()
            line = (
                f"{i:>10}{1:>10}{self.drudename:>6}{atom_name:>9}"
                f"{pos.x:>25.9f}{pos.y:>20.9f}{pos.z:>20.9f}  "
                f"{self.het}      1               0.0000000000"
            )
            atoms.append(line)
        with open(crd_file, "w") as f:
            f.write("\n".join(atoms))
#%%
results = {}
pair="ETOH_MBZ"
sdfpath = "/pubhome/lzeng/data/pair25/del_CA_THR_ETOH/delCA_ETOH_MBZ_01_noproxim.sdf"
sdf = Chem.SDMolSupplier(sdfpath, removeHs=False)
mol = sdf[0]

import parmed as pm

# 假设已通过 Modeller 删除分子 A，得到仅含分子 B 的拓扑和坐标
# fraga_system = Modeller(psf.topology, pdb.positions)
# fraga_system.delete(fragb_idx)
new_topology = fraga_system.topology
new_positions = fraga_system.positions

# ------------------------------------------------
# 关键步骤：将 OpenMM 的 Topology 转换为 ParmEd Structure
# ------------------------------------------------

# 1. 创建 ParmEd 的 Structure 对象
struct = pm.openmm.load_topology(
    topology=new_topology,        # 新拓扑（仅含分子 B）
    positions=new_positions       # 新坐标
)

# 2. 从原始 PSF 中提取分子 B 的参数（原子类型、电荷等）
# 假设原始 PSF 已加载到 ParmEd 的原始结构 orig_struct
orig_struct = pm.load_file('original.psf')  # 读取原始 PSF

# 筛选分子 B 的原子索引（例如根据残基名）
b_indices = [i for i, atom in enumerate(orig_struct.atoms)
             if atom.residue.name == 'B']  # 替换为实际筛选条件

# 将原始参数复制到新 Structure
for idx in b_indices:
    orig_atom = orig_struct.atoms[idx]
    # 找到新 Structure 中对应的原子（通过名称和残基索引匹配）
    new_atom = next(atom for atom in struct.atoms
                    if atom.name == orig_atom.name
                    and atom.residue.number == orig_atom.residue.number)
    # 复制原子类型、电荷等属性
    new_atom.type = orig_atom.type
    new_atom.charge = orig_atom.charge

# 3. 保存为 PSF 和 PDB 文件
struct.save('moleculeB.psf')  # 保存 PSF
struct.save('moleculeB.pdb')  # 保存 PDB

# 可选：检查 PSF 完整性
struct = pm.load_file('moleculeB.psf')  # 重新加载验证
print(struct)