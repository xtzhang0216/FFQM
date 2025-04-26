
#!/pubhome/xtzhang/mambaforge/envs/htmd/bin/python
#%%
from rdkit import Chem
from fullspace import fragment, tools
import argparse
import os

Name2Charge = {
    "ACEM": 0,
    "ACET": -1,
    # "ACEH": 0,
    "MBZ": 0,
    "MIMD": 0,
    "MIME": 0,
    "MIMM": 1,
    "MIND": 0,
    "ETAM": 1,
    "ETOH": 0,
    "ETSH": 0,
    "MGDM": 1,
    "MSM": 0,
    "NMA": 0,
    "MTYR": 0,
    "PRPA": 0,
    "N1PA": 0,
    "HOH" : 0
}


# 从/pubhome/xtzhang/interaction/test.pair.sdf读取mol
def bsse_mol(mol_a, mol_b, charge_a=0, charge_b=0):
    mols = []
    species_fa = list(mol_a.species)
    species_fb = list(mol_b.species)
    new_gto_fa = list(mol_a.new_gto)
    new_gto_fb = list(mol_b.new_gto)
    new_gto = new_gto_fa + new_gto_fb

    if charge_a != -1 and charge_b != -1:
        new_gto_fa = new_gto
        new_gto_fb = new_gto
    elif charge_a == -1:
        diff_fuc = 'newgto "ma-def2-TZVP" end'
        # 给氧原子加弥散函数
        O_idx = [i for i, ele in enumerate(species_fa) if ele == 'O']
        for i in O_idx:
            new_gto[i] = diff_fuc
            new_gto_fa = new_gto
            new_gto_fb = new_gto
    elif charge_b == -1:
        diff_fuc = 'newgto "ma-def2-TZVP" end'
        # 给氧原子加弥散函数
        O_idx = [i for i, ele in enumerate(species_fb) if ele == 'O']
        O_idx = [i + len(species_fa) for i in O_idx]
        for i in O_idx:
            new_gto[i] = diff_fuc
            new_gto_fa = new_gto
            new_gto_fb = new_gto
        

    coord_fa = list(mol_a.refCoord)
    coord_fb = list(mol_b.refCoord)
    charge_fa = mol_a.charge
    charge_fb = mol_b.charge
    fa_num = len(species_fa)
    name_fa = mol_a.name
    name = name_fa.split('_')[0]
    name_fb = mol_b.name
    species = species_fa + species_fb
    coord = coord_fa + coord_fb
    charge = charge_fa + charge_fb
    species_fa = species[:fa_num] + [f'{ele} :' for ele in species[fa_num:]]
    species_fb = [f'{ele} :' for ele in species[:fa_num]] + species[fa_num:]

    coord_fa = coord
    coord_fb = coord
    f = fragment.Mol(
        name,
        species=species,
        refCoord=coord,
        charge=charge,
        new_gto=new_gto,
        reset_com=False,
        )
    mols.append(f)
    fa = fragment.Mol(
            name_fa,
            species=species_fa,
            refCoord=coord_fa,
            charge=charge_fa,
            new_gto=new_gto_fa,
            reset_com=False,
        )
    mols.append(fa)
    fb = fragment.Mol(
            name_fb,
            species=species_fb,
            refCoord=coord_fb,
            charge=charge_fb,
            new_gto=new_gto_fb,
            reset_com=False,
        )
    mols.append(fb)
    return mols

def rdkit2frag(mol, name, charge=0):
    # 从rdkit mol中得到fragment mol
    species = [a.GetSymbol() for a in mol.GetAtoms()]
    new_gto = [" " for i in mol.GetAtoms()]
    traj = [
        mol.GetConformer().GetAtomPosition(a) for a in range(len(mol.GetAtoms()))
    ]
    coord = [[pos.x, pos.y, pos.z] for pos in traj]

    frag_m = fragment.Mol(
            name,
            species=species,
            refCoord=coord,
            charge=charge,
            new_gto=new_gto,
            reset_com=False,
        )
    return frag_m
#%%
if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--sdf", required=True, help="contact smiles"
    )
    args = parser.parse_args(
        ["--sdf", "/pubhome/xtzhang/interaction/FF/scripts/1.sdf"]
    )
    # mols = Chem.rdmolfiles.SDMolSupplier(f'/pubhome/xtzhang/interaction/data/{args.central_smiles}_{args.contact_smiles}.pair.sdf', removeHs=False)
    mol = Chem.rdmolfiles.SDMolSupplier(args.sdf, removeHs=False)[0]

    frag_mols = []
    # name = 'H2O_NH3'
    name = args.sdf.split('/')[-1].split('.')[0]
    # frag_name = mol.GetProp("FRAG")
    # frag_name_1 = frag_name.split('_')[0]
    # frag_name_2 = frag_name.split('_')[1]
    frag_name_1, frag_name_2 = mol.GetProp("FRAG_NAME").split(' ')[0], mol.GetProp("FRAG_NAME").split('   ')[1]
    charge_1 = Name2Charge[frag_name_1]
    charge_2 = Name2Charge[frag_name_2]
    # 寻找是否已经含有该文件
    if os.path.exists(f'/pubhome/xtzhang/interaction/data/csd_pairs/modified_orca/{frag_name_1}_{frag_name_2}/{name}.json.gz'):
        exit(0)
    match_frags = Chem.GetMolFrags(
                mol, asMols=True, sanitizeFrags=True
            )
    frag_mols.extend(bsse_mol(rdkit2frag(match_frags[0], name, charge_1), rdkit2frag(match_frags[1], name, charge_2),charge_1,charge_2) )
    # print(name)
    data = tools.orca_runall(
        mols=frag_mols,
        parse=True,
        setting="! engrad wB97X-D3BJ  def2-TZVPP def2/J def2-TZVPP/C RIJCOSX\n",
        )
    orca_energy = data[0]["energy"] - data[1]["energy"] - data[2]["energy"]
    from ase.units import eV,mol,kcal,Hartree
    orca_energy = orca_energy * Hartree/(kcal/mol)
    print(f'{name} {round(orca_energy,2)}')
    folder = f'/pubhome/xtzhang/interaction/data/csd_pairs/modified_orca/'
    # folder = f'/pubhome/xtzhang/interaction/data/{args.central_smiles}_{args.contact_smiles}_orca'
    if not os.path.exists(folder):
        os.makedirs(folder)
    # num = len([file for file in os.listdir(folder) if file.startswith(name)])
    # print(f'{folder}/{name}_{num}.json.gz')

    # 检查文件.json.gz是否存在，如果已经存在，则其他后缀
    # if os.path.exists(f'{folder}/{name}.json.gz'):
        # 已经同样name的数量

    tools.json_gzip_dump(data,f'{folder}/{frag_name}/{name}.json.gz')