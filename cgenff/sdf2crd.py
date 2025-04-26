from rdkit import Chem

# 读取SDF文件
suppl = Chem.SDMolSupplier('/pubhome/xtzhang/interaction/FF/data/delCA_ETOH_MBZ.sdf')
mol = next(suppl)

# 生成CRD内容
crd_content1 = ["* Generated from SDF\n", f"9\n"]
crd_content2 = ["* Generated from SDF\n", f"15\n"]
for i, atom in enumerate(mol.GetAtoms()):
    if i <= 8:
        pos = mol.GetConformer().GetAtomPosition(i)
        line = f"{i+1:9d}         1  ETOH      {atom.GetSymbol()}{i+1:4d} {atom.GetSymbol():2s}         {pos.x:8.3f}        {pos.y:8.3f}         {pos.z:8.3f} HETA 1 0.0\n"
        crd_content1.append(line)
    else:
        pos = mol.GetConformer().GetAtomPosition(i)
        line = f"{i+1:9d}         1  TOLU      {atom.GetSymbol()}{i+1:4d} {atom.GetSymbol():2s}         {pos.x:8.3f}        {pos.y:8.3f}         {pos.z:8.3f} HETA 1 0.0\n"
        crd_content1.append(line)

# 写入文件
with open('output.crd', 'w') as f:
    f.writelines(crd_content)