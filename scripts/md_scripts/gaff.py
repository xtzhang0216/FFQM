#%%
from parmed import load_file, Atom, Structure
from parmed.amber import AmberParm
from parmed.tools import addLJType
from simtk.openmm.app import ForceField
from simtk.openmm import XmlSerializer

# 1. 加载mol2文件
# mol = load_file(mol2_file)

# 2. 使用GAFF力场分配参数
# 需要安装openff-toolkit并调用antechamber
from openff.toolkit import Molecule
from openmmforcefields.generators import (
    GAFFTemplateGenerator,
)
from openff.toolkit.typing.engines.smirnoff import ForceField as OFFForceField
#%%
# 1. 创建小分子对象（示例为苯）
molecule = Molecule.from_smiles("c1ccccc1")
# 2. 生成GAFF力场模板生成器
gaff = GAFFTemplateGenerator(molecules=molecule)
# 3. 创建OpenMM力场对象，并注册GAFF生成器
forcefield = ForceField()
forcefield.registerTemplateGenerator(gaff.generator)
# 4. 生成分子的拓扑和参数
#   注意：需要分子坐标！先生成三维结构或加载PDB文件
molecule.generate_conformers(n_conformers=1)
topology = molecule.to_topology().to_openmm()
positions = molecule.conformers[0].to_openmm()
# 5. 创建OpenMM系统（触发参数生成）
system = forcefield.createSystem(topology, molecules=molecule)
# 6. 将系统序列化为XML文件
with open("/pubhome/xtzhang/interaction/FF/charmm/para", "w") as f:
    f.write(XmlSerializer.serialize(system))

# %%
