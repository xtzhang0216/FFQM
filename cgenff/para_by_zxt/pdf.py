#%%
from parmed import load_file, charmm
import os

os.chdir("/pubhome/xtzhang/interaction/FF/charmm/para")
# 1. 加载坐标文件（如PDB）
struct = load_file('/pubhome/xtzhang/interaction/FF/charmm/para/n1pa.pdb')

# 2. 加载CHARMM力场拓扑和参数
#    需要指定拓扑文件、参数文件、流文件（str）
charmm_top = charmm.CharmmParameterSet(
    'toppar/top_all36_prot.rtf',     # 拓扑文件
    'toppar/par_all36m_prot.prm',     # 参数文件
    'N1PA.str'         # 自定义流文件（可选）
)

# 3. 生成PSF和坐标
#    假设结构已包含正确的残基和连接性
psf = charmm.CharmmPsfFile.from_structure(struct)

# 4. 保存PSF和修正后的坐标
psf.save('output.psf')
struct.save('output.crd')

# %%
