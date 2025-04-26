from parmed.charmm import CharmmParameterSet
from simtk.openmm.app import ForceField

# 1. 读取str文件
def read_params(filename):
    extlist = ["rtf", "prm", "str"]

    parFiles = ()
    for line in open(filename, "r"):
        if "!" in line:
            line = line.split("!")[0]
        parfile = line.strip()
        if len(parfile) != 0:
            ext = parfile.lower().split(".")[-1]
            if not ext in extlist:
                continue
            parFiles += (parfile,)

    params = CharmmParameterSet(*parFiles)
    return params

# 2. 读取str文件并生成参数集
params = read_params('/pubhome/xtzhang/interaction/FF/charmm/ChgParam/ACEM.str')

# 3. 将参数集保存为OpenMM的XML格式
params.write(str='/pubhome/xtzhang/interaction/FF/charmm/ChgParam/ACEM.xml')
