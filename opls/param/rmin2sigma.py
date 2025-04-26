import math

def charmm_to_openmm_lj(epsilon_charmm, rmin_over_2_angstrom):
    """
    将CHARMM的LJ参数转换为OpenMM的sigma和epsilon。
    
    参数:
    epsilon_charmm (float): CHARMM的epsilon参数（单位：kcal/mol）
    rmin_over_2_angstrom (float): CHARMM的Rmin/2参数（单位：Å）
    
    返回:
    (sigma, epsilon): OpenMM的LJ参数，sigma（单位：nm），epsilon（单位：kJ/mol）
    """
    # 转换常数
    kcal_to_kj = 4.184  # 1 kcal/mol = 4.184 kJ/mol
    angstrom_to_nm = 0.1  # 1 Å = 0.1 nm
    
    # 计算Rmin并转换为nm
    rmin_angstrom = rmin_over_2_angstrom * 2
    rmin_nm = rmin_angstrom * angstrom_to_nm
    
    # 计算sigma（2^(1/6) ≈ 1.122462048）
    sigma = rmin_nm / (2 ** (1.0/6.0))
    
    # 计算epsilon（取绝对值并转换单位）
    epsilon = abs(epsilon_charmm) * kcal_to_kj
    
    return sigma, epsilon

# 示例用法

charmm_params = "C245  0.00  -0.066000   1.9643086 0.00  -0.033000   1.9643086"


# 解析CHARMM参数行
parts = charmm_params.split()
atom_type = parts[0]
charge = parts[1]
epsilon_charmm = float(parts[2])
rmin_over_2 = float(parts[3])

# 转换参数
sigma, epsilon = charmm_to_openmm_lj(epsilon_charmm, rmin_over_2)

# 输出OpenMM格式
print(f'<Atom type="{atom_type}" sigma="{sigma:.6f}" epsilon="{epsilon:.6f}" />')