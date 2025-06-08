
def run_antechamber(file, prefix, netcharge=0):
    subprocess.run([
        "/pubhome/xtzhang/mambaforge/envs/htmd/bin/antechamber", "-i", file, "-fi", "sdf",
        "-o", f"{prefix}.mol2", "-fo", "mol2", "-c", "bcc", "-s", "2", "-nc", f"{netcharge} -at gaff2"
    ])  # generate mol2 with bcc charge
    subprocess.run([
        "/pubhome/xtzhang/mambaforge/envs/htmd/bin/parmchk2", "-i", f"{prefix}.mol2", "-f", "mol2", "-o", f"{prefix}.frcmod -s gaff2"
    ]) # generate mainly dihedral param

def write_leap(prefix):
    with open("leap.in", "w") as f:
        f.write("source leaprc.gaff2\n")
        f.write(f"LIG = loadmol2 {prefix}.mol2\n")
        f.write(f"loadamberparams {prefix}.frcmod\n")
        f.write(f"saveamberparm LIG {prefix}.prmtop {prefix}.inpcrd\n")# generate topology and coordinate for specific conformation 
        f.write("quit\n")
        f.write("quit\n")
    subprocess.run(["/pubhome/xtzhang/mambaforge/envs/htmd/bin/tleap", "-f", "leap.in"])
