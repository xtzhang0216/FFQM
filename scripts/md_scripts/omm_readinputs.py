from simtk.unit import *
from simtk.openmm import *
from simtk.openmm.app import *


class _OpenMMReadInputs:
    def __init__(self):

        self.pdbid = None  # PDBID of protein
        self.toppar = "toppar_drude_openmm.str"  # Path to toppar.str file
        self.psffile = "3kwj.psf"  # Path to .psf file
        self.pdbfile = "3kwj.pdb"  # Path to .pdb file

        self.mini_nstep = 5000  # Number of steps for minimization
        self.mini_Tol = 10.0  # Minimization energy tolerance
        self.gen_vel = "yes"  # Generate initial velocities
        self.gen_temp = 298.15  # Temperature for generating initial velocities (K)
        self.gen_seed = None  # Seed for generating initial velocities

        self.nvt_nstep = 50000  # Number of steps for NVT equilibrium
        self.nvt_nstout = 100  # Writing output frequency (steps)
        self.nvt_nstdcd = 50000  # Writing coordinates trajectory frequency (steps)
        self.nvt_dt = 0.0001  # Time-step of NVT (ps)

        self.npt_nstep = 50000  # Number of steps for NPT equilibrium
        self.npt_nstout = 100  # Writing output frequency (steps)
        self.npt_nstdcd = 50000  # Writing coordinates trajectory frequency
        self.npt_dt = 0.0005  # Time-step of N

        self.reduce_nstep = 500000  # Number of steps for reduce restraint NPT
        self.reduce_nstout = 100  # Writing output frequency (steps)
        self.reduce_nstdcd = 50000  # Writing coordinates trajectory frequency
        self.reduce_dt = 0.002  # Time-step of reduce restraint
        self.turn = 5

        self.start = 0  # Start number of simulation
        self.end = 10  # End number of simulations
        self.nstep = 50000  # Number of steps to run in every cycle
        self.dt = 0.001  # Time-step (ps)
        self.nstout = 1000  # Writing output frequency (steps)
        self.nstdcd = 1000  # Wrtiing coordinates trajectory frequency (steps)

        self.coulomb = PME  # Electrostatic cut-off method
        self.ewald_Tol = 0.0001  # Ewald error tolerance
        self.vdw = "Switch"  # vdW cut-off method
        self.r_on = 1.0  # Switch-on distance (nm)
        self.r_off = 1.2  # Switch-off distance (nm)

        self.temp = 298.15  # Temperature (K)
        self.fric_coeff = 5  # Friction coefficient for Langevin dynamics

        self.drude_temp = 1.0  # Drude Temperature (K)
        self.drude_fric_coeff = 20  # Drude Friction coefficient for Langevin dynamics
        self.drude_hardwall = 0.025  # Drude Hardwall

        self.pcouple = "yes"  # Turn on/off pressure coupling
        self.p_ref = 1.0  # Pressure (Pref or Pxx, Pyy, Pzz; bar)
        self.p_xx = 1.0  # Pressure Pxx (bar)
        self.p_yy = 1.0  # Pressure Pyy (bar)
        self.p_zz = 1.0  # Pressure Pzz (bar)
        self.p_type = "anisotropic"  # MonteCarloBarotat type
        self.p_scalex = (
            True  # whether to allow the X dimension of the periodic box to change size
        )
        self.p_scaley = (
            True  # whether to allow the Y dimension of the periodic box to change size
        )
        self.p_scalez = (
            True  # whether to allow the Z dimension of the periodic box to change size
        )
        self.p_XYMode = (
            MonteCarloMembraneBarostat.XYIsotropic
        )  # For MonteCarloMembraneBarostat
        self.p_ZMode = (
            MonteCarloMembraneBarostat.ZFree
        )  # For MonteCarloMembraneBarostat
        self.p_tens = 0.0  # Sulface tension for MonteCarloMembraneBarostat (dyne/cm)
        self.p_freq = 15  # Pressure coupling frequency (steps)

        self.cons = HBonds  # Constraints method

        self.rest = "no"  # Turn on/off restraints
        self.rest_type = "CA"  # restrint type
        self.fc_pos = 50.208  # Positional restraint force constant
        self.ftype = "drude"  # Type of force field

        self.rest_lig = "no"  # Turn on/off range restraints on ligand
        self.fc_pos_lig = 1000.0  # Positional restraint force constant on ligand when the ligand run out of pocket
        self.pocket = 5  # The radius of pocket(A)

    def read(self, inputFile):
        for line in open(inputFile, "r"):
            if line.find("#") >= 0:
                line = line.split("#")[0]
            line = line.strip()
            if len(line) > 0:
                segments = line.split("=")
                input_param = segments[0].upper().strip()
                try:
                    input_value = segments[1].strip()
                except:
                    input_value = None
                if input_value:
                    if input_param == "START":
                        self.start = int(input_value)
                    if input_param == "END":
                        self.end = int(input_value)
                    if input_param == "PDBID":
                        self.pdbid = str(input_value)
                    if input_param == "TOPPAR":
                        self.toppar = str(input_value)
                    if input_param == "PSFFILE":
                        self.psffile = str(input_value)
                    if input_param == "PDBFILE":
                        self.pdbfile = str(input_value)
                    if input_param == "MINI_NSTEP":
                        self.mini_nstep = int(input_value)
                    if input_param == "MINI_TOL":
                        self.mini_Tol = float(input_value)
                    if input_param == "GEN_VEL":
                        if input_value.upper() == "YES":
                            self.gen_vel = "yes"
                        if input_value.upper() == "NO":
                            self.gen_vel = "no"
                    if input_param == "GEN_TEMP":
                        self.gen_temp = float(input_value)
                    if input_param == "GEN_SEED":
                        if input_value.upper() == "NONE":
                            self.gen_seed = None
                        else:
                            self.gen_seed = int(input_value)
                    if input_param == "NVT_NSTEP":
                        self.nvt_nstep = int(input_value)
                    if input_param == "NVT_NSTOUT":
                        self.nvt_nstout = int(input_value)
                    if input_param == "NVT_NSTDCD":
                        self.nvt_nstdcd = int(input_value)
                    if input_param == "NVT_DT":
                        self.nvt_dt = float(input_value)
                    if input_param == "NPT_NSTEP":
                        self.npt_nstep = int(input_value)
                    if input_param == "NPT_NSTOUT":
                        self.npt_nstout = int(input_value)
                    if input_param == "NPT_NSTDCD":
                        self.npt_nstdcd = int(input_value)
                    if input_param == "NPT_DT":
                        self.npt_dt = float(input_value)
                    if input_param == "REDUCE_NSTEP":
                        self.reduce_nstep = int(input_value)
                    if input_param == "REDUCE_NSTOUT":
                        self.npt_nstout = int(input_value)
                    if input_param == "REDUCE_NSTDCD":
                        self.npt_nstdcd = int(input_value)
                    if input_param == "REDUCE_DT":
                        self.npt_dt = float(input_value)
                    if input_param == "TURN":
                        self.turn = int(input_value)
                    if input_param == "NSTEP":
                        self.nstep = int(input_value)
                    if input_param == "DT":
                        self.dt = float(input_value)
                    if input_param == "NSTOUT":
                        self.nstout = int(input_value)
                    if input_param == "NSTDCD":
                        self.nstdcd = int(input_value)
                    if input_param == "COULOMB":
                        if input_value.upper() == "NOCUTOFF":
                            self.coulomb = NoCutoff
                        if input_value.upper() == "CUTOFFNONPERIODIC":
                            self.coulomb = CutoffNonPeriodic
                        if input_value.upper() == "CUTOFFPERIODIC":
                            self.coulomb = CutoffPeriodic
                        if input_value.upper() == "EWALD":
                            self.coulomb = Ewald
                        if input_value.upper() == "PME":
                            self.coulomb = PME
                    if input_param == "EWALD_TOL":
                        self.ewald_Tol = float(input_value)
                    if input_param == "VDW":
                        if input_value.upper() == "FORCE-SWITCH":
                            self.vdw = "Force-switch"
                        if input_value.upper() == "SWITCH":
                            self.vdw = "Switch"
                    if input_param == "R_ON":
                        self.r_on = float(input_value)
                    if input_param == "R_OFF":
                        self.r_off = float(input_value)
                    if input_param == "TEMP":
                        self.temp = float(input_value)
                    if input_param == "FRIC_COEFF":
                        self.fric_coeff = float(input_value)
                    if input_param == "DRUDE_TEMP":
                        self.drude_temp = float(input_value)
                    if input_param == "DRUDE_FRIC_COEFF":
                        self.drude_fric_coeff = float(input_value)
                    if input_param == "DRUDE_HARDWALL":
                        self.drude_hardwall = float(input_value)
                    if input_param == "PCOUPLE":
                        if input_value.upper() == "YES":
                            self.pcouple = "yes"
                        if input_value.upper() == "NO":
                            self.pcouple = "no"
                    if input_param == "P_REF":
                        if input_value.find(",") < 0:
                            self.p_ref = float(input_value)
                        else:
                            Pxx = float(input_value.split(",")[0])
                            Pyy = float(input_value.split(",")[1])
                            Pzz = float(input_value.split(",")[2])
                            self.p_ref = Pxx, Pyy, Pzz
                    if input_param == "P_TYPE":
                        if input_value.upper() == "ISOTROPIC":
                            self.p_type = "isotropic"
                        if input_value.upper() == "MEMBRANE":
                            self.p_type = "membrane"
                        if input_value.upper() == "ANISOTROPIC":
                            self.p_type = "anisotropic"
                    if input_param == "P_SCALEX":
                        self.p_scalex = input_value
                    if input_param == "P_SCALEY":
                        self.p_scaley = input_value
                    if input_param == "P_SCALEZ":
                        self.p_scalez = input_value
                    if input_param == "P_XYMODE":
                        if input_value.upper() == "XYISOTROPIC":
                            self.p_XYMode = MonteCarloMembraneBarostat.XYIsotropic
                        if input_value.upper() == "XYANISOTROPIC":
                            self.p_XYMode = MonteCarloMembraneBarostat.XYAnisotropic
                    if input_param == "P_ZMODE":
                        if input_value.upper() == "ZFREE":
                            self.p_ZMode = MonteCarloMembraneBarostat.ZFree
                        if input_value.upper() == "ZFIXED":
                            self.p_ZMode = MonteCarloMembraneBarostat.ZFixed
                        if input_value.upper() == "CONSTANTVOLUME":
                            self.p_ZMode = MonteCarloMembraneBarostat.ConstantVolume
                    if input_param == "P_TENS":
                        self.p_tens = float(input_value)
                    if input_param == "P_FREQ":
                        self.p_freq = int(input_value)
                    if input_param == "CONS":
                        if input_value.upper() == "NONE":
                            self.cons = None
                        if input_value.upper() == "HBONDS":
                            self.cons = HBonds
                        if input_value.upper() == "ALLBONDS":
                            self.cons = AllBonds
                        if input_value.upper() == "HANGLES":
                            self.cons = HAngles
                    if input_param == "REST":
                        if input_value.upper() == "YES":
                            self.rest = "yes"
                        if input_value.upper() == "NO":
                            self.rest = "no"
                    if input_param == "REST_TYPE":
                        self.rest_type = str(input_value)
                    if input_param == "FC_POS":
                        self.fc_pos = float(input_value)
                    if input_param == "FTYPE":
                        self.ftype = str(input_value)
                    if input_param == "REST_LIG":
                        if input_value.upper() == "YES":
                            self.rest_lig = "yes"
                        if input_value.upper() == "NO":
                            self.rest_lig = "no"
                    if input_param == "FC_POS_LIG":
                        self.fc_pos_lig = float(input_value)
                    if input_param == "POCKET":
                        self.pocket = float(input_value)
        return self


def read_inputs(inputFile):
    return _OpenMMReadInputs().read(inputFile)
