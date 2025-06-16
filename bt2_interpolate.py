import BoltzTraP2.dft
from pymatgen.io.vasp import Vasprun
import numpy as np
from BoltzTraP2 import fite
from BoltzTraP2 import serialization

import numpy as np
import scipy as sp
import scipy.linalg as la
import matplotlib
import matplotlib.pyplot as plt

import ase
import spglib

import BoltzTraP2.sphere
import BoltzTraP2.fite
import BoltzTraP2.bandlib
from BoltzTraP2.units import *
import ase.dft.kpoints as asekp
from BoltzTraP2.misc import TimerContext


# Number of k points along each direction in the input grid
NK = 25
# Effective mass of the only band
EFFM = 1.
# Minimum energy of the band
OFFSET = .5
# Amplification factor for the interpolation
FACTOR =  50

HARTREE = 27.2114  # EV

# from pymatgen.io.vasp import Vasprun

filename = "./mp_data/vasprun.xml"
# vasprun = Vasprun(filename, parse_projected_eigen=True)

import BoltzTraP2
from BoltzTraP2.io import parse_vasprunxml

# out = parse_vasprunxml(filename)

# data = BoltzTraP2.dft.VASPLoader("./mp_data/")

def get_dft_from_vasp():
    # Will read vasprun.xml from the folder. 
    # You can use BoltzTraP2.dft.VASPLoader or parse_vasprunxml or Vasprun but I'm using 
    # BoltzTraP2.dft.DFTData
    dft_vasp = BoltzTraP2.dft.DFTData("./mp_data/")
    # Set energy limit here in Hartree unit


    # emin = -0.3/HARTREE
    # emax = 0.3/HARTREE
    # emin = emin + dft_vasp.fermi/HARTREE
    # emax = emax + dft_vasp.fermi/HARTREE
    
    # dft_vasp.bandana(emin, emax)
    # print(data.atoms)
    # print(data.magmom)
    # print("data.kpoints.shape ", data.kpoints.shape)
    lattvec = dft_vasp.get_lattvec()
    # print(lattvec)
    kpoints = dft_vasp.kpoints
    # print("kpoints ", kpoints.shape)
    # print("kpoints ", kpoints)
    magmom = dft_vasp.magmom
    atoms = dft_vasp.atoms
    # print("magmom ", magmom)
    # print("atoms ", atoms)
    # data.kpoints = np.array([
    #     [0.0, 0.0, 0.0],
    #     [0.0, 0.0, 0.5]
    # ])
    return dft_vasp







from pymatgen.electronic_structure.core import Spin
from pymatgen.io.ase import AseAtomsAdaptor
from mp_api.client import MPRester


with open("./apikey", 'r') as f:
    api_key=f.readline()[:-1]
    print(len(api_key))
    pass
material_id="mp-12627"

class DFTdata_Obj_from_mp:
    """
    creates DFTdata object form materials project api
    """
    def __init__(self, api_key, material_id):
        with MPRester(api_key) as mpr:
            bs = mpr.get_bandstructure_by_material_id(material_id)
            structure = mpr.get_structure_by_material_id(material_id)
            self.atoms = AseAtomsAdaptor.get_atoms(structure)
            self.magmom = structure.site_properties['magmom']
            self.lattice_vec = structure.lattice.matrix
            pass
        self.kpoints = np.array(list(map(lambda x: x.frac_coords, bs.kpoints)))
        self.ebands = bs.bands[Spin.up]
        self.mommat = None   # Momentum Matrix
        self.magmom = None   # Magnetic Moment
        self.fermi = bs.efermi
        self.efermi = bs.efermi

        pass

    def get_lattvec(self):
        return self.lattice_vec




def get_seebeck():
    dft = get_dft_from_vasp()
    erange = (-0.02, 0.02) # in Hartree

    # dft = DFTdata_Obj_from_mp(api_key, material_id)
    # erange = (-0.5, 0.5) # in ev

    lattvec = dft.get_lattvec()

    print("dft.fermi ", dft.fermi)
    # emin = emin/HARTREE + dft.fermi
    # emax = emax/HARTREE + dft.fermi
    
    nkpt = int(FACTOR * dft.kpoints.shape[0])
    nkpt = 10000
    print("nkpt =", nkpt)
    equivalences = BoltzTraP2.sphere.get_equivalences(dft.atoms, dft.magmom, nkpt)

    coeffs = BoltzTraP2.fite.fitde3D(dft, equivalences)
    # print("coeffs ", coeffs[0])

    eband, vvband, cband = BoltzTraP2.fite.getBTPbands(equivalences, coeffs, lattvec, curvature=True, nworkers=2)

    dose, dos, vvdos, cdos = BoltzTraP2.bandlib.BTPDOS(eband, vvband, erange=erange, npts=2000)
    
    volume = np.linalg.det(dft.get_lattvec())

    tau = 1e-14
    # Define the temperatures and chemical potentials we are interested in
    Tr = np.arange(4, 300, 20)
    # Tr = np.array([500.])
    margin = 9. * BOLTZMANN * Tr.max()
    mur_indices = np.logical_and(dose > dose.min() + margin,
                                 dose < dose.max() - margin)
    mur = dose[mur_indices]

    mur  = np.array([0])  # Chemical potential
    # Obtain the Fermi integrals required to get the Onsager coefficients
    N, L0, L1, L2, Lm11 = BoltzTraP2.bandlib.fermiintegrals(
        dose, dos, vvdos, mur=mur, Tr=Tr)
    # Translate those into Onsager coefficients
    sigma, seebeck, kappa, Hall = BoltzTraP2.bandlib.calc_Onsager_coefficients(L0, L1, L2, mur, Tr, volume)
    
    print("sigma.shape ", sigma.shape)
    print("seebeck ", seebeck.shape)

    # # Rescale the carrier count into a volumetric density in cm**(-3)
    # N = -N[0, ...] / (volume / (Meter / 100.)**3)
    # # Obtain the scalar conductivity and Seebeck coefficient
    # sigma = tau * sigma[0, ...].trace(axis1=1, axis2=2) / 3.
    # seebeck = seebeck[0, ...].trace(axis1=1, axis2=2) / 3.
    

    # Compute the scalar power factor
    P = sigma * seebeck * seebeck
    # Transform these quantities to more convenient units
    sigma *= 1e-3  # kS / m
    seebeck *= 1e6  # microvolt / K
    P *= 1e6  # microwatt / m / K**2

    # print("seebeck ", seebeck)

    fig,  axes= plt.subplots(2,1,figsize=(4,6), dpi=200)

    i, j = 0, 0
    S_ij = seebeck[:,0,i,j]
    sigma_ij = sigma[:,0,i,j] * tau
    axes[0].plot(Tr, sigma_ij)
    axes[1].plot(Tr, S_ij)
    print("S_{}{}".format(i, j))
    axes[1].set_ylabel("S_{}{}".format(i, j))
    plt.tight_layout()
    plt.savefig("fig/Seebeck-nkpt{}-vasp.png".format(nkpt))
    # plt.show()

    pass

def plot_bands():
    kpaths = [[
        [0.0, 0.0, 0.0],
        [0.0, 0.0, 0.5]
    ]]

    dft = get_dft_from_vasp()
    erange = (-0.02, 0.02) # in Hartree

    # dft = DFTdata_Obj_from_mp(api_key, material_id)
    # erange = (-0.5, 0.5) # in ev

    lattvec = dft.get_lattvec()

    print("dft.fermi ", dft.fermi)
    # emin = emin/HARTREE + dft.fermi
    # emax = emax/HARTREE + dft.fermi
    
    nkpt = int(FACTOR * dft.kpoints.shape[0])
    nkpt = 1000
    print("nkpt =", nkpt)
    equivalences = BoltzTraP2.sphere.get_equivalences(dft.atoms, dft.magmom, nkpt)

    coeffs = BoltzTraP2.fite.fitde3D(dft, equivalences)


    lattvec = dft.get_lattvec()
    bands_list = []
    dkp_list = []
    ticks = []
    dividers = []
    offset = 0.0
    efermi = dft.fermi

    fig, axes = plt.subplots(2,1, figsize=(3,4), dpi=200)
    for ikpath, kpath in enumerate(kpaths):
        print("k path #{}".format(ikpath + 1))
        # Generate the explicit point list.
        band_path = asekp.bandpath(kpath, dft.atoms.cell, nkpt)
        if isinstance(band_path, asekp.BandPath):
            # For newer versions of ASE.
            kp = band_path.kpts
            dkp, dcl = band_path.get_linear_kpoint_axis()[:2]
        else:
            # For older versions of ASE.
            kp, dkp, dcl = band_path
        dkp += offset
        dcl += offset
        # Compute the band energies
        with TimerContext() as timer:
            eband, vband, cband = fite.getBands(
                kp, equivalences, lattvec, coeffs, curvature=True
            )
            deltat = timer.get_deltat()
            print("rebuilding the bands took {:.3g} s".format(deltat))
            pass
        

        axes[0].plot(dkp, eband.T-efermi, lw=2.0)
        
        # break
        # vband_scalar = list(map(np.linalg.norm, vband[:, i, :].T))
        # print(dkp.shape)
        # print(type(vband_scalar))
        # print(len(vband_scalar))
        # axes[1].plot(dkp, vband_scalar, lw=2.0)
        
        dkp_list.append(dkp)
        # Create the plot
        # nbands = egrid.shape[0]
        # for i in range(nbands):
        #     plt.plot(dkp, egrid[i, :], lw=2.0)
        ticks += dcl.tolist()
        dividers += [dcl[0], dcl[-1]]
        offset = dkp[-1]
        pass

    
    
    axes[0].set_xticks(ticks)
    axes[0].set_xticklabels([])
    for d in ticks:
        axes[0].axvline(x=d, ls="--", lw=0.5, color='k')
    for d in dividers:
        axes[0].axvline(x=d, ls="-", lw=2.0, color='k')
    axes[0].axhline(y=0.0, lw=1.0, color='k')
    axes[0].set_ylabel(r"$\varepsilon - \varepsilon_F\;\left[\mathrm{Ha}\right]$")
    axes[0].set_ylim(erange)
    axes[1].set_ylabel(r"$v$")
    plt.tight_layout()
    plt.show()


        
    pass


if __name__ == "__main__":  
    # get_seebeck()
    plot_bands()



