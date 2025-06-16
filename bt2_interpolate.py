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

# Number of k points along each direction in the input grid
NK = 25
# Effective mass of the only band
EFFM = 1.
# Minimum energy of the band
OFFSET = .5
# Amplification factor for the interpolation
FACTOR =  5

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


    emin = -0.3/HARTREE
    emax = 0.3/HARTREE
    emin = emin + dft_vasp.fermi/HARTREE
    emax = emax + dft_vasp.fermi/HARTREE
    dft_vasp.bandana(emin, emax)
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







if __name__ == "__main__":  
    dft = get_dft_from_vasp()
    # dft = DFTdata_Obj_from_mp(api_key, material_id)
    
    lattvec = dft.get_lattvec()

    emin, emax = -0.1, 0.1
    emin = emin + dft.fermi/HARTREE
    emax = emax + dft.fermi/HARTREE


    nkpt = int(FACTOR * dft.kpoints.shape[0])
    print("nkpt =", nkpt)
    equivalences = BoltzTraP2.sphere.get_equivalences(dft.atoms, dft.magmom, nkpt)

    coeffs = BoltzTraP2.fite.fitde3D(dft, equivalences)

    eband, vvband, cband = BoltzTraP2.fite.getBTPbands(equivalences, coeffs, lattvec, curvature=True, nworkers=2)

    dose, dos, vvdos, cdos = BoltzTraP2.bandlib.BTPDOS(eband, vvband, erange=(emin, emax), npts=2000)
    
    volume = np.linalg.det(dft.get_lattvec())

    tau = 1e-14
    # Define the temperatures and chemical potentials we are interested in
    Tr = np.array([10., 20., 40., 80., 100., 150., 200., 300.])
    Tr = np.array([500.])
    margin = 9. * BOLTZMANN * Tr.max()
    mur_indices = np.logical_and(dose > dose.min() + margin,
                                 dose < dose.max() - margin)
    mur = dose[mur_indices]

    # Obtain the Fermi integrals required to get the Onsager coefficients
    N, L0, L1, L2, Lm11 = BoltzTraP2.bandlib.fermiintegrals(
        dose, dos, vvdos, mur=mur, Tr=Tr)
    # Translate those into Onsager coefficients
    sigma, seebeck, kappa, Hall = BoltzTraP2.bandlib.calc_Onsager_coefficients(L0, L1, L2, mur, Tr, volume)
    print("sigma.shape ", sigma.shape)
    # Rescale the carrier count into a volumetric density in cm**(-3)
    N = -N[0, ...] / (volume / (Meter / 100.)**3)
    # Obtain the scalar conductivity and Seebeck coefficient
    sigma = tau * sigma[0, ...].trace(axis1=1, axis2=2) / 3.
    seebeck = seebeck[0, ...].trace(axis1=1, axis2=2) / 3.
    # Compute the scalar power factor
    P = sigma * seebeck * seebeck
    # Transform these quantities to more convenient units
    sigma *= 1e-3  # kS / m
    seebeck *= 1e6  # microvolt / K
    P *= 1e6  # microwatt / m / K**2

    fig,  axes= plt.subplots(2,1,figsize=(3,4), dpi=200)
    print(sigma.shape)
    print(seebeck.shape)
    axes[0].plot(Tr, sigma)
    axes[1].plot(Tr, seebeck)
    plt.show()
    pass



