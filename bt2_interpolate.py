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
FACTOR = 5

HARTREE = 27.2114  # EV

# from pymatgen.io.vasp import Vasprun

filename = "./mp_data/vasprun.xml"
# vasprun = Vasprun(filename, parse_projected_eigen=True)

import BoltzTraP2
from BoltzTraP2.io import parse_vasprunxml

# out = parse_vasprunxml(filename)

# data = BoltzTraP2.dft.VASPLoader("./mp_data/")

# Will read vasprun.xml from the folder. 
# You can use BoltzTraP2.dft.VASPLoader or parse_vasprunxml or Vasprun but I'm using 
# BoltzTraP2.dft.DFTData
data = BoltzTraP2.dft.DFTData("./mp_data/")
# Set energy limit here in Hartree unit
e_lim = 0.3/HARTREE
# print(e_lim)
data.bandana(-.5, .5)
# print(data.atoms)
# print(data.magmom)
# print("data.kpoints.shape ", data.kpoints.shape)
lattvec = data.get_lattvec()
kpoints = data.kpoints
magmom = data.magmom
atoms = data.atoms

if __name__ == "__main__":  

    equivalences = BoltzTraP2.sphere.get_equivalences(atoms, magmom, FACTOR * kpoints.shape[0])
    # print(equivalences)

    coeffs = BoltzTraP2.fite.fitde3D(data, equivalences)

    
    serialization.save_calculation("./mp_data/interpolation.bt2", data, equivalences, coeffs,
                                   serialization.gen_bt2_metadata(
                                       data, data.mommat is not None))


    # print(coeffs.shape)

    eband, vvband, cband = BoltzTraP2.fite.getBTPbands(equivalences, coeffs, lattvec, curvature=True, nworkers=2)

    dose, dos, vvdos, cdos = BoltzTraP2.bandlib.BTPDOS(
        eband, vvband, erange=(OFFSET - .1, np.max(eband)), npts=2000)

    volume = np.linalg.det(data.get_lattvec())

    tau = 1e-14
    # Define the temperatures and chemical potentials we are interested in
    Tr = np.array([500.])
    margin = 9. * BOLTZMANN * Tr.max()
    mur_indices = np.logical_and(dose > dose.min() + margin,
                                 dose < dose.max() - margin)
    mur = dose[mur_indices]

    # Obtain the Fermi integrals required to get the Onsager coefficients
    N, L0, L1, L2, Lm11 = BoltzTraP2.bandlib.fermiintegrals(
        dose, dos, vvdos, mur=mur, Tr=Tr)
    # Translate those into Onsager coefficients
    sigma, seebeck, kappa, Hall = BoltzTraP2.bandlib.calc_Onsager_coefficients(
        L0, L1, L2, mur, Tr, volume)
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
    pass



