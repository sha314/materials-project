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


class DFTdata_MP:
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
        print("subtracting efermi in __init__()")
        self.efermi = bs.efermi
        self.ebands = bs.bands[Spin.up] - self.efermi
        self.mommat = None   # Momentum Matrix
        self.magmom = None   # Magnetic Moment
        self.fermi = bs.efermi
        self.source = "Custome DFT from materials project api"
        self.erange = None

        pass

    def bands_by_erange(self, emin, emax):
        """
        emin, emax: energy in eV relative to fermi level
        """
        # emin += self.efermi
        # emax += self.efermi
        self.erange = (emin, emax)
        selected_indices = [i for i, band in enumerate(self.ebands) if any(emin <= e <= emax for e in band)]
        self.bands_by_index(selected_indices)
        pass

    def bands_by_index(self, selected_indices):
        print("selected bands ", selected_indices)
        # print("total band count ", len(selected_indices))
        self.ebands = self.ebands[selected_indices,:]
        # print("self.ebands.shape ", self.ebands.shape)
        pass

    def plot_bands(self):
        x = np.linspace(0, 1, self.ebands.shape[1])
        for i in range(self.ebands.shape[0]):
            plt.plot(x, self.ebands[i])
            pass
        plt.show()

    def get_lattvec(self):
        return self.lattice_vec
    




def get_seebeck(material_id):
    # dft = get_dft_from_vasp()
    # erange = (-0.02, 0.02) # in Hartree

    dft = DFTdata_MP(api_key, material_id)
    erange = (dft.efermi - 3 , dft.efermi + 3) # in ev

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

    print("does ", dose.shape)
    print("dos ", dos.shape)

    # Extract energy grid and DOS values
    energies = dose  # shape: (npts,)
    # total_DOS = dos.sum(axis=0)  # sum over all bands
    plt.plot(energies, dos, label="Total DOS")
    plt.axvline(0, color='k', linestyle='--', label='Fermi level')  # Fermi level at 0 eV
    plt.xlabel("Energy (eV)")
    plt.ylabel("DOS (states/eV)")
    plt.title("Density of States from BoltzTraP2")
    plt.legend()
    plt.tight_layout()
    plt.show()
        

    
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
    plt.savefig("fig/Seebeck-xx-nkpt{}-mp.png".format(nkpt))
    # plt.show()

    pass

def get_seebeck_data(material_id):
    # dft = get_dft_from_vasp()
    # erange = (-0.02, 0.02) # in Hartree

    dft = DFTdata_MP(api_key, material_id)
    dft.bands_by_erange(-.3,.3)
    # dft.plot_bands()

    # erange = (dft.efermi - 3 , dft.efermi + 3) # in ev
    erange = None

    lattvec = dft.get_lattvec()

    print("dft.fermi ", dft.fermi)
    # emin = emin/HARTREE + dft.fermi
    # emax = emax/HARTREE + dft.fermi
    
    nkpt = int(FACTOR * dft.kpoints.shape[0])
    nkpt = 100_000
    print("nkpt =", nkpt)
    equivalences = BoltzTraP2.sphere.get_equivalences(dft.atoms, dft.magmom, nkpt)

    coeffs = BoltzTraP2.fite.fitde3D(dft, equivalences)
    # print("coeffs ", coeffs[0])

    eband, vvband, cband = BoltzTraP2.fite.getBTPbands(equivalences, coeffs, lattvec, curvature=True, nworkers=2)

    dataToFile = dict()
    for binNum in [2000, 5000, 10_000, 20_000]:
        print("binNum ", binNum)
        data = dict()
        dose, dos, vvdos, cdos = BoltzTraP2.bandlib.BTPDOS(eband, vvband, erange=erange, npts=binNum)
        
        volume = np.linalg.det(dft.get_lattvec())

        tau = 1e-14
        # Define the temperatures and chemical potentials we are interested in
        Tr = np.array([5, 10, 15, 20, 30, 40, 55, 80, 110, 130, 160, 190, 220, 260, 300, 350, 400, 450, 500])
        data['temp'] = Tr
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

        i, j = 0, 0
        S_ij = seebeck[:,0,i,j]
        print("S max ", np.max(np.abs(seebeck[:,0,0,0])))
        data['seebeck_xx'] = seebeck[:,0,0,0]
        data['seebeck_yy'] = seebeck[:,0,1,1]
        data['seebeck_zz'] = seebeck[:,0,2,2]
        data['nkpt'] = nkpt
        data['npts'] = binNum
        data['dos_x'] = dos[0]
        data['dos_y'] = dos[1]
        data['erange'] = dft.erange

        dataToFile[binNum] = data
        
        # plt.show()
        pass
    import json
    import pickle
    # Write to JSON file
    with open("seebeck-mp-nkpt{}.pkl".format(nkpt), "wb") as f:
        pickle.dump(dataToFile, f)  # indent for pretty-printing


    pass

def plot_bands(material_id):
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
    nkpt = 2000
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

    fig, axes = plt.subplots(1,1, figsize=(5,3), dpi=200)
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
        

        plt.plot(dkp, eband.T-efermi, lw=2.0)
        
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

    
    
    axes.set_xticks(ticks)
    axes.set_xticklabels([])
    for d in ticks:
        plt.axvline(x=d, ls="--", lw=0.5, color='k')
    for d in dividers:
        plt.axvline(x=d, ls="-", lw=2.0, color='k')
    plt.axhline(y=0.0, lw=1.0, color='k')
    axes.set_ylabel(r"$\varepsilon - \varepsilon_F\;\left[\mathrm{Ha}\right]$")
    axes.set_ylim(erange)
    # axes[1].set_ylabel(r"$v$")
    plt.tight_layout()
    
    plt.savefig("fig/bands-vasp-nkpt{}.png".format(nkpt))
    # plt.show()

        
    pass

def plot_bands_and_dos(material_id):
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
    nkpt = 2000
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

    fig, axes = plt.subplots(1,1, figsize=(5,3), dpi=200)
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
        

        plt.plot(dkp, eband.T-efermi, lw=2.0)
        
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

    
    
    axes.set_xticks(ticks)
    axes.set_xticklabels([])
    for d in ticks:
        plt.axvline(x=d, ls="--", lw=0.5, color='k')
    for d in dividers:
        plt.axvline(x=d, ls="-", lw=2.0, color='k')
    plt.axhline(y=0.0, lw=1.0, color='k')
    axes.set_ylabel(r"$\varepsilon - \varepsilon_F\;\left[\mathrm{Ha}\right]$")
    axes.set_ylim(erange)
    # axes[1].set_ylabel(r"$v$")
    plt.tight_layout()
    
    plt.savefig("fig/bands-vasp-nkpt{}.png".format(nkpt))
    # plt.show()

        
    pass


def get_interpolation(material_id, out_file="interpolation.bt2"):
    # dft = get_dft_from_vasp()
    # erange = (-0.02, 0.02) # in Hartree

    dft = DFTdata_MP(api_key, material_id)
    erange = (-0.5, 0.5) # in ev

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

    metadata = BoltzTraP2.serialization.gen_bt2_metadata(
        dft, True
    )

    BoltzTraP2.serialization.save_calculation(
            out_file, dft, equivalences, coeffs, metadata
        )
    
    

    pass

if __name__ == "__main__":  
    material_id="mp-12627"  # Nb3S4
    # material_id="mp-924130"  # TiNiSn

    # get_seebeck(material_id)
    get_seebeck_data(material_id)
    # plot_bands(material_id)
    # get_interpolation(material_id)



