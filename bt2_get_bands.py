#
# Get individual bands from bt2
#
import logging
import os

import matplotlib.pylab as pl
import matplotlib.pyplot as plt
import numpy as np

import contextlib
import copy
import os.path
import sys
import ast
import itertools


from BoltzTraP2 import bandlib as BL
from BoltzTraP2 import dft as BTP
from BoltzTraP2 import fite, serialization, sphere, units
import matplotlib
import ase.dft.kpoints as asekp

from BoltzTraP2.misc import TimerContext, dir_context, info, lexit, warning




def plot_bannds_along_kpath(bt2_file, kpaths, nkpoints):
    band_list, ticks, dividers = get_bands_from_bt2(bt2_file, kpaths, nkpoints)
    for i, kband in enumerate(band_list):
        egrid = kband[1]
        nbands = egrid.shape[0]
        dkp = kband[0]
        for i in range(nbands):
            plt.plot(dkp, egrid[i, :], lw=2.0)
            pass
        pass

    ax = plt.gca()
    ax.set_xticks(ticks)
    ax.set_xticklabels([])
    for d in ticks:
        ax.axvline(x=d, ls="--", lw=0.5)
    for d in dividers:
        ax.axvline(x=d, ls="-", lw=2.0)
    ax.axhline(y=0.0, lw=1.0)
    plt.ylabel(r"$\varepsilon - \varepsilon_F\;\left[\mathrm{Ha}\right]$")
    plt.tight_layout()
    plt.show()


def get_bands_from_bt2(bt2_file, kpaths, nkpoints):
    data, equivalences, coeffs, metadata = serialization.load_calculation(
        bt2_file
        )
    lattvec = data.get_lattvec()
    bands_list = []
    ticks = []
    dividers = []
    offset = 0.0
    for ikpath, kpath in enumerate(kpaths):
        print("k path #{}".format(ikpath + 1))
        # Generate the explicit point list.
        band_path = asekp.bandpath(kpath, data.atoms.cell, nkpoints)
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
            egrid = fite.getBands(
                kp, equivalences, lattvec, coeffs
            )[0]
            deltat = timer.get_deltat()
            print("rebuilding the bands took {:.3g} s".format(deltat))
        egrid -= data.fermi
        bands_list.append([dkp, egrid])
        # Create the plot
        nbands = egrid.shape[0]
        # for i in range(nbands):
        #     plt.plot(dkp, egrid[i, :], lw=2.0)
        ticks += dcl.tolist()
        dividers += [dcl[0], dcl[-1]]
        offset = dkp[-1]
    return bands_list, ticks, dividers




if __name__ == "__main__":

    nkpoints=100
    bt2_file =   "./mp_data/" + "interpolation.bt2"

    # path Gamma-A
    kpaths = [[
        [0.0, 0.0, 0.0],
        [0.0, 0.0, 0.5]
    ]]


    # kpaths = [[
    #     [0.0, 0.0, 0.0],
    #     [0.0, 0.0, 0.5],
    #     [0.5, 0.0, 0.5]
    # ]]


    plot_bannds_along_kpath(bt2_file, kpaths, nkpoints)

