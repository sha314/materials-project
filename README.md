# materials-project
scripts for using materials-project API

# Steps
1. Create free account at materialsproject.org
2. Install the api (It's done using conda here)
3.  

# Official Documentation
https://docs.materialsproject.org/downloading-data/using-the-api/getting-started


# Installing using conda
1. Install conda. You can install anaconda or miniconda from https://anaconda.org
2. I'm gonna create a new environment in conda and install everything there
```
conda create -n band python=3.8
conda activate band
conda install jupyter mp-api -c conda-forge
conda install pymatgen=2024.3 -c conda-forge
```

To see the version information
```
pip show pymatgen
```


3. Run jupyter notebook file

4. If you have your api saved in a file named "apikey" then the program will read api-key from file and use it.



# Indexing of bandstructure

branch   : A particular direction in k-space

branches : Contains information about 'start' and 'end' index of all branches. 

kpoints  : List of 3D vectors in k-space. Index range specified in branches for a particular branch corresponds to particular direction or band.

bands    : an array with indexing [band_index, kpoint_index]

projections (dict[Spin, NDArray]): Orbital projections as {spin: array}.

 |              The indices of the array are [band_index, kpoint_index, orbital_index, ion_index].
 |               If the band structure is not spin polarized, we only store one data set under Spin.up.
 

# Installing BoltzTraP2
It solves the Boltzmann transport equation for electrons under the constant relaxation time approximation (CRTA), using electronic band structures (E vs. k) and optionally the density of states (DOS) and crystal symmetry.

https://gitlab.com/sousaw/BoltzTraP2.git


```
conda activate band
conda install numpy scipy matplotlib cython spglib ase gfortran cmake make -c conda-forge
conda install -c conda-forge phonopy
```

Make sure f2py is working:
```
f2py -v
```

```
conda install pyfftw jupyter vtk -c conda-forge
```

or 
```
pip install vtk
```


```
git clone https://gitlab.com/sousaw/BoltzTraP2.git
cd BoltzTraP2
export CMAKE_POLICY_VERSION_MINIMUM=3.5
pip install .
```

Or, try this if pip doesn't work
```
python setup.py install
```




https://gitlab.com/yiwang62/BoltzTraP2/-/blob/20210126/BoltzTrap2Y.ipynb



Terminal commands
```
btp2 interpolate -m 5 .
btp2 interpolate -e -0.1 -E 0.1 -m 5 .
btp2 dope -b 1000 interpolation.bt2 30:300:10 0
btp2 -vv integrate interpolation.bt2 30:300:10 0
```

Fermisurface
```
btp2 -n 10 fermisurface interpolation.bt2 0
```


# For regression tasks
```
conda install scikit-learn
```


```
environment
```

# BoltzTrap2 as a module



# Export and Install from .yml file
To export
```
conda env export > environment.yml
```

To install
```
conda env create -f environment.yml
```
