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

$ conda create -n band

$ conda activate band

$ conda install conda-forge::mp-api

$ conda install jupyter

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
 

