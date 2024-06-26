# SCEMENT : Large scale Single-cell RNA-seq Integration

SCEMENT, a SCalablE and Memory-Efficient iNTegration method for large single-cell RNA-seq datasets.
Our parallel algorithm builds upon and extends the linear regression model applied in ComBat, 
to an unsupervised sparse matrix setting to enable accurate integration of diverse and 
large collections of single cell RNA-sequencing data.

## Dependencies

The software has been only tested on Linux. It requires python 3.11 or above and 
C++ compiler with OpenMP support (tested with gcc 7.5).

## Python Libraries

Depends on the following python libraries

  - psutil
  - numpy
  - scipy
  - patsy
  - anndata
  - scanpy
  - formulaic


## C/C++ Libraries

Eigen and pybind11 libraries are provided as submodules in the ext/ folder.
Other C/C++ libraries needs  

   - CMake 
   - BLAS 
   - LAPACK
   - LAPACKE
   - Armadillo

The above libraries can be installed in Ubuntu/Debian with the following command:

    sudo apt-get install cmake libopenblas-dev liblapack-dev libarpack2-dev libsuperlu-dev armadillo

In Fedora/Red-hat the libraries are:  

    cmake, openblas-devel, lapack-devel, arpack-devel, SuperLU-devel armadillo


## Installation

The package can be installed from source using pip as follows:

    pip install git+https://github.com/AluruLab/scement.git


