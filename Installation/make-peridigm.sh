#!/bin/bash

# Define paths to dependencies
Trilinos=/cluster/stages/2024.0/software/trilinos/install-peridigm-13-0
OPENMPI=/cluster/stages/2024.0/spack-0.22/opt/spack/linux-rocky8-cascadelake/gcc-11.4.0/openmpi-4.1.6-5voubczsrbixdifmpczll66tpm34th25
HDF5=/cluster/stages/2024.0/spack-0.22/opt/spack/linux-rocky8-cascadelake/gcc-11.4.0/hdf5-1.12.0-sxmurkczktjr357tjvp3bdoq5j5bhprs

# Set your home directory paths
PERIDIGM_SOURCE=~/peridigm/source
INSTALL_DIR=~/peridigm/install
BUILD_DIR=~/peridigm/build

# Remove any previous build cache
rm -f ${BUILD_DIR}/CMakeCache.txt

# Configure the Peridigm build
cmake \
  -D CMAKE_INSTALL_PREFIX:PATH=${INSTALL_DIR} \
  -D CMAKE_BUILD_TYPE:STRING=Release \
  -D Trilinos_DIR:PATH=${Trilinos}/lib/cmake/Trilinos \
  -D CMAKE_C_COMPILER:STRING=${OPENMPI}/bin/mpicc \
  -D CMAKE_CXX_COMPILER:STRING=${OPENMPI}/bin/mpicxx \
  -D CMAKE_Fortran_COMPILER:FILEPATH=$(which mpif90) \
  -D BOOST_ROOT=${BOOST_ROOT}/ \
  -D HDF5_DIR=${HDF5} \
  -D CMAKE_CXX_FLAGS:STRING="-O2 -Wall -pedantic -Wno-long-long -ftrapv -Wno-deprecated" \
  ${PERIDIGM_SOURCE}

make -j4
make install