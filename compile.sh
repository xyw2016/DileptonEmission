#!/usr/bin/env bash

module load nixpkgs/16.09  intel/2018.3  impi/2018.3.222
module load gsl
module load hdf5-mpi/1.8.18
module load StdEnv/2020 intel/2020.1.217 hdf5/1.12.1
module load gcc

module load cmake

num_of_cores=$1
if [ -z "$num_of_cores" ]
then
    num_of_cores=4
fi

mkdir -p build
cd build
rm -fr *
cmake ..
make -j$num_of_cores
make install
cd ..
rm -fr build
