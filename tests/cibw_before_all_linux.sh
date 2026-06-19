#!/bin/bash
# cibuildwheel before-all step (runs once inside the manylinux container).
# Stage the MPI headers where the CMake build expects them (see the NETGEN_MPI
# override in pyproject.toml). Netgen is built with the MPI wrapper, so only the
# headers are needed at build time -- the actual MPI library is dispatched at
# run time.
set -e

mkdir -p /opt/openmpi/include /opt/mpich/include
cp -a /usr/include/openmpi-x86_64/* /opt/openmpi/include/
cp -a /usr/include/mpich-x86_64/* /opt/mpich/include/
