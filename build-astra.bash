#!/bin/bash
# ARM_64 processors with gcc and clang compilers

# re-initialize modules to purge stuff from the parent shell
/etc/profile.d/lmod.sh
echo "Using arm-20 to build for astra"
# module swap gnu7 arm/20.1
module swap devpack-gnu7 devpack-arm/20200529

# Clang
#CC=mpicc
#CFLAGS="-g -O2"
#LDFLAGS="-g -O2"

# gcc
#CC=mpicc
#CFLAGS=-O
#LDFLAGS=-O

LIBS=

# Always make clean because the code builds so quickly. 
make clean

# make CC=mpicc CFLAGS="-g -O2" LDFLAGS="-g -O2" fftpf3d
make CC=mpicc CFLAGS="-g -O1" LDFLAGS="-g -O1" fftpf3d


