**PF3DCOMM**
=======

PF3DCOMM is an MPI benchmark that tests the performance of communication
patterns used in pF3D, a laser-plasma simulation code developed at LLNL.
The benchmark is intended for use in evaluating the effectiveness
of HPC interconnects.

Quick Start
-----------

The pF3DCOMM repo contains build and run scripts for HPC systems at LLNL.
To get started, make a copy of a build and a run script and modify
them to match your system. You can type something like this to run
on an LLNL Broadwell cluster:

```shell
# build the code
./build-llnl.kshsh

# allocate a node to run on
salloc -N 2 -n 72 --exclusive

# run the benchmark to measure performance
./run-4x4-toss3-one.bash
```

The next section briefly describes pF3D and the message passing
patterns it uses.
If you are just interested in building and running the benchmark,
jump down to the Building and Running sections.


Overview of pF3D
-----------

pF3D is used to simulate the interaction between a high intensity laser
and a plasma (ionized gas) in experiments performed at LLNL's National
Ignition Facility. pF3D simulations can consume large amounts of computer
time, so it is important that it is well optimized.

pF3D normally haveuses float precision. pF3D decomposes the grid
in x, y, and z. An MPI domain has nxl-by-nyl-by-nzl zones
and there are mp_p-b-mp_q-by-mp_r MPI domains. 

The PF3DCOMM benchmark mimics three message passing patterns used by pF3D.

1) Some physics modules in pF3D use FFTs ("spectral methods").
The laser propagates in the z-direction, so pF3D usesperforms
2D FFTs over xy-planes. MPI messages are used to pass data
between processes sharing an xy-plane. FFT messages account
for most of the message bytes sent by pF3D in typical runs.

pF3D uses MPI_Alltoall to pass messages and operates
on 1D communicators when using the default message passing
package. There are mp_q*mp_r 1D x-communicators and
mp_p*mp_r 1D y-communicators. Some MPI implementations have
special optimizations for MPI_Alltoall operating
on MPI_commworld. pF3DCOMM canmay benefit from MPI_Alltoall
optimizations on small communicators. 

PF3DCOMM has a message pattern that mimics 2D FFT communication
in pF3D, but it does not perform any FFT computation. 

2) Advection of light in the z-direction is responsible for the second
highest number of message bytes. This pattern involves the exchange
of 2 or 3 planes between nearest neighbors in z.

3) Hydrodynamics is treated using a finite difference approach. The message 
pattern is a standard hydrodynamics halo exchange with the 6 nearest
neighbor MPI domains.


Performance Optimization
-----------

The choice of MPI decomposition can have a strong impact on the
performance of pF3D and PF3DCOMM. Messages between processes on the
same node can be passed via shared memory. The DDR bandwidth
of a node is usually significantly higher than the per node
interconnect bandwidth. Performance is improved if communicating
processes are on the same node. 2D FFT performance can be improved
by decomposing so that nodes have complete xy-planes.
The 2D FFT performance on systems with GPUs can be improved
by fitting an entire xy-plane into a single GPU because 2D FFT
data motion happens in the high-bandwidth GPU memory. There
may not be enough memory on the GPU to fit all the necessary arrays. 
The performance of the advection messages
becomes a bottleneck when nzl is small. Using a 1-by-1
decompostion in xy leads to small values for nzl or else
restricts a run to using a small number of nodes. 
pF3DCOMM makes it easy to explore the performance impact
of different MPI decompositions.


Building PF3DCOMM
-----------

PF3DCOMM has very few external dependencies so it is easy
to build. An MPI wrapper (e.g. mpicc) for the C compiler is assumed
to be available. If it isn't, you will need to add the appropriate
include files to the compile lines and the appropriate libraries
to the load line.

Sample build scripts for LLNL systems (build-llnl.ksh) and for
an ARM cluster at Sandia National Laboratory (build-astra.bash)
are provided. 


Running PF3DCOMM
-----------

A single run of the benchmark takes a few minutes. 

Run scripts using different mappings between MPI rank and node
are provided to make it easy to check the impact of shread memory
message passing. Some scripts run a single test. Other scripts
run a weak scaling study. 

The scripts for "toss3" are intended for 36 core Intel Broadwell nodes.
The scripts for "blueos" are intended for 40 core IBM Power 9 nodes.
The scripts for "astra" are intended for 56 core ARM nodes with
Cavium Thunder X2 processors. 

The build and run scripts rely on the SYS_TYPE variable defined
on all LLNL clusters. To create your own build and run scripts, start with
the sample scripts and adjust them to match your system.


License
-----------

PF3DCOMM is distributed under the terms of the BSD-3 license. All new contributions
must be made under the BSD-3 license.

See LICENSE and NOTICE for details.

SPDX-License-Identifier: BSD-3-Clause

LLNL-CODE-815620
