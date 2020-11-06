#!/bin/bash
#

umask 027

########################################################################
# IBM Power 9
#   bsub -x -nnodes 2 -Is -XF -W 60 -G guests -R "span[ptile=20]" /usr/bin/tcsh

echo "Starting FFT message passing test"

date

# WARNING!! Manually set the number of nodes and processes
export TYPE=direct
export NODES=2
export PPN=32
let tasks=${NODES}*${PPN}
lrun -n ${PPN} -N ${NODES} ./fftpf3d 4 4 1024 768 20 100 ${NODES} ${TYPE} > log_pwr9_blueos_${NODES}_node_${TYPE}_${PPN}_proc.txt

date

exit
