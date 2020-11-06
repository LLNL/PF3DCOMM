#!/bin/bash
#

umask 027

export NODES=$SLURM_NNODES
# Output from the following command will be voluminous for large runs
echo $SLURM_NODELIST

# This is typically run like this on 36 core Intel Broadwell nodes
# sbatch -N 4 -n 144 -t 45 --exclusive run-4x4-toss3-cyclic.bash

echo "Starting fft message passing test"

date

# WARNING!! Manually set the number of nodes and processes
export TYPE=cyclic
export PPN=36
let tasks=${NODES}*${PPN}
srun -N ${NODES} --ntasks-per-node=${PPN} -n ${tasks} ./fftpf3d 4 4 1024 768 20 100 ${NODES} ${TYPE} > log_x86_64_toss3_${NODES}_node_${TYPE}_${PPN}_proc.txt

##############################################################################
# All done

date

exit
