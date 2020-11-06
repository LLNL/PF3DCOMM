#!/bin/ksh
#

umask 027

# re-initialize modules to purge stuff from the parent shell
. /opt/atse/lmod/lmod/init/bash
# /etc/profile.d/lmod.sh
echo "Using devpack-arm to run on astra/stria"
# module swap gnu7 arm/20.1
module swap devpack-gnu7 devpack-arm/20200529
# module load totalview

export NODES=$SLURM_NNODES
# Output from the following command will be voluminous for large runs
echo $SLURM_NODELIST

##############################################################################
# Run the program

# This is typically run with a command like:
# sbatch -N 4  --account XXXXXX  -t 20:00 run-4x4-astra-direct-one-48p.bash

echo "Starting fft message passing test"

date

# WARNING!! Manually set the number of processes
export TYPE=direct
export PPN=48
let tasks=${NODES}*${PPN}
let ppn_socket=`expr $PPN/2`
mpirun -n ${tasks} --npersocket ${ppn_socket} --bind-to core ./fftpf3d 4 4 1024 768 20 100 ${NODES} ${TYPE} > log_stria_${NODES}_node_${TYPE}_${PPN}_proc.txt

date

exit
