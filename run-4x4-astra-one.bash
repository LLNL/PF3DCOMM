#!/bin/bash
#

umask 027

# re-initialize modules to purge stuff from the parent shell
. /opt/atse/lmod/lmod/init/bash
# /etc/profile.d/lmod.sh
echo "Using devpack-arm to run on  astra/stria"
# module swap gnu7 arm/20.1
module swap devpack-gnu7 devpack-arm/20200529

export NODES=$SLURM_NNODES
echo $SLURM_NODELIST

##############################################################################
# Run the program

# Astra ARM Thunder X2
# This is typically run after an salloc command.
# salloc --account=fy180039 -N 2 -n 112 -t 30 --exclusive
# 48 nodes are 2 edge switches on stria
# salloc -N 48 -t 48:00:00 --account fy180039 --reservation=shlange

echo "Starting fft message passing test"

date

export TYPE=direct
export PPN=48
let tasks=${NODES}*${PPN}
# srun -N ${NODES} --ntasks-per-node=${PPN} -n ${tasks} ./fftpf3d 4 4 1024 768 20 100 ${NODES} ${TYPE} > log_stria_${NODES}_node_${TYPE}_${PPN}_proc.txt
mpirun -n ${tasks} --npersocket ${ppn_socket} --bind-to core ./fftpf3d 4 4 1024 768 20 100 ${NODES} ${TYPE} > log_stria_${NODES}_node_${TYPE}_${PPN}_proc.txt

date
export PPN=56
let tasks=${NODES}*${PPN}
# srun -N ${NODES} --ntasks-per-node=${PPN} -n ${tasks} ./fftpf3d 4 4 1024 768 20 100 ${NODES} ${TYPE} > log_stria_${NODES}_node_${TYPE}_${PPN}_proc.txt
mpirun -n ${tasks} --npersocket ${ppn_socket} --bind-to core ./fftpf3d 4 4 1024 768 20 100 ${NODES} ${TYPE} > log_stria_${NODES}_node_${TYPE}_${PPN}_proc.txt

date
export TYPE=direct
let tasks=${NODES}*${PPN}
# srun -N ${NODES} --ntasks-per-node=${PPN} -n ${tasks} --distribution=cyclic ./fftpf3d 4 4 1024 768 20 100 ${NODES} ${TYPE} > log_stria_${NODES}_node_${TYPE}_${PPN}_proc.txt
mpirun -n ${tasks} --npersocket ${ppn_socket} --bind-to core --rank-by node ./fftpf3d 4 4 1024 768 20 100 ${NODES} ${TYPE} > log_stria_${NODES}_node_${TYPE}_${PPN}_proc.txt

##############################################################################
# All done

date

exit
