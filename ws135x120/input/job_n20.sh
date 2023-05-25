#!/bin/bash
#SBATCH --account=clidyn.clidyn
#SBATCH -J ws00
#SBATCH -n 20
#SBATCH -p smp
##SBATCH --qos=12h
##SBATCH -t 12:00:00
#SBATCH --qos=30min
#SBATCH -t 00:30:00
##SBATCH --tasks-per-node=32
#SBATCH --cpus-per-task=6
#SBATCH -o ./output.txt

# list of hosts that you are running on
hostlist=$(scontrol show hostnames | tr '\n' ',' | rev | cut -c 2- | rev)
echo "hosts: $hostlist"
umask 022
ulimit -s
# maximum possible stacksize
#ulimit -s 1048576
#ulimit -s 1500000
ulimit -s unlimited
echo "ulimit -s"
ulimit -s

module purge
module load intel-oneapi-compilers load intel-oneapi-mpi
module load netcdf-fortran/4.5.4-oneapi2022.1.0
module list

echo "SLURM_CPUS_PER_TASK: $SLURM_CPUS_PER_TASK"
export SRUN_CPUS_PER_TASK=$SLURM_CPUS_PER_TASK

# even though we do not run an OpenMP code it is still a good idea to
# always set this
export OMP_NUM_THREADS=1
#export OMP_PROC_BIND=close
cd ${SLURM_SUBMIT_DIR}

ntx=1
nty=1
cat > eedata<<EOF
 &EEPARMS
 nTx=${ntx},
 nTy=${nty},
 printMapIncludesZeros = .TRUE.
 &
EOF

bdir="/isibhv/projects/p_sochic/SOCHIC_WG6_PACKAGE/LOCAL"
idir="../input"
ln -s ${idir}/* . >/dev/null 2>&1
ln -s ${bdir}/input/* . >/dev/null 2>&1
ln -s ${bdir}/domain_data/*.bin . >/dev/null 2>&1
ln -s ${bdir}/bathymetry/* . >/dev/null 2>&1
ln -s ${bdir}/initial_ECCO/pickup* . >/dev/null 2>&1

echo "SLURM_NTASKS             $SLURM_NTASKS"
echo "OMP_NUM_THREADS          $OMP_NUM_THREADS"
echo "OMP_PROC_BIND            $OMP_PROC_BIND"

olditer0=`grep nIter0 data | awk -F= '{print $NF}' | awk -F, '{print $1}'`
newiter0=$olditer0

# mpi specific environment variables:
export OMPI_MCA_osc="ucx"
export OMPI_MCA_pml="ucx"
export OMPI_MCA_btl="self"
export UCX_HANDLE_ERRORS="bt"
export OMPI_MCA_pml_ucx_opal_mem_hooks=1
#
export OMPI_MCA_io="romio321"          # basic optimisation of I/O
export UCX_TLS="shm,rc_mlx5,rc_x,self" # for jobs using LESS than 150 nodes
export UCX_UNIFIED_MODE="y"            # JUST for homogeneous jobs on CPUs, do not use for GPU nodes

# intel MPI
export I_MPI_PMI=pmi
export I_MPI_PMI_LIBRARY=/usr/lib64/libpmi.so

exe=../build/mitgcmuv
sopt="--cpu_bind=cores"
#sopt="--cpu_bind=cores --distribution=block:block"

echo "run the model starting at iter=$olditer0 with"
echo "srun ${sopt} $exe"
#srun --cpu_bind=cores --distribution=block:block $exe
srun ${sopt} $exe

\cp STDOUT.0000 stdout.0000

#\rm -rf STD*
#find . -type l -delete

sstat -j $SLURM_JOB_ID.batch --format=JobID,MaxVMSize
