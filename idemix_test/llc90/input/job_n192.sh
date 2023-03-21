#!/bin/bash
#SBATCH --account=clidyn.clidyn
#SBATCH -J llc16
#SBATCH -n 192
#SBATCH -p mpp
#SBATCH --qos=12h
##SBATCH --qos=30min
#SBATCH -t 12:00:00
##SBATCH -t 00:30:00
#SBATCH --tasks-per-node=64
#SBATCH --cpus-per-task=2
##SBATCH --cpus-per-task=1
#SBATCH -o ./output.txt

echo "SLURM_CPUS_PER_TASK: $SLURM_CPUS_PER_TASK"
export SRUN_CPUS_PER_TASK=$SLURM_CPUS_PER_TASK

module purge
module load intel-oneapi-compilers
# module load intel-oneapi-mkl
module load intel-oneapi-mpi
module load netcdf-fortran/4.5.4-oneapi2022.1.0
module list

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
# even though we do not run an OpenMP code it is still a good idea to
# always set this
export OMP_NUM_THREADS=1
#export OMP_PROC_BIND=close
# very important if you want to run on more than one node!!!
#export I_MPI_FABRICS=shm:tmi
cd ${SLURM_SUBMIT_DIR}
echo ${SLURM_NTASKS}

ntx=1
nty=1
cat > eedata<<EOF
 &EEPARMS
 nTx=${ntx},
 nTy=${nty},
 useCubedSphereExchange = .TRUE.,
 printMapIncludesZeros = .TRUE.
 &
EOF

bdir="/albedo/home/mlosch/MITgcm/nils/llc90"
idir="../input"
ln -s ${idir}/* . >/dev/null 2>&1
ln -s ${bdir}/input/* . >/dev/null 2>&1
#ln -s ${bdir}/build_fwd/mitgcmuv . >/dev/null 2>&1
#cp ${bdir}/build_fwd_192/mitgcmuv . >/dev/null 2>&1

exe=./mitgcmuv

echo "SLURM_NTASKS             $SLURM_NTASKS"
echo "OMP_NUM_THREADS          $OMP_NUM_THREADS"
echo "OMP_PROC_BIND            $OMP_PROC_BIND"

olditer0=`grep nIter0 data | awk -F= '{print $NF}' | awk -F, '{print $1}'`
newiter0=$olditer0

#echo "copy forcing data to local disk"
#rsync -aur /albedo/work/projects/p_forcing_ieee/JRA55-do-v1.4.0_bin /tmp/
#echo "... done"

sopt="--cpu_bind=cores"
#sopt="--cpu_bind=cores --distribution=block:block"

echo "run the model starting at iter=$olditer0 with"
echo "srun ${sopt} $exe"
#srun --cpu_bind=cores --distribution=block:block $exe
srun ${sopt} $exe

lastline=`tail -1 STDOUT.0000`
echo $lastline
testline=`echo $lastline  | grep 'Execution ended Normally'`
echo $testline

echo "looking for pickups"

# filesize=381607200
# last_pickup=`ls -lt pickup.*.data | awk '$5~/'${filesize}'/{print $NF}'|head -1 `
# we do not need the filesize, because we make sure that the job terminates cleanly
# so it is enough to make sure that we link the latest pickup
last_pickup=`ls -t pickup.*.data | head -1 `
if [ -n "${last_pickup}" ]; then
    echo "linking last_pickup ${last_pickup}"
    last_meta=`echo $last_pickup | sed 's/data/meta/'`
    timestep=`awk '$1~/timeStepNumber/{printf("%010d",$4)}' $last_meta`
    \ln -sf ${last_pickup} pickup.${timestep}
    last_pickupmeta=`echo $last_pickup|awk '{split($0,a,".");print a[1]"."a[2]".meta"}'`
    \ln -sf ${last_pickupmeta} pickup.${timestep}.meta
fi

# filesize=10108800
# last_pickup_seaice=`ls -lt pickup_seaice.*.data | awk '$5~/'${filesize}'/{print $NF}'|head -1 `
last_pickup_seaice=`ls -t pickup_seaice.*.data | head -1 `
if [ -n "${last_pickup_seaice}" ]; then
    echo "linking last_pickup_seaice ${last_pickup_seaice}"
    last_meta=`echo $last_pickup_seaice | sed 's/data/meta/'`
    timestep=`awk '$1~/timeStepNumber/{printf("%010d",$4)}' $last_meta`
    \ln -sf ${last_pickup_seaice} pickup_seaice.${timestep}
    last_pickupmeta=`echo $last_pickup_seaice|awk '{split($0,a,".");print a[1]"."a[2]".meta"}'`
    \ln -sf ${last_pickupmeta} pickup_seaice.${timestep}.meta
fi

# filesize=42120000
# last_pickup_ggl90=`ls -lt pickup_ggl90.*.data | awk '$5~/'${filesize}'/{print $NF}'|head -1 `
last_pickup_ggl90=`ls -t pickup_ggl90.*.data | head -1 `
if [ -n "${last_pickup_ggl90}" ]; then
    echo "linking last_pickup_ggl90 ${last_pickup_ggl90}"
    last_meta=`echo $last_pickup_ggl90 | sed 's/data/meta/'`
    timestep=`awk '$1~/timeStepNumber/{printf("%010d",$4)}' $last_meta`
    \ln -sf ${last_pickup_ggl90} pickup_ggl90.${timestep}
    last_pickupmeta=`echo $last_pickup_ggl90|awk '{split($0,a,".");print a[1]"."a[2]".meta"}'`
    \ln -sf ${last_pickupmeta} pickup_ggl90.${timestep}.meta
fi

# save old data and STDOUT
\cp data data.${olditer0}
\cp STDOUT.0000 stdout_rerun.${olditer0}

if [[ "x$timestep" == "x" ]]; then
    echo "did not find appropriate pickup files"
    exitcode=1
else
    # adjust the startpoint in data
    newiter0=`expr $timestep \* 1`
    echo "new niter0 = "$newiter0
fi

lastline=`tail -1 STDOUT.0000`
echo "last line of STDOUT.0000"
echo $lastline
testline=`echo $lastline | grep '======================================================='`
if test "x$testline" = "x" ; then
    # check if we have terminated properly
    testline=`echo $lastline | grep 'Execution ended Normally'`
else
    # this usually means that we ran out of forcing, so it is time to
    # start a new cycle. For this we need to change the start date in
    # data.cal
    #
    # startDate_1=19580101, for the first cycle
    # startDate_1=18960101, for the second cycle
    # startDate_1=18340101, for the third cycle
    # startDate_1=17720101, for the fourth cycle
    # startDate_1=17100101, for the fifth cycle
    nyears=`\ls diags2D.*.meta | wc -l`
    echo "Found $nyears of output,"
    if (( $nyears == 62 )); then
	ncycle=2
	startdate=18960101
    elif (( $nyears == 124 )); then
	ncycle=3
	startdate=18340101
    elif (( $nyears == 186 )); then
	ncycle=4
	startdate=17720101
	echo "it is time for cycle 4"
    elif (( $nyears == 248 )); then
	ncycle=5
	startdate=17100101
	echo "it is time for cycle 5"
    else
	echo "number years not consistent with expectation"
	ncycle=0
	startdate=000000000
        echo "triggering termination by resetting testline"
        testline=
    fi
    echo "it is time for cycle $ncycle, setting start date to $startdate"
    oldcycle="cycle_$((${ncycle}-1))"
    \cp data.cal data.cal.${oldcycle}
    cat data.cal.${oldcycle} | sed "/ startDate_1/s/.*/ startDate_1=${startdate},/" >| data.cal
fi

if test "x$testline" = "x" ; then
    echo "something is wrong; exiting"
    exitcode=1
else
    # 42 years with deltaT = 1800sec : 1956528000./1800. = 1086960 timesteps
    # if (( $newiter0 < 1086960 )); then
        cat data.${olditer0} | sed "/ nIter0/s/.*/ nIter0=${newiter0},/" >| data
        \rm STD*
        echo "starting next job"
        exitcode=0

        echo "sbatch $0"
        sbatch $0
    # else
    #     echo "simulation reached the end, cleaning up"
    #     find . -type l -delete
    #     \rm -rf *.log
    # fi
fi

sstat -j $SLURM_JOB_ID.batch --format=JobID,MaxVMSize
