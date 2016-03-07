#!/bin/bash
#
# Convenient script to run the program IMCMC run
# But use the scripts in examples/ instead to get
# started
#
###################################################

# Recompile all ('yes' or 'no' or 'justWhatHasChanged')
MAKEALL='justWhatHasChanged'
# MPI ('yes' or 'no')
MPI='no'
# number of processes
NPROC=10

##################################################

# Compiles executables in root directory
if [ $MAKEALL == 'yes' ]; then 
  echo
  echo "Compilation and linking in progress ..."
  echo
  make clean
  if [ $MPI == 'yes' ]; then 
    make PAR=yes
  else
    make
  fi
  echo
  echo "Compilation and linking done !"
  echo
fi

if [ $MAKEALL == 'justWhatHasChanged' ]; then 
  echo
  echo "Compilation and linking in progress ..."
  echo
  if [ $MPI == 'yes' ]; then 
    make PAR=yes
  else
    make
  fi
  echo
  echo "Compilation and linking done !"
  echo
fi

# Runs program
echo
echo "Starting Program..."
echo
if [ $MPI == 'yes' ]; then 
  mpiexec -np $NPROC ./bin/RealisticDataMCMC # or mpirun -n $NPROC
else
  ./bin/RealisticDataMCMC
fi

echo
echo "Done!"
echo
echo "See results in directory: OUTPUT_FILES/ using the tools in utils"
date


