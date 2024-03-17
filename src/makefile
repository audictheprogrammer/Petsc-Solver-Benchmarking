export PETSC_DIR=/mnt/c/Users/xu/petsc
export wPETSC_DIR=/mnt/c/Users/xu/petsc/share/petsc/datafiles/matrices
export wPETSC_DIR=/mnt/c/Users/xu/petsc/share/petsc/datafiles/matrices/MYMAT
export wPETSC_DIR=../data
-include ../../../petscdir.mk

CLEANFILES   = bench_spmv

ex72: ex72.o mmloader.o mmio.o
bench_spmv: bench_spmv.o mmloader.o mmio.o cJSON.o

include ${PETSC_DIR}/lib/petsc/conf/variables
include ${PETSC_DIR}/lib/petsc/conf/rules
