# diy_tests
Steps to reproduce:

git clone --recursive https://github.com/mclarsen/diy_tests.git

cd diy_tests

mkdir build

cd build 

cmake ..

srun -n 189 bin/reduce 

Fatal error in PMPI_Isend: Invalid rank, error stack:

PMPI_Isend(160): MPI_Isend(buf=0x817560, count=16, MPI_BYTE, dest=189, tag=0, MPI_COMM_WORLD, request=0x7fffffffb810) failed

PMPI_Isend(108): Invalid rank has value 189 but must be nonnegative and less than 189

srun -n 10 bin/all_to_all 

srun: error: rztopaz12: task 8: Segmentation fault (core dumped)
