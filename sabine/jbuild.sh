module purge
# module load cmake fenics PETSc MPFR/4.0.2-GCCcore-6.4.0

module load FENICS 

echo "module list"
module list


echo $PETSC_DIR

touch ~/eQ/CMakeLists.txt  #re-scans git branch and hash
touch ~/eQ/src/inputOutput.cpp

cd ~/eQ/build

make

