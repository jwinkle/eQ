module purge

# module load cmake fenics PETSc MPFR/4.0.2-GCCcore-6.4.0

# module load FENICS/2019.1.0-foss-2019b-Python-3.7.4

module load fenics/2019.1.0-intel-2019-Python-3.8


echo "module list: $(module list)"

echo "PETSC_DIR = $PETSC_DIR"

cd ..
touch ./CMakeLists.txt  #re-scans git branch and hash
touch ./src/inputOutput.cpp

cd ./build

make

