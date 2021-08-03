
git clone https://github.com/slembcke/Chipmunk2D.git
cd Chipmunk*
git checkout Chipmunk-7.0.1 
cd ..
mv Chipmunk* Chipmunk-7.0.1 

mkdir -p obj_c99
make
make archive

cd fenics
./ufl*.sh
cd ..

mkdir build
cd build
cmake ..
make

./eQ test
