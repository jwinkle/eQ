# eQ     &nbsp;  :Agent-based modeling for E.coli

  
eQ is a C++ software for 2D Agent-Based Modeling (ABM) of rod-shaped bacteria growing in microfluidic traps with an integrated 
high-resolution quorum sensing (QS) PDE solver using the finite-element software _Fenics_.

## Software overview:

1. programmable growth laws with mechanical inhibition feedback using the Chipmunk 2D physics engine: https://chipmunk-physics.net/
2. integrated small-molecule diffusion solver with implicit numerics using Fenics finite-element method (FEM) computing platform: https://fenicsproject.org/
3. data recording to JSON format (https://github.com/nlohmann/json), 
    with decoding scripts in Matlab (currently) to plot time series and colony dynamics images
4. MPI parallel model with multiple QS layers and physics engine time-stepped on separate cores
    
## Installation:

The simplest option to install and run eQ requires installing the Docker desktop engine to easily install Fenics.

### Installation via Docker

###  Overview

1. Install the Docker desktop for Mac, Windows, Linux:  https://www.docker.com/products/docker-desktop
2. Install the latest _Fenics_ image as found here: https://fenics.readthedocs.io/projects/containers/en/latest/introduction.html#running-fenics-in-docker
3. Run the ```eQinstall.sh``` script to fetch/compile/install dependencies and to compile eQ and test eQ within the Docker image for _Fenics_
4. Checkout an example branch (e.g. the latest PLoS article simulation code) to run and verify functionality

### Detailed Installation Instructions

#### I. INSTALL DOCKER DESKTOP
* Follow instructions here to install Docker Desktop Engine:
https://www.docker.com/products/docker-desktop



#### II. INSTALL FENICS IMAGE INTO DOCKER DESKTOP
* Create a convenient directory to install the Fenics/stable image and to clone eQ repo.
For example, create a folder: ```docker```
``` bash
winkle$: pwd 
/Users/winkle
winkle$: mkdir docker; cd docker; pwd
/Users/winkle/docker
```

* See detailed (current) instructions here: https://fenics.readthedocs.io/projects/containers/en/latest/introduction.html#running-fenics-in-docker
* The ```docker run ...``` command can be run in a terminal or via the CLI button on Docker Desktop once the image is installed into Docker Desktop
* Example when run from Mac terminal:
``` bash
winkle$:docker run -ti -v $(pwd):/home/fenics/shared -w /home/fenics/shared quay.io/fenicsproject/stable:current
Unable to find image 'quay.io/fenicsproject/stable:current' locally
current: Pulling from fenicsproject/stable
c64513b74145: Pull complete 
01b8b12bad90: Pull complete 
c5d85cf7a05f: Pull complete 
b6b268720157: Pull complete 
e12192999ff1: Pull complete 
d39ece66b667: Pull complete 
65599be66378: Pull complete 
04de8bc2d500: Pull complete 
abb684b96e3d: Pull complete 
5fb302170f66: Pull complete 
56d9f5e23832: Pull complete 
2362411179ee: Pull complete 
f0f1cc16c840: Downloading [=======>                                           ]  37.74MB/255.6MB
c5bca2c84e5e: Downloading [===============================>                   ]  50.51MB/80.97MB
6f80705d1d37: Download complete 
91e158461d4c: Download complete 
4b299c049e4f: Download complete 
4642a6d46aeb: Downloading [==============================================>    ]   64.5MB/68.66MB
0526f57f6bb3: Waiting 
6bda00aab163: Waiting 
599370ebdfdb: Waiting 
3debe21df7b8: Waiting 
bd6da2b7ad00: Waiting 
3111ae43df8e: Waiting 
dd2003bd419e: Waiting 
49e25fdb34c3: Waiting 
d95fd6a8326e: Waiting 
```

...  [after finishing:] ...

``` bash
c64513b74145: Pull complete 
01b8b12bad90: Pull complete 
c5d85cf7a05f: Pull complete 
b6b268720157: Pull complete 
e12192999ff1: Pull complete 
d39ece66b667: Pull complete 
65599be66378: Pull complete 
04de8bc2d500: Pull complete 
abb684b96e3d: Pull complete 
5fb302170f66: Pull complete 
56d9f5e23832: Pull complete 
2362411179ee: Pull complete 
f0f1cc16c840: Pull complete 
c5bca2c84e5e: Pull complete 
6f80705d1d37: Pull complete 
91e158461d4c: Pull complete 
4b299c049e4f: Pull complete 
4642a6d46aeb: Pull complete 
0526f57f6bb3: Pull complete 
6bda00aab163: Pull complete 
599370ebdfdb: Pull complete 
3debe21df7b8: Pull complete 
bd6da2b7ad00: Pull complete 
3111ae43df8e: Pull complete 
dd2003bd419e: Pull complete 
49e25fdb34c3: Pull complete 
d95fd6a8326e: Pull complete 
Digest: sha256:7aa3dad185d43fcf3d32abdf62fec86dab31b5159d66b2228275edc702181bb1
Status: Downloaded newer image for quay.io/fenicsproject/stable:current
# FEniCS stable version image

Welcome to FEniCS/stable!

This image provides a full-featured and optimized build of the stable
release of FEniCS.

To help you get started this image contains a number of demo
programs. Explore the demos by entering the 'demo' directory, for
example:

    cd ~/demo/python/documented/poisson
    python3 demo_poisson.py
fenics@312b2d51ad9b:~/shared$

fenics@312b2d51ad9b:~/shared$ pwd
/home/fenics/shared

```
You should now be able to see the image in **RUNNING** mode by opening Docker Desktop Application:
![](/images/dockerContainers.png)


**Notes:**  
* The Container name seen here ```ecstatic_ishizaka``` is auto-generated and used as a mnemonic for the hash of the container,
and will be different for each user/instance.
See an explanation here:  https://frightanic.com/computers/docker-default-container-names/

* The ```-v``` option installs a shared volume with the Docker image;
  you can now transfer files as: ```/home/fenics/shared``` in the Docker image and (here) ```/Users/winkle/docker``` map to the same folder.
  
  
  
#### III. CLONE THE eQ REPOSITORY INTO THE DOCKER/FENICS IMAGE
* Here, the ```CLI``` was used and ```bash``` was run, which gives the different prompt used below.
* Clone the eQ github repository by running: ``` git clone https://github.com/jwinkle/eQ.git```, as below
* ```cd``` into the eQ folder, then update the permissions to the shell script ```chmod +x eQsetup1.sh```
``` bash
root@e85404d46ca8:/home/fenics/shared# pwd
/home/fenics/shared
root@e85404d46ca8:/home/fenics/shared# ls
root@e85404d46ca8:/home/fenics/shared# git clone https://github.com/jwinkle/eQ.git
Cloning into 'eQ'...
remote: Enumerating objects: 1150, done.
remote: Counting objects: 100% (38/38), done.
remote: Compressing objects: 100% (32/32), done.
remote: Total 1150 (delta 15), reused 10 (delta 4), pack-reused 1112
Receiving objects: 100% (1150/1150), 8.07 MiB | 5.86 MiB/s, done.
Resolving deltas: 100% (842/842), done.

root@e85404d46ca8:/home/fenics/shared# ls -l
total 0
drwxr-xr-x 15 root root 480 Aug  3 18:38 eQ

root@e85404d46ca8:/home/fenics/shared# cd eQ;  ls -al
total 76
drwxr-xr-x 17 root root   544 Aug  3 18:43 .
drwxr-xr-x  4 root root   128 Aug  3 18:38 ..
-rw-r--r--  1 root root  6148 Aug  3 18:38 .DS_Store
drwxr-xr-x 13 root root   416 Aug  3 18:38 .git
-rw-r--r--  1 root root   558 Aug  3 18:38 .gitignore
-rw-r--r--  1 root root  2983 Aug  3 18:38 CMakeLists.txt
-rw-r--r--  1 root root  1072 Aug  3 18:38 LICENSE
-rw-r--r--  1 root root   551 Aug  3 18:38 Makefile
-rw-r--r--  1 root root 37564 Aug  3 18:38 diffuclass.cpp
-rw-r--r--  1 root root  3426 Aug  3 18:38 diffuclass.h
-rw-r--r--  1 root root   241 Aug  3 18:43 eQsetup1.sh
drwxr-xr-x 20 root root   640 Aug  3 18:38 fenics
drwxr-xr-x  7 root root   224 Aug  3 18:38 matlab
drwxr-xr-x  3 root root    96 Aug  3 18:38 nlohmann
drwxr-xr-x 16 root root   512 Aug  3 18:38 qtensor
drwxr-xr-x 18 root root   576 Aug  3 18:38 src
-rw-r--r--  1 root root   231 Aug  3 18:38 version.h.in

root@e85404d46ca8:/home/fenics/shared/eQ# chmod +x eQsetup1.sh 
root@e85404d46ca8:/home/fenics/shared/eQ# ls -al
total 76
drwxr-xr-x 17 root root   544 Aug  3 18:43 .
drwxr-xr-x  4 root root   128 Aug  3 18:38 ..
-rw-r--r--  1 root root  6148 Aug  3 18:38 .DS_Store
drwxr-xr-x 13 root root   416 Aug  3 18:38 .git
-rw-r--r--  1 root root   558 Aug  3 18:38 .gitignore
-rw-r--r--  1 root root  2983 Aug  3 18:38 CMakeLists.txt
-rw-r--r--  1 root root  1072 Aug  3 18:38 LICENSE
-rw-r--r--  1 root root   551 Aug  3 18:38 Makefile
-rw-r--r--  1 root root 37564 Aug  3 18:38 diffuclass.cpp
-rw-r--r--  1 root root  3426 Aug  3 18:38 diffuclass.h
-rwxr-xr-x  1 root root   241 Aug  3 18:43 eQsetup1.sh
drwxr-xr-x 20 root root   640 Aug  3 18:38 fenics
drwxr-xr-x  7 root root   224 Aug  3 18:38 matlab
drwxr-xr-x  3 root root    96 Aug  3 18:38 nlohmann
drwxr-xr-x 16 root root   512 Aug  3 18:38 qtensor
drwxr-xr-x 18 root root   576 Aug  3 18:38 src
-rw-r--r--  1 root root   231 Aug  3 18:38 version.h.in
```


* run the shell script;  there is a lot of output to compile files, shortened version follows
``` bash
root@e85404d46ca8:/home/fenics/shared/eQ# ./eQsetup1.sh 
Cloning into 'Chipmunk2D'...
remote: Enumerating objects: 11552, done.
remote: Counting objects: 100% (10/10), done.
remote: Compressing objects: 100% (8/8), done.
remote: Total 11552 (delta 0), reused 4 (delta 0), pack-reused 11542
Receiving objects: 100% (11552/11552), 6.07 MiB | 2.02 MiB/s, done.
Resolving deltas: 100% (7013/7013), done.
Note: checking out 'Chipmunk-7.0.1'.

You are in 'detached HEAD' state. You can look around, make experimental
changes and commit them, and you can discard any commits you make in this
state without impacting any branches by performing another checkout.

If you want to create a new branch to retain commits you create, you may
do so (now or later) by using -b with the checkout command again. Example:

  git checkout -b <new-branch-name>

HEAD is now at e98a7ce Updating version notes.
gcc -c -DNDEBUG -IChipmunk-7.0.1/include -Wall -O -std=c99 -pthread  -o obj_c99/cpCollision.o Chipmunk-7.0.1/src/cpCollision.c
gcc -c -DNDEBUG -IChipmunk-7.0.1/include -Wall -O -std=c99 -pthread  -o obj_c99/cpPolyline.o Chipmunk-7.0.1/src/cpPolyline.c
gcc -c -DNDEBUG -IChipmunk-7.0.1/include -Wall -O -std=c99 -pthread  -o obj_c99/cpRatchetJoint.o Chipmunk-7.0.1/src/cpRatchetJoint.c
gcc -c -DNDEBUG -IChipmunk-7.0.1/include -Wall -O -std=c99 -pthread  -o obj_c99/cpRotaryLimitJoint.o Chipmunk-7.0.1/src/cpRotaryLimitJoint.c
gcc -c -DNDEBUG -IChipmunk-7.0.1/include -Wall -O -std=c99 -pthread  -o obj_c99/cpSpaceDebug.o Chipmunk-7.0.1/src/cpSpaceDebug.c
gcc -c -DNDEBUG -IChipmunk-7.0.1/include -Wall -O -std=c99 -pthread  -o obj_c99/cpSpaceQuery.o Chipmunk-7.0.1/src/cpSpaceQuery.c
gcc -c -DNDEBUG -IChipmunk-7.0.1/include -Wall -O -std=c99 -pthread  -o obj_c99/cpSpaceComponent.o Chipmunk-7.0.1/src/cpSpaceComponent.c
gcc -c -DNDEBUG -IChipmunk-7.0.1/include -Wall -O -std=c99 -pthread  -o obj_c99/cpConstraint.o Chipmunk-7.0.1/src/cpConstraint.c
gcc -c -DNDEBUG -IChipmunk-7.0.1/include -Wall -O -std=c99 -pthread  -o obj_c99/cpSpaceHash.o Chipmunk-7.0.1/src/cpSpaceHash.c
gcc -c -DNDEBUG -IChipmunk-7.0.1/include -Wall -O -std=c99 -pthread  -o obj_c99/cpShape.o Chipmunk-7.0.1/src/cpShape.c
gcc -c -DNDEBUG -IChipmunk-7.0.1/include -Wall -O -std=c99 -pthread  -o obj_c99/cpMarch.o Chipmunk-7.0.1/src/cpMarch.c
gcc -c -DNDEBUG -IChipmunk-7.0.1/include -Wall -O -std=c99 -pthread  -o obj_c99/cpHastySpace.o Chipmunk-7.0.1/src/cpHastySpace.c
gcc -c -DNDEBUG -IChipmunk-7.0.1/include -Wall -O -std=c99 -pthread  -o obj_c99/chipmunk.o Chipmunk-7.0.1/src/chipmunk.c
gcc -c -DNDEBUG -IChipmunk-7.0.1/include -Wall -O -std=c99 -pthread  -o obj_c99/cpBody.o Chipmunk-7.0.1/src/cpBody.c
gcc -c -DNDEBUG -IChipmunk-7.0.1/include -Wall -O -std=c99 -pthread  -o obj_c99/cpSpace.o Chipmunk-7.0.1/src/cpSpace.c
gcc -c -DNDEBUG -IChipmunk-7.0.1/include -Wall -O -std=c99 -pthread  -o obj_c99/cpSimpleMotor.o Chipmunk-7.0.1/src/cpSimpleMotor.c
gcc -c -DNDEBUG -IChipmunk-7.0.1/include -Wall -O -std=c99 -pthread  -o obj_c99/cpBBTree.o Chipmunk-7.0.1/src/cpBBTree.c
gcc -c -DNDEBUG -IChipmunk-7.0.1/include -Wall -O -std=c99 -pthread  -o obj_c99/cpRobust.o Chipmunk-7.0.1/src/cpRobust.c
gcc -c -DNDEBUG -IChipmunk-7.0.1/include -Wall -O -std=c99 -pthread  -o obj_c99/cpArray.o Chipmunk-7.0.1/src/cpArray.c
gcc -c -DNDEBUG -IChipmunk-7.0.1/include -Wall -O -std=c99 -pthread  -o obj_c99/cpGearJoint.o Chipmunk-7.0.1/src/cpGearJoint.c
gcc -c -DNDEBUG -IChipmunk-7.0.1/include -Wall -O -std=c99 -pthread  -o obj_c99/cpPinJoint.o Chipmunk-7.0.1/src/cpPinJoint.c
gcc -c -DNDEBUG -IChipmunk-7.0.1/include -Wall -O -std=c99 -pthread  -o obj_c99/cpDampedRotarySpring.o Chipmunk-7.0.1/src/cpDampedRotarySpring.c
gcc -c -DNDEBUG -IChipmunk-7.0.1/include -Wall -O -std=c99 -pthread  -o obj_c99/cpHashSet.o Chipmunk-7.0.1/src/cpHashSet.c
gcc -c -DNDEBUG -IChipmunk-7.0.1/include -Wall -O -std=c99 -pthread  -o obj_c99/cpPivotJoint.o Chipmunk-7.0.1/src/cpPivotJoint.c
gcc -c -DNDEBUG -IChipmunk-7.0.1/include -Wall -O -std=c99 -pthread  -o obj_c99/cpGrooveJoint.o Chipmunk-7.0.1/src/cpGrooveJoint.c
gcc -c -DNDEBUG -IChipmunk-7.0.1/include -Wall -O -std=c99 -pthread  -o obj_c99/cpDampedSpring.o Chipmunk-7.0.1/src/cpDampedSpring.c
gcc -c -DNDEBUG -IChipmunk-7.0.1/include -Wall -O -std=c99 -pthread  -o obj_c99/cpSpatialIndex.o Chipmunk-7.0.1/src/cpSpatialIndex.c
gcc -c -DNDEBUG -IChipmunk-7.0.1/include -Wall -O -std=c99 -pthread  -o obj_c99/cpSlideJoint.o Chipmunk-7.0.1/src/cpSlideJoint.c
gcc -c -DNDEBUG -IChipmunk-7.0.1/include -Wall -O -std=c99 -pthread  -o obj_c99/cpPolyShape.o Chipmunk-7.0.1/src/cpPolyShape.c
gcc -c -DNDEBUG -IChipmunk-7.0.1/include -Wall -O -std=c99 -pthread  -o obj_c99/cpSpaceStep.o Chipmunk-7.0.1/src/cpSpaceStep.c
gcc -c -DNDEBUG -IChipmunk-7.0.1/include -Wall -O -std=c99 -pthread  -o obj_c99/cpSweep1D.o Chipmunk-7.0.1/src/cpSweep1D.c
gcc -c -DNDEBUG -IChipmunk-7.0.1/include -Wall -O -std=c99 -pthread  -o obj_c99/cpArbiter.o Chipmunk-7.0.1/src/cpArbiter.c
mkdir -p lib
ar rcs ./lib/libcp.a ./obj_c99/*.o 
This is FFC, the FEniCS Form Compiler, version 2019.1.0.post0.
UFC backend version 2018.1.0, signature 1e5e8d5476f82af05ab9c93365d99f971e0284d6.
For further information, visit https://bitbucket.org/fenics-project/ffc/.

Python 3.6.7 (default, Oct 22 2018, 11:32:17) 
[GCC 8.2.0] on linux

Compiling form AD1Dss

Compiler stage 1: Analyzing form(s)
-----------------------------------
  
  Geometric dimension:       1
  Number of cell subdomains: 0
  Rank:                      2
  Arguments:                 '(v_0, v_1)'
  Number of coefficients:    2
  Coefficients:              '[w_1, w_2]'
  Unique elements:           'CG1(?,?), R0(?,?), Vector<1 x CG1(?,?)>'
  Unique sub elements:       'CG1(?,?), R0(?,?), Vector<1 x CG1(?,?)>'
  
  representation:    auto --> uflacs
  quadrature_rule:   auto --> default
  quadrature_degree: auto --> 1
  quadrature_degree: 1
  
  Geometric dimension:       1
  Number of cell subdomains: 0
  Rank:                      1
  Arguments:                 '(v_0)'
  Number of coefficients:    1
  Coefficients:              '[w_0]'
  Unique elements:           'CG1(?,?), Vector<1 x CG1(?,?)>'
  Unique sub elements:       'CG1(?,?), Vector<1 x CG1(?,?)>'
  
  representation:    auto --> uflacs
  quadrature_rule:   auto --> default
  quadrature_degree: auto --> 2
  quadrature_degree: 2
  
Compiler stage 1 finished in 0.0328395 seconds.

Compiler stage 2: Computing intermediate representation
-------------------------------------------------------
  Computing representation of 3 elements
  Computing representation of 3 dofmaps
  Computing representation of 1 coordinate mappings
  Computing representation of integrals
  Computing uflacs representation
  Computing uflacs representation
  Computing representation of forms
  
Compiler stage 2 finished in 0.0315979 seconds.

Compiler stage 3: Optimizing intermediate representation
--------------------------------------------------------
  Optimizing uflacs representation
  Optimizing uflacs representation
  
Compiler stage 3 finished in 0.000434637 seconds.

Compiler stage 4: Generating code
---------------------------------
  Generating code for 3 finite_element(s)
  Generating code for 3 dofmap(s)
  Generating code for 1 coordinate_mapping(s)
  Generating code for integrals
  Generating code from ffc.uflacs representation
  Generating code from ffc.uflacs representation
  Generating code for forms
  
Compiler stage 4 finished in 0.0342863 seconds.

Compiler stage 4.1: Generating additional wrapper code
------------------------------------------------------
  Generating wrapper code for DOLFIN
  
Compiler stage 4.1 finished in 0.000582933 seconds.

Compiler stage 5: Formatting code
---------------------------------
  
Compiler stage 5 finished in 0.000532389 seconds.

FFC finished in 0.100855 seconds.
Output written to ./AD1Dss.h.
```

[ ... ]

``` bash
FFC finished in 0.0749705 seconds.
Output written to ./project.h.
-- The C compiler identification is GNU 7.4.0
-- The CXX compiler identification is GNU 7.4.0
-- Check for working C compiler: /usr/bin/cc
-- Check for working C compiler: /usr/bin/cc -- works
-- Detecting C compiler ABI info
-- Detecting C compiler ABI info - done
-- Detecting C compile features
-- Detecting C compile features - done
-- Check for working CXX compiler: /usr/bin/c++
-- Check for working CXX compiler: /usr/bin/c++ -- works
-- Detecting CXX compiler ABI info
-- Detecting CXX compiler ABI info - done
-- Detecting CXX compile features
-- Detecting CXX compile features - done
-- Boost version: 1.65.1
-- Found PkgConfig: /usr/bin/pkg-config (found version "0.29.1") 
-- Checking for one of the modules 'craypetsc_real;PETSc'
-- Looking for sys/types.h
-- Looking for sys/types.h - found
-- Looking for stdint.h
-- Looking for stdint.h - found
-- Looking for stddef.h
-- Looking for stddef.h - found
-- Check size of PetscInt
-- Check size of PetscInt - failed
-- Checking for one of the modules 'crayslepc_real;SLEPc'
-- Performing Test HAVE_NO_MULTLINE
-- Performing Test HAVE_NO_MULTLINE - Success
CMAKE_CXX_FLAGS = -Wno-comment  -D_GLIBCXX_USE_CXX11_ABI=1
CPM_DIR = ./Chipmunk-7.0.1/include
HEADERS = ./Chipmunk-7.0.1/include/chipmunk/chipmunk.h./Chipmunk-7.0.1/include/chipmunk/chipmunk_private.h
CMAKE_SOURCE_DIR = /home/fenics/shared/eQ
CMAKE_BINARY_DIR = /home/fenics/shared/eQ/build
GIT_BRANCH = master
GIT_COMMIT_HASH = 62bf999dee491a918981c0f8e380ccf184e1e6c7
-- Configuring done
-- Generating done
-- Build files have been written to: /home/fenics/shared/eQ/build
Scanning dependencies of target eQ
[  8%] Building CXX object CMakeFiles/eQ.dir/src/main.cpp.o
[ 16%] Building CXX object CMakeFiles/eQ.dir/src/fHSL.cpp.o
[ 25%] Building CXX object CMakeFiles/eQ.dir/src/inputOutput.cpp.o
[ 33%] Building CXX object CMakeFiles/eQ.dir/src/simulation.cpp.o
[ 41%] Building CXX object CMakeFiles/eQ.dir/src/Strain.cpp.o
[ 50%] Building CXX object CMakeFiles/eQ.dir/src/abm/eQabm.cpp.o
[ 58%] Building CXX object CMakeFiles/eQ.dir/src/abm/Ecoli.cpp.o
[ 66%] Building CXX object CMakeFiles/eQ.dir/src/abm/cpmHabitat.cpp.o
[ 75%] Building CXX object CMakeFiles/eQ.dir/src/abm/cpmTrap.cpp.o
[ 83%] Building CXX object CMakeFiles/eQ.dir/src/abm/cpmEcoli.cpp.o
[ 91%] Building CXX object CMakeFiles/eQ.dir/diffuclass.cpp.o
[100%] Linking CXX executable eQ
[100%] Built target eQ
MPI node #0 initially reads: time(nullptr) = 1628016433
MPI node #0 after MPI transfer = 1628016433
MPI node #0 reads: std::thread::hardware_concurrency() = 3

TESTING ACKNOWLEDGED
Testing, localArrayIndex = 0
root@e85404d46ca8:/home/fenics/shared/eQ# 
```



