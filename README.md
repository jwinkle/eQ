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

#### Overview

1. Install the Docker desktop for Mac, Windows, Linux:  https://www.docker.com/products/docker-desktop
2. Install the latest _Fenics_ image as found here: https://fenics.readthedocs.io/projects/containers/en/latest/introduction.html#running-fenics-in-docker
3. Run the ```eQinstall.sh``` script to fetch/compile/install dependencies and to compile eQ and test eQ within the Docker image for _Fenics_
4. Checkout an example branch (e.g. the latest PLoS article simulation code) to run and verify functionalith

#### Detailed Installation Instructions

