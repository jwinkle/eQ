# How to Run eQ code for the ABM oscillator

RUNNING OSCILLATOR FROM PLOS COMPUTATIONAL BIOLOGY PAPER (2021)

1. setup the eQ environment from the **master** branch of eQ (see `README`)
2. from the PLosCompBio21 branch (where this file is located), run the following 3 scripts:
  - [approx. 3 minutes] `eQcompileBranch.sh`
  - [approx. 3 hours] `eQrunBranchN3.sh`
  - [approx. 30 minutes] `eQdecodeDataMatlab.sh`
3. summary:
  - the first re-builds the Fenics header files, runs cmake, and then compiles the eQ executable by running make.  
  - the second runs the ABM oscillator using 3 MPI nodes; generates data files for the ABM oscillator from PLoSCompBio (2021; ref here)
  - the third script decodes the data using Matlab (requires matlab executable in path)
4. data is generated in the folder `/build/images/0/` image files are generated and can be visualized using ImageJ, QuickTime, or other image-sequencer playback software.
