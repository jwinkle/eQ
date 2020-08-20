

# USAGE:  from the directory above eQ, run (for example to run "--array=0-3") :
# $>  ./jlaunch.sh 0 3



JOB_DIR="/project/josic/winkle/job"

./jbuild.sh

mkdir -p ~/eQ/run
cp ~/eQ/build/eQ ~/eQ/run/
cp ./winkle.job ${JOB_DIR}/


# echo "exporting: modulusOption=$3"
# export modulusOption="$3"
# echo "option is: ${modulusOption}"

echo "running:  sbatch --array=$1-$2 ${JOB_DIR}/winkle.job"

cd ${JOB_DIR}

sbatch --array=$1-$2 ./winkle.job

echo "sbatch job launched...reading squeue -u jjwinkle..."

sleep 2
squeue -u jjwinkle -v

