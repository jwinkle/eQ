

# USAGE:  from the directory above eQ, run (for example to run "--array=0-3") :
# $>  ./jlaunch.sh 0 3


TIME_SINCE_EPOCH=$(date +%s)
DATE_TIME_STRING=$(date +%F)_$(date +%H)-$(date +%M)-$(date +%S)

BASE_DIR="/project/josic/winkle/job/"
JOB_DIR="${BASE_DIR}${TIME_SINCE_EPOCH}_${DATE_TIME_STRING}"

./jbuild.sh

mkdir -p ${JOB_DIR}
cp ~/eQ/build/eQ ${JOB_DIR}
cp ./winkle.job ${JOB_DIR}/


# echo "exporting: modulusOption=$3"
# export modulusOption="$3"
# echo "option is: ${modulusOption}"

echo "running:  sbatch --array=$1-$2 ${JOB_DIR}/winkle.job"

cd ${JOB_DIR}

sbatch --array=$1-$2 ./winkle.job $3

echo "sbatch job launched...reading squeue -u jjwinkle..."

sleep 2
squeue -u jjwinkle -v

echo "job directory contents:  ls -l ${BASE_DIR}"
ls -l ${BASE_DIR}
