#!/bin/bash
# change SLURM options below as needed before generating job submission scripts.
#SBATCH -p condo
#SBATCH -q condo
#SBATCH -A sio141
#SBATCH -N 1 -n 1 -c 1
#SBATCH --mem=8gb
#SBATCH -J make_phylodb
#SBATCH -t 1-00:00:00
#SBATCH -o grid_logs/makedb-%j.out-%N
#SBATCH -e grid_logs/makedb-%j.err-%N
#SBATCH --mail-type END,FAIL

### load fastp module and dependencies cleanly after purge
module purge
module load cpu slurm gcc
#export MODULEPATH=/projects/builder-group/jpg/modulefiles/applications:$MODULEPATH
export OMP_NUM_THREADS=$SLURM_JOB_CPUS_PER_NODE
export PATH=$PATH:/tscc/projects/ps-allenlab/projdata/common/bin

#set TMPDIR to scratch on execution host
export TMPDIR=/scratch/$USER/job_$SLURM_JOB_ID

### Change to Run directory
cd $SLURM_SUBMIT_DIR
mkdir -p grid_logs


# Path to the file
file_path="PhyloDB_2.0c"

DATETIME=`date '+Date: %F Time: %T'`;
echo "Time: make_phylodb STARTED for ${file_path}: ${DATETIME}";

cat database/*/*faa database/viruses.faa > $file_path.faa || exit 1
diamond makedb --in $file_path.faa --db $file_path

DATETIME=`date '+Date: %F Time: %T'`;
echo "Time: make_phylodb FINISHED for ${file_path}: ${DATETIME}";

# Loop through each line in the file
#while IFS= read -r line; do
#    # Pass the line to the Python script as an argument
#    python get_ranks_pr2.py "$line"
#done < "$file_path"