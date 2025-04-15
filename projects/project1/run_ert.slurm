#!/bin/bash --login
#SBATCH --job-name=ERT_Benchmark
#SBATCH --output=ERT_output_%j.txt
#SBATCH --error=ERT_error_%j.txt
#SBATCH --time=01:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=100G
#SBATCH --constraint=amd20  # Change this to test different architectures

module purge
module load gnuplot 

export REPO_DIR=~/cmse822/cmse822-codex-private
cd $SLURM_SUBMIT_DIR
cp $REPO_DIR/projects/project1/config.hpcc.msu.edu.01 ./config.local
# change the name of the results directory to match the architecture in the constraint
sed -i "s|ERT_RESULTS .*|ERT_RESULTS amd20_Results|" config.local
./ert config.local
