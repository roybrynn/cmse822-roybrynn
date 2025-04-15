#!/bin/bash --login
#SBATCH --job-name=scaling_study
#SBATCH --nodes=4
#SBATCH --ntasks=512
#SBATCH --mem=493G
#SBATCH --time=3:00:00
#SBATCH --constraint=amd20

# Communication mode (blocking, nonblocking, or both)
# Default: both
COMM_MODE=${1:-"both"}

# Load necessary modules
module purge
module load HDF5

# Change to the directory containing the scripts
cd /mnt/home/scouch/cmse822/cmse822-codex-private/agoge

# Create timestamped directory for this run
TIMESTAMP=$(date +%Y%m%d_%H%M%S)
RESULTS_DIR="../results/scaling_study_${TIMESTAMP}"
mkdir -p ${RESULTS_DIR}

# Save job information
echo "Job ID: ${SLURM_JOB_ID}" > ${RESULTS_DIR}/job_info.txt
echo "Nodes: ${SLURM_JOB_NUM_NODES}" >> ${RESULTS_DIR}/job_info.txt
echo "Tasks: ${SLURM_NTASKS}" >> ${RESULTS_DIR}/job_info.txt
echo "Communication mode: ${COMM_MODE}" >> ${RESULTS_DIR}/job_info.txt
echo "Start time: $(date)" >> ${RESULTS_DIR}/job_info.txt

# Environment setup (for tracing path to the current executable)
echo "Current agoge_run path: $(which agoge_run)" >> ${RESULTS_DIR}/job_info.txt
agoge_run --version >> ${RESULTS_DIR}/job_info.txt 2>&1 || echo "agoge_run version not available" >> ${RESULTS_DIR}/job_info.txt

# Run scaling studies based on communication mode
if [[ "${COMM_MODE}" == "blocking" || "${COMM_MODE}" == "both" ]]; then
    echo "Running blocking communication scaling studies..."
    
    # Run the weak scaling study for blocking
    bash ../utils/weak_scaling.sh --mode blocking --output-dir ${RESULTS_DIR}/blocking
    
    # Run the strong scaling study for blocking
    bash ../utils/strong_scaling.sh --mode blocking --output-dir ${RESULTS_DIR}/blocking
    
    # Extract and plot the results
    bash ../utils/extract_zone_updates.sh --input-dir ${RESULTS_DIR}/blocking/logs --prefix blocking
    python3 ../utils/plot_performance.py --input-dir ${RESULTS_DIR}/blocking --prefix blocking
fi

if [[ "${COMM_MODE}" == "nonblocking" || "${COMM_MODE}" == "both" ]]; then
    echo "Running non-blocking communication scaling studies..."
    
    # Run the weak scaling study for non-blocking
    bash ../utils/weak_scaling.sh --mode nonblocking --output-dir ${RESULTS_DIR}/nonblocking
    
    # Run the strong scaling study for non-blocking
    bash ../utils/strong_scaling.sh --mode nonblocking --output-dir ${RESULTS_DIR}/nonblocking
    
    # Extract and plot the results
    bash ../utils/extract_zone_updates.sh --input-dir ${RESULTS_DIR}/nonblocking/logs --prefix nonblocking
    python3 ../utils/plot_performance.py --input-dir ${RESULTS_DIR}/nonblocking --prefix nonblocking
fi

# Generate comparative plots if both modes were run
if [[ "${COMM_MODE}" == "both" ]]; then
    echo "Generating comparative plots..."
    python3 ../utils/plot_performance.py --comparative --input-dir ${RESULTS_DIR}
fi

# Record completion time
echo "End time: $(date)" >> ${RESULTS_DIR}/job_info.txt
echo "Results available at: ${RESULTS_DIR}"
