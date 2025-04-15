#!/bin/bash
#
# This script conducts a weak scaling study with agoge_run.
# For a fixed number of zones per process (baseline), the global grid dimensions
# are computed from the processor count by factoring np into Px, Py, and Pz.
#
# Usage: ./weak_scaling.sh [options]
#
# Options:
#   --mode <mode>            Communication mode: 'blocking' or 'nonblocking' (default: blocking)
#   --output-dir <dir>       Output directory (default: ./scaling_results)
#   --verify                 Compare results with single-proc reference solution
#   --baselines <list>       Comma-separated list of baseline resolutions (default: 16,32,64)
#   --procs <list>           Comma-separated list of processor counts (default: 1,4,16,64,128,256,512)
#   --executable <path>      Path to agoge_run executable (default: ./agoge_run)
#   --problem <path>         Path to problem file (default: ./problems/Sedov.yaml)
#   --help                   Show this help message

# Parse command line arguments
MODE="blocking"
OUTPUT_DIR="./scaling_results"
VERIFY=false
BASELINES="16,32,64"
PROCS="1,4,16,64,128,256,512"
EXECUTABLE="./agoge_run"
PROBLEM_FILE="./problems/Sedov.yaml"

while [[ $# -gt 0 ]]; do
    case $1 in
        --mode)
            MODE="$2"
            shift 2
            ;;
        --output-dir)
            OUTPUT_DIR="$2"
            shift 2
            ;;
        --verify)
            VERIFY=true
            shift
            ;;
        --baselines)
            BASELINES="$2"
            shift 2
            ;;
        --procs)
            PROCS="$2"
            shift 2
            ;;
        --executable)
            EXECUTABLE="$2"
            shift 2
            ;;
        --problem)
            PROBLEM_FILE="$2"
            shift 2
            ;;
        --help)
            echo "Usage: ./weak_scaling.sh [options]"
            echo ""
            echo "Options:"
            echo "  --mode <mode>            Communication mode: 'blocking' or 'nonblocking' (default: blocking)"
            echo "  --output-dir <dir>       Output directory (default: ./scaling_results)"
            echo "  --verify                 Compare results with single-proc reference solution"
            echo "  --baselines <list>       Comma-separated list of baseline resolutions (default: 16,32,64)"
            echo "  --procs <list>           Comma-separated list of processor counts (default: 1,4,16,64,128,256,512)"
            echo "  --executable <path>      Path to agoge_run executable (default: ./agoge_run)"
            echo "  --problem <path>         Path to problem file (default: ./problems/Sedov.yaml)"
            echo "  --help                   Show this help message"
            exit 0
            ;;
        *)
            echo "Unknown option: $1"
            exit 1
            ;;
    esac
done

# Set up directories
LOG_DIR="${OUTPUT_DIR}/logs"
DATA_DIR="${OUTPUT_DIR}/data"
REFERENCE_DIR="${OUTPUT_DIR}/reference"

# Create necessary directories
mkdir -p "${LOG_DIR}"
mkdir -p "${DATA_DIR}"
if [ "$VERIFY" = true ]; then
    mkdir -p "${REFERENCE_DIR}"
fi

# Export MODE as environment variable for agoge_run to detect
export AGOGE_COMM_MODE="${MODE}"

# Print configuration
echo "=== Weak Scaling Study Configuration ==="
echo "Communication Mode: ${MODE}"
echo "Output Directory: ${OUTPUT_DIR}"
echo "Verify with Reference: ${VERIFY}"
echo "Baseline Resolutions: ${BASELINES}"
echo "Process Counts: ${PROCS}"
echo "Executable: ${EXECUTABLE}"
echo "Problem File: ${PROBLEM_FILE}"
echo "======================================="

# Backup the original problem file
cp "${PROBLEM_FILE}" "${PROBLEM_FILE}.bak"

# Convert comma-separated lists to arrays
IFS=',' read -ra BASELINE_ARRAY <<< "$BASELINES"
IFS=',' read -ra PROC_ARRAY <<< "$PROCS"

# Run the scaling study
for baseline in "${BASELINE_ARRAY[@]}"; do
    echo "====== Baseline zones per process: ${baseline} ======"
    
    # Create a reference solution if verification is requested
    if [ "$VERIFY" = true ] && [ ! -d "${REFERENCE_DIR}/Sedov_${baseline}_single" ]; then
        echo "Creating reference solution for baseline ${baseline}"
        REF_DIR="${REFERENCE_DIR}/Sedov_${baseline}_single"
        mkdir -p "${REF_DIR}"
        
        # Set grid dimensions for single processor
        sed -i -E "s/^(nx:)[[:space:]]*[0-9]+/\1 ${baseline}/" "${PROBLEM_FILE}"
        sed -i -E "s/^(ny:)[[:space:]]*[0-9]+/\1 ${baseline}/" "${PROBLEM_FILE}"
        sed -i -E "s/^(nz:)[[:space:]]*[0-9]+/\1 ${baseline}/" "${PROBLEM_FILE}"
        
        # Save a copy of the modified problem file for reference
        cp "${PROBLEM_FILE}" "${LOG_DIR}/Sedov_${baseline}_single_${MODE}.yaml"
        
        # Run with a single processor for reference
        mpiexec -n 1 "${EXECUTABLE}" "${PROBLEM_FILE}" --output-dir="${REF_DIR}" > "${LOG_DIR}/reference_${baseline}_single.log" 2>&1
    fi
    
    for np in "${PROC_ARRAY[@]}"; do
        # Skip single processor case if it has already been used as reference
        if [ "$np" -eq 1 ] && [ "$VERIFY" = true ]; then
            echo "Skipping np=1 (already used as reference)"
            continue
        fi
        
        # Compute a factorization of np into (Px, Py, Pz) using a Python snippet.
        factors=$(python3 -c "import math, sys
np_val = int(sys.argv[1])
cube = int(round(np_val ** (1/3)))
Px = 1
for i in range(cube, 0, -1):
    if np_val % i == 0:
        Px = i
        break
rem = np_val // Px
sq = int(round(math.sqrt(rem)))
Py = 1
for i in range(sq, 0, -1):
    if rem % i == 0:
        Py = i
        break
Pz = rem // Py
print(Px, Py, Pz)" "$np")
        read Px Py Pz <<< "$factors"
        
        # Compute global grid dimensions so that each process receives exactly baseline zones
        NX=$(( baseline * Px ))
        NY=$(( baseline * Py ))
        NZ=$(( baseline * Pz ))
        
        # Define output locations
        LOG_FILE="${LOG_DIR}/weak_${MODE}_${baseline}_${np}procs.log"
        RUN_DIR="${DATA_DIR}/Sedov_${baseline}_${np}procs_${MODE}"
        
        echo "Running grid ${NX}x${NY}x${NZ} on ${np} procs. Output will be saved to ${LOG_FILE}"
        
        # Record run information in the log file
        {
            echo "-----------------------------------------------"
            echo "Weak Scaling Study: ${MODE} communication"
            echo "Global Grid: ${NX}x${NY}x${NZ}, Processors: ${np}"
            echo "Baseline zones/process: ${baseline}"
            echo "Decomposition: Px=${Px}, Py=${Py}, Pz=${Pz}"
            echo "Timestamp: $(date)"
            echo "-----------------------------------------------"
        } > "${LOG_FILE}"
        
        # Update the problem file with the new global grid dimensions
        sed -i -E "s/^(nx:)[[:space:]]*[0-9]+/\1 ${NX}/" "${PROBLEM_FILE}"
        sed -i -E "s/^(ny:)[[:space:]]*[0-9]+/\1 ${NY}/" "${PROBLEM_FILE}"
        sed -i -E "s/^(nz:)[[:space:]]*[0-9]+/\1 ${NZ}/" "${PROBLEM_FILE}"
        
        # Save a copy of the modified problem file
        cp "${PROBLEM_FILE}" "${LOG_DIR}/Sedov_${baseline}_${np}procs_${MODE}.yaml"
        
        # Create the output directory for this run
        mkdir -p "${RUN_DIR}"
        
        # Run agoge with the specified number of processors
        mpiexec -n "${np}" "${EXECUTABLE}" "${PROBLEM_FILE}" --output-dir="${RUN_DIR}" >> "${LOG_FILE}" 2>&1
        
        # Verify the result if requested and not the reference case
        if [ "$VERIFY" = true ]; then
            echo "Verifying ${np} proc result against reference solution..."
            python3 ../utils/agoge_compare.py "${REFERENCE_DIR}/Sedov_${baseline}_single" "${RUN_DIR}" \
                --output-dir="${LOG_DIR}/compare_${baseline}_${np}procs_${MODE}" \
                --threshold=1e-12 >> "${LOG_FILE}" 2>&1
        fi
    done
done

# Restore the original problem file
mv "${PROBLEM_FILE}.bak" "${PROBLEM_FILE}"

echo "Weak scaling study complete. Results are in ${OUTPUT_DIR}."
