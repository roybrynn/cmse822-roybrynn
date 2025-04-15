#!/bin/bash
#
# This script conducts a strong scaling study with agoge_run.
# It updates problems/Sedov.yaml for various grid resolutions
# and runs mpiexec jobs for different process counts.
#
# Usage: ./strong_scaling.sh [options]
#
# Options:
#   --mode <mode>            Communication mode: 'blocking' or 'nonblocking' (default: blocking)
#   --output-dir <dir>       Output directory (default: ./scaling_results)
#   --verify                 Compare results with single-proc reference solution
#   --resolutions <list>     Comma-separated list of grid resolutions (default: 64,128,256,512)
#   --procs <list>           Comma-separated list of processor counts (default: 1,2,4,8,16,32,64,128,256,512)
#   --executable <path>      Path to agoge_run executable (default: ./agoge_run)
#   --problem <path>         Path to problem file (default: ./problems/Sedov.yaml)
#   --help                   Show this help message

# Parse command line arguments
MODE="blocking"
OUTPUT_DIR="./scaling_results"
VERIFY=false
RESOLUTIONS="64,128,256,512"
PROCS="1,2,4,8,16,32,64,128,256,512"
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
        --resolutions)
            RESOLUTIONS="$2"
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
            echo "Usage: ./strong_scaling.sh [options]"
            echo ""
            echo "Options:"
            echo "  --mode <mode>            Communication mode: 'blocking' or 'nonblocking' (default: blocking)"
            echo "  --output-dir <dir>       Output directory (default: ./scaling_results)"
            echo "  --verify                 Compare results with single-proc reference solution"
            echo "  --resolutions <list>     Comma-separated list of grid resolutions (default: 64,128,256,512)"
            echo "  --procs <list>           Comma-separated list of processor counts (default: 1,2,4,8,16,32,64,128,256,512)"
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
echo "=== Strong Scaling Study Configuration ==="
echo "Communication Mode: ${MODE}"
echo "Output Directory: ${OUTPUT_DIR}"
echo "Verify with Reference: ${VERIFY}"
echo "Grid Resolutions: ${RESOLUTIONS}"
echo "Process Counts: ${PROCS}"
echo "Executable: ${EXECUTABLE}"
echo "Problem File: ${PROBLEM_FILE}"
echo "======================================="

# Backup the original problem file
cp "${PROBLEM_FILE}" "${PROBLEM_FILE}.bak"

# Convert comma-separated lists to arrays
IFS=',' read -ra RES_ARRAY <<< "$RESOLUTIONS"
IFS=',' read -ra PROC_ARRAY <<< "$PROCS"

# Run the scaling study
for res in "${RES_ARRAY[@]}"; do
    echo "====== Grid Resolution: ${res}x${res}x${res} ======"
    
    # Update the problem file with the current resolution
    sed -i -E "s/^(nx:)[[:space:]]*[0-9]+/\1 ${res}/" "${PROBLEM_FILE}"
    sed -i -E "s/^(ny:)[[:space:]]*[0-9]+/\1 ${res}/" "${PROBLEM_FILE}"
    sed -i -E "s/^(nz:)[[:space:]]*[0-9]+/\1 ${res}/" "${PROBLEM_FILE}"
    
    # Save a copy of the modified problem file for reference
    cp "${PROBLEM_FILE}" "${LOG_DIR}/Sedov_${res}_${MODE}.yaml"
    
    # Create a reference solution if verification is requested
    if [ "$VERIFY" = true ] && [ ! -d "${REFERENCE_DIR}/Sedov_${res}_single" ]; then
        echo "Creating reference solution for ${res}x${res}x${res}"
        REF_DIR="${REFERENCE_DIR}/Sedov_${res}_single"
        mkdir -p "${REF_DIR}"
        
        # Run with a single processor for reference
        mpiexec -n 1 "${EXECUTABLE}" "${PROBLEM_FILE}" --output-dir="${REF_DIR}" > "${LOG_DIR}/reference_${res}_single.log" 2>&1
    fi
    
    for np in "${PROC_ARRAY[@]}"; do
        # Skip single processor case if it has already been used as reference
        if [ "$np" -eq 1 ] && [ "$VERIFY" = true ]; then
            echo "Skipping np=1 (already used as reference)"
            continue
        fi
        
        # Define output locations
        LOG_FILE="${LOG_DIR}/strong_${MODE}_${res}_${np}procs.log"
        RUN_DIR="${DATA_DIR}/Sedov_${res}_${np}procs_${MODE}"
        
        echo "Running grid ${res}x${res}x${res} with ${np} procs. Output will be saved to ${LOG_FILE}"
        
        # Record run information in the log file
        {
            echo "-----------------------------------------------"
            echo "Strong Scaling Study: ${MODE} communication"
            echo "Grid: ${res}x${res}x${res}, Processors: ${np}"
            echo "Timestamp: $(date)"
            echo "Baseline zones/process: $((res*res*res/np))"
            echo "-----------------------------------------------"
        } > "${LOG_FILE}"
        
        # Create the output directory for this run
        mkdir -p "${RUN_DIR}"
        
        # Run agoge with the specified number of processors
        mpiexec -n "${np}" "${EXECUTABLE}" "${PROBLEM_FILE}" --output-dir="${RUN_DIR}" >> "${LOG_FILE}" 2>&1
        
        # Verify the result if requested and not the reference case
        if [ "$VERIFY" = true ]; then
            echo "Verifying ${np} proc result against reference solution..."
            python3 ../utils/agoge_compare.py "${REFERENCE_DIR}/Sedov_${res}_single" "${RUN_DIR}" \
                --output-dir="${LOG_DIR}/compare_${res}_${np}procs_${MODE}" \
                --threshold=1e-12 >> "${LOG_FILE}" 2>&1
        fi
    done
done

# Restore the original problem file
mv "${PROBLEM_FILE}.bak" "${PROBLEM_FILE}"

echo "Strong scaling study complete. Results are in ${OUTPUT_DIR}."
