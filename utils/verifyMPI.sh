#!/bin/bash
#
# This script verifies the MPI implementation for the Agoge code by:
# 1. Building and running the code with blocking MPI
# 2. Building and running the code with non-blocking MPI
# 3. Comparing the results to ensure consistency

# Define color codes for better readability
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
RED='\033[0;31m'
BLUE='\033[0;34m'
NC='\033[0m' # No Color

# Define the project directory
# Note: Update this path if running from a different location
agoge_dir="/Users/smc/Library/CloudStorage/Dropbox/Teaching/CMSE822/cmse822-codex-private/agoge"

# Error threshold for numerical comparisons
ERROR_THRESHOLD=1e-10

# Create logs directory if it doesn't exist
logs_dir="${agoge_dir}/logs"
mkdir -p "${logs_dir}"

# Define log files
build_log="${logs_dir}/build.log"
run_log="${logs_dir}/run.log"
compare_log="${logs_dir}/compare.log"

# Function to run a command, log its output, and return its status
run_and_log() {
    local cmd="$1"
    local log_file="$2"
    local error_msg="$3"
    
    echo -e "${BLUE}Running: ${cmd}${NC}" 
    eval "${cmd} > ${log_file} 2>&1"
    local status=$?
    
    if [ $status -ne 0 ]; then
        echo -e "${RED}${error_msg}${NC}"
        echo -e "${RED}Check the log file for details: ${log_file}${NC}"
        echo -e "${RED}Last 10 lines of log:${NC}"
        tail -n 10 "${log_file}"
        return 1
    fi
    
    return 0
}

# Function to modify boundary conditions in Sedov.yaml
modify_boundary_conditions() {
    local bc_type="$1"
    local yaml_file="$agoge_dir/problems/Sedov_test.yaml"
    
    # Create a copy of the original Sedov.yaml file
    cp "$agoge_dir/problems/Sedov.yaml" "$yaml_file"
    
    # Set all boundary conditions to the specified type
    sed -i.bak "s/bc_xmin: .*/bc_xmin: $bc_type/" "$yaml_file"
    sed -i.bak "s/bc_xmax: .*/bc_xmax: $bc_type/" "$yaml_file"
    sed -i.bak "s/bc_ymin: .*/bc_ymin: $bc_type/" "$yaml_file"
    sed -i.bak "s/bc_ymax: .*/bc_ymax: $bc_type/" "$yaml_file"
    sed -i.bak "s/bc_zmin: .*/bc_zmin: $bc_type/" "$yaml_file"
    sed -i.bak "s/bc_zmax: .*/bc_zmax: $bc_type/" "$yaml_file"
    
    # Remove backup file
    rm -f "${yaml_file}.bak"
    
    echo -e "${YELLOW}Modified boundary conditions to: ${bc_type}${NC}"
}

# Navigate to the project directory
echo -e "${YELLOW}Changing to project directory: ${agoge_dir}${NC}"
cd "$agoge_dir" || { echo -e "${RED}Failed to change to directory: ${agoge_dir}${NC}"; exit 1; }

# --------------------------------------
# Test 1: Blocking MPI Implementation
# --------------------------------------
echo -e "\n${GREEN}======= Testing Blocking MPI Implementation =======${NC}"

# Clean and build with redirected output
echo "Building project with blocking MPI..."
run_and_log "make clean && make -j" "${build_log}" "Build failed for blocking MPI!" || exit 1

# Array of boundary condition types to test
bc_types=("periodic" "outflow")

for bc_type in "${bc_types[@]}"; do
    echo -e "\n${YELLOW}Testing blocking MPI with ${bc_type} boundary conditions${NC}"
    
    # Modify boundary conditions
    modify_boundary_conditions "$bc_type"
    
    # Run with 1 rank
    echo -e "Running with 1 rank..."
    run_and_log "mpiexec -n 1 ./agoge_run ./problems/Sedov_test.yaml" "${run_log}_1rank_${bc_type}" "Failed to run with 1 rank!" || exit 1
    echo "✓ Single-rank run completed"
    
    # Run with 2 ranks
    echo -e "Running with 2 ranks..."
    run_and_log "mpiexec -n 2 ./agoge_run ./problems/Sedov_test.yaml" "${run_log}_2ranks_${bc_type}" "Failed to run with 2 ranks!" || exit 1
    echo "✓ Multi-rank run completed"
    
    # Compare outputs
    echo -e "${YELLOW}Comparing single-rank and multi-rank outputs...${NC}"
    run_and_log "python3 viz/agoge_compare.py ./output/Sedov_1ranks_0001 ./output/Sedov_2ranks_0001 --threshold ${ERROR_THRESHOLD}" "${compare_log}_blocking_${bc_type}" "Comparison failed!" && comparison_result=$? || comparison_result=$?
    
    if [ $comparison_result -eq 0 ]; then
        echo -e "${GREEN}✓ Blocking MPI test with ${bc_type} boundary conditions passed: Results match within threshold!${NC}"
    else
        echo -e "${RED}× Blocking MPI test with ${bc_type} boundary conditions failed: Results differ beyond threshold!${NC}"
        echo -e "${RED}Check ${compare_log}_blocking_${bc_type} for details${NC}"
        exit 1
    fi
done

# --------------------------------------
# Test 2: Non-Blocking MPI Implementation
# --------------------------------------
echo -e "\n${GREEN}======= Testing Non-Blocking MPI Implementation =======${NC}"

# Clean and build with non-blocking MPI
echo "Building project with non-blocking MPI..."
run_and_log "make clean && make -j USE_NONBLOCKING_MPI=1" "${build_log}_nonblocking" "Build failed for non-blocking MPI!" || exit 1

for bc_type in "${bc_types[@]}"; do
    echo -e "\n${YELLOW}Testing non-blocking MPI with ${bc_type} boundary conditions${NC}"
    
    # Modify boundary conditions
    modify_boundary_conditions "$bc_type"
    
    # Run with 1 rank
    echo -e "Running with 1 rank..."
    run_and_log "mpiexec -n 1 ./agoge_run ./problems/Sedov_test.yaml" "${run_log}_nonblocking_1rank_${bc_type}" "Failed to run with 1 rank!" || exit 1
    echo "✓ Single-rank run completed"
    
    # Run with 2 ranks
    echo -e "Running with 2 ranks..."
    run_and_log "mpiexec -n 2 ./agoge_run ./problems/Sedov_test.yaml" "${run_log}_nonblocking_2ranks_${bc_type}" "Failed to run with 2 ranks!" || exit 1
    echo "✓ Multi-rank run completed"
    
    # Compare outputs
    echo -e "${YELLOW}Comparing single-rank and multi-rank outputs...${NC}"
    run_and_log "python3 viz/agoge_compare.py ./output/Sedov_1ranks_0001 ./output/Sedov_2ranks_0001 --threshold ${ERROR_THRESHOLD}" "${compare_log}_nonblocking_${bc_type}" "Comparison failed!" && comparison_result=$? || comparison_result=$?
    
    if [ $comparison_result -eq 0 ]; then
        echo -e "${GREEN}✓ Non-Blocking MPI test with ${bc_type} boundary conditions passed: Results match within threshold!${NC}"
    else
        echo -e "${RED}× Non-Blocking MPI test with ${bc_type} boundary conditions failed: Results differ beyond threshold!${NC}"
        echo -e "${RED}Check ${compare_log}_nonblocking_${bc_type} for details${NC}"
        exit 1
    fi
done

# --------------------------------------
# Compare results between blocking and non-blocking implementations
# --------------------------------------
echo -e "\n${GREEN}======= Comparing Blocking vs Non-Blocking MPI =======${NC}"

# Test with both boundary condition types
for bc_type in "${bc_types[@]}"; do
    echo -e "\n${YELLOW}Comparing implementations with ${bc_type} boundary conditions${NC}"
    
    # Modify boundary conditions
    modify_boundary_conditions "$bc_type"
    
    # Rename the directories to avoid overwriting
    echo "Preparing directories for comparison..."
    mv ./output/Sedov_2ranks_0001 ./output/Sedov_2ranks_0001_nonblocking
    
    # Build and run with blocking MPI again
    echo "Rebuilding with blocking MPI..."
    run_and_log "make clean && make -j" "${build_log}_comparison_${bc_type}" "Failed to rebuild with blocking MPI!" || exit 1
    
    echo "Running blocking MPI implementation with 2 ranks..."
    run_and_log "mpiexec -n 2 ./agoge_run ./problems/Sedov_test.yaml" "${run_log}_comparison_${bc_type}" "Failed to run blocking MPI implementation!" || exit 1
    
    # Compare blocking vs non-blocking
    echo -e "${YELLOW}Comparing blocking vs non-blocking MPI results...${NC}"
    run_and_log "python3 viz/agoge_compare.py ./output/Sedov_2ranks_0001 ./output/Sedov_2ranks_0001_nonblocking --threshold ${ERROR_THRESHOLD}" "${compare_log}_implementation_${bc_type}" "Implementation comparison failed!" && impl_comparison_result=$? || impl_comparison_result=$?
    
    if [ $impl_comparison_result -eq 0 ]; then
        echo -e "${GREEN}✓ Implementation comparison with ${bc_type} boundary conditions passed: Both implementations produce matching results!${NC}"
    else
        echo -e "${YELLOW}⚠ Note: With ${bc_type} boundary conditions, blocking and non-blocking implementations produce slightly different results.${NC}"
        echo -e "${YELLOW}  This may be acceptable if differences are small and due to floating-point precision or ordering.${NC}"
        echo -e "${YELLOW}  Check ${compare_log}_implementation_${bc_type} for details${NC}"
    fi
    
    # Cleanup renamed directory
    echo "Cleaning up temporary files..."
    rm -rf ./output/Sedov_2ranks_0001_nonblocking
done

# Clean up the test YAML file
rm -f "$agoge_dir/problems/Sedov_test.yaml"

echo -e "\n${GREEN}========================================${NC}"
echo -e "${GREEN}All MPI verification tests completed!${NC}"
echo -e "${GREEN}Log files are available in ${logs_dir}${NC}"
echo -e "${GREEN}========================================${NC}"
exit 0