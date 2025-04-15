#!/bin/bash
#
# This script extracts performance metrics from scaling log files
# and writes the results to CSV files for later analysis.
#
# Usage: ./extract_zone_updates.sh [options]
#
# Options:
#   --input-dir <dir>        Directory containing log files (default: ./logs)
#   --output-dir <dir>       Directory for CSV output (default: same as input-dir)
#   --prefix <prefix>        Prefix for output CSV files (default: "")
#   --scaling <type>         Type of scaling to process: 'strong', 'weak', or 'both' (default: both)
#   --help                   Show this help message

# Parse command line arguments
INPUT_DIR="./logs"
OUTPUT_DIR=""
PREFIX=""
SCALING="both"

while [[ $# -gt 0 ]]; do
    case $1 in
        --input-dir)
            INPUT_DIR="$2"
            shift 2
            ;;
        --output-dir)
            OUTPUT_DIR="$2"
            shift 2
            ;;
        --prefix)
            PREFIX="$2"
            shift 2
            ;;
        --scaling)
            SCALING="$2"
            shift 2
            ;;
        --help)
            echo "Usage: ./extract_zone_updates.sh [options]"
            echo ""
            echo "Options:"
            echo "  --input-dir <dir>        Directory containing log files (default: ./logs)"
            echo "  --output-dir <dir>       Directory for CSV output (default: same as input-dir)"
            echo "  --prefix <prefix>        Prefix for output CSV files (default: \"\")"
            echo "  --scaling <type>         Type of scaling to process: 'strong', 'weak', or 'both' (default: both)"
            echo "  --help                   Show this help message"
            exit 0
            ;;
        *)
            echo "Unknown option: $1"
            exit 1
            ;;
    esac
done

# If no output directory specified, use the input directory
if [ -z "$OUTPUT_DIR" ]; then
    OUTPUT_DIR="$INPUT_DIR"
fi

# Create output directory if it doesn't exist
mkdir -p "$OUTPUT_DIR"

# Add prefix to output filenames if specified
if [ -n "$PREFIX" ]; then
    PREFIX="${PREFIX}_"
fi

# Process strong scaling logs if requested
if [ "$SCALING" = "strong" ] || [ "$SCALING" = "both" ]; then
    STRONG_OUTPUT="${OUTPUT_DIR}/${PREFIX}strong_scaling.csv"
    
    # CSV header
    echo "GridSize,Processors,ZonesPerProcess,ZoneUpdates(M),TimeStepSeconds,CommunicationTime" > "${STRONG_OUTPUT}"
    
    # Find and process strong scaling log files
    for log_file in "${INPUT_DIR}"/strong_*_*_*procs.log; do
        if [ ! -f "$log_file" ]; then
            echo "No strong scaling logs found in ${INPUT_DIR}"
            break
        fi
        
        # Extract metadata from filename
        # Format expected: strong_<mode>_<gridsize>_<num>procs.log
        base=$(basename "${log_file}")
        mode=$(echo "${base}" | cut -d'_' -f2)
        gridsize=$(echo "${base}" | cut -d'_' -f3)
        procs=$(echo "${base}" | cut -d'_' -f4 | sed 's/procs\.log//')
        
        # Calculate zones per process
        zones_per_proc=$(awk "/Baseline zones\/process:/ {print \$3}" "${log_file}")
        
        # If not explicitly listed in the file, calculate it
        if [ -z "$zones_per_proc" ]; then
            # Calculate as (gridsize^3)/procs if gridsize is available
            if [ -n "$gridsize" ]; then
                zones_per_proc=$(( (gridsize * gridsize * gridsize) / procs ))
            else
                zones_per_proc="N/A"
            fi
        fi
        
        # Extract performance metrics
        zone_updates=$(grep -E "Zone Updates per Second \(M\):" "${log_file}" | tail -1 | sed -E 's/.*Zone Updates per Second \(M\):\s*//')
        timestep_sec=$(grep -E "Time per timestep \(seconds\):" "${log_file}" | tail -1 | sed -E 's/.*Time per timestep \(seconds\):\s*//')
        comm_time=$(grep -E "Communication time \(seconds\):" "${log_file}" | tail -1 | sed -E 's/.*Communication time \(seconds\):\s*//')
        
        # Write CSV record if zone_updates was found
        if [ -n "$zone_updates" ]; then
            echo "${gridsize},${procs},${zones_per_proc},${zone_updates},${timestep_sec:-N/A},${comm_time:-N/A}" >> "${STRONG_OUTPUT}"
        fi
    done
    
    echo "Strong scaling CSV file created at: ${STRONG_OUTPUT}"
fi

# Process weak scaling logs if requested
if [ "$SCALING" = "weak" ] || [ "$SCALING" = "both" ]; then
    WEAK_OUTPUT="${OUTPUT_DIR}/${PREFIX}weak_scaling.csv"
    
    # CSV header
    echo "Baseline,Processors,GlobalGridSize,ZoneUpdates(M),TimeStepSeconds,CommunicationTime" > "${WEAK_OUTPUT}"
    
    # Find and process weak scaling log files
    for log_file in "${INPUT_DIR}"/weak_*_*_*procs.log; do
        if [ ! -f "$log_file" ]; then
            echo "No weak scaling logs found in ${INPUT_DIR}"
            break
        fi
        
        # Extract metadata from filename
        # Format expected: weak_<mode>_<baseline>_<num>procs.log
        base=$(basename "${log_file}")
        mode=$(echo "${base}" | cut -d'_' -f2)
        baseline=$(echo "${base}" | cut -d'_' -f3)
        procs=$(echo "${base}" | cut -d'_' -f4 | sed 's/procs\.log//')
        
        # Extract grid dimensions from the log file
        grid_size=$(grep -E "Global Grid:" "${log_file}" | sed -E 's/.*Global Grid: ([0-9]+)x([0-9]+)x([0-9]+).*/\1x\2x\3/')
        
        # Extract performance metrics
        zone_updates=$(grep -E "Zone Updates per Second \(M\):" "${log_file}" | tail -1 | sed -E 's/.*Zone Updates per Second \(M\):\s*//')
        timestep_sec=$(grep -E "Time per timestep \(seconds\):" "${log_file}" | tail -1 | sed -E 's/.*Time per timestep \(seconds\):\s*//')
        comm_time=$(grep -E "Communication time \(seconds\):" "${log_file}" | tail -1 | sed -E 's/.*Communication time \(seconds\):\s*//')
        
        # Write CSV record if zone_updates was found
        if [ -n "$zone_updates" ]; then
            echo "${baseline},${procs},${grid_size:-N/A},${zone_updates},${timestep_sec:-N/A},${comm_time:-N/A}" >> "${WEAK_OUTPUT}"
        fi
    done
    
    echo "Weak scaling CSV file created at: ${WEAK_OUTPUT}"
fi

# If both types were processed, generate a combined comparative CSV
if [ "$SCALING" = "both" ]; then
    COMBINED_OUTPUT="${OUTPUT_DIR}/${PREFIX}combined_scaling.csv"
    
    # CSV header
    echo "ScalingType,Resolution,Processors,ZonesPerProcess,ZoneUpdates(M),TimeStepSeconds,CommunicationTime" > "${COMBINED_OUTPUT}"
    
    # Add strong scaling data
    if [ -f "${STRONG_OUTPUT}" ]; then
        # Skip header
        tail -n +2 "${STRONG_OUTPUT}" | while IFS=, read -r gridsize procs zones_per_proc updates timestep commtime; do
            echo "Strong,${gridsize},${procs},${zones_per_proc},${updates},${timestep},${commtime}" >> "${COMBINED_OUTPUT}"
        done
    fi
    
    # Add weak scaling data
    if [ -f "${WEAK_OUTPUT}" ]; then
        # Skip header
        tail -n +2 "${WEAK_OUTPUT}" | while IFS=, read -r baseline procs grid_size updates timestep commtime; do
            echo "Weak,${baseline},${procs},${baseline},${updates},${timestep},${commtime}" >> "${COMBINED_OUTPUT}"
        done
    fi
    
    echo "Combined scaling CSV file created at: ${COMBINED_OUTPUT}"
    
    # Generate efficiency metrics for both scaling types
    EFFICIENCY_OUTPUT="${OUTPUT_DIR}/${PREFIX}scaling_efficiency.csv"
    echo "ScalingType,Resolution,Processors,ParallelEfficiency,CommunicationOverhead" > "${EFFICIENCY_OUTPUT}"
    
    # Process strong scaling efficiency (relative to single processor case)
    if [ -f "${STRONG_OUTPUT}" ]; then
        # Group by resolution and find the baseline performance for each resolution
        awk -F, 'NR > 1 {if (!($1 in baseline_procs) || $2 < baseline_procs[$1]) {
            baseline_procs[$1] = $2;
            baseline_perf[$1] = $4 * $2; # zone_updates * procs for comparison
        }}
        END {
            for (res in baseline_procs) {
                print res "," baseline_procs[res] "," baseline_perf[res];
            }
        }' "${STRONG_OUTPUT}" > "${OUTPUT_DIR}/temp_baselines.csv"
        
        # Calculate efficiency for each data point
        awk -F, 'BEGIN {
            # Load baseline performance data
            while (getline < "'${OUTPUT_DIR}'/temp_baselines.csv") {
                split($0, parts, ",");
                base_perf[parts[1]] = parts[3]; # res -> perf mapping
            }
        }
        NR > 1 {
            # Skip header and calculate efficiency
            res = $1;
            procs = $2;
            perf = $4 * procs; # zone_updates * procs
            if (res in base_perf) {
                efficiency = perf / base_perf[res];
                comm_overhead = $6 / $5; # comm_time / timestep_time
                print "Strong," res "," procs "," efficiency "," comm_overhead;
            }
        }' "${STRONG_OUTPUT}" >> "${EFFICIENCY_OUTPUT}"
    fi
    
    # Process weak scaling efficiency (should remain close to 1.0)
    if [ -f "${WEAK_OUTPUT}" ]; then
        # Group by baseline and find the baseline performance for each resolution
        awk -F, 'NR > 1 {if (!($1 in baseline_procs) || $2 < baseline_procs[$1]) {
            baseline_procs[$1] = $2;
            baseline_perf[$1] = $4; # zone_updates for comparison
        }}
        END {
            for (res in baseline_procs) {
                print res "," baseline_procs[res] "," baseline_perf[res];
            }
        }' "${WEAK_OUTPUT}" > "${OUTPUT_DIR}/temp_baselines_weak.csv"
        
        # Calculate efficiency for each data point
        awk -F, 'BEGIN {
            # Load baseline performance data
            while (getline < "'${OUTPUT_DIR}'/temp_baselines_weak.csv") {
                split($0, parts, ",");
                base_perf[parts[1]] = parts[3]; # baseline -> perf mapping
            }
        }
        NR > 1 {
            # Skip header and calculate efficiency
            baseline = $1;
            procs = $2;
            perf = $4; # zone_updates
            if (baseline in base_perf) {
                efficiency = perf / base_perf[baseline];
                comm_overhead = $6 / $5; # comm_time / timestep_time
                print "Weak," baseline "," procs "," efficiency "," comm_overhead;
            }
        }' "${WEAK_OUTPUT}" >> "${EFFICIENCY_OUTPUT}"
    fi
    
    # Clean up temporary files
    rm -f "${OUTPUT_DIR}/temp_baselines.csv"
    rm -f "${OUTPUT_DIR}/temp_baselines_weak.csv"
    
    echo "Scaling efficiency metrics saved to: ${EFFICIENCY_OUTPUT}"
fi

echo "Data extraction complete. Use these CSV files with plotting tools to visualize performance."
