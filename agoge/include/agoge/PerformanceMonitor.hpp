#pragma once

#include <mpi.h>

#include <chrono>
#include <iomanip>
#include <iostream>
#include <stack>
#include <string>
#include <unordered_map>
#include <vector>

/**
 * @file PerformanceMonitor.hpp
 * @brief Singleton class for timing various parts of the Agoge solver.
 */

namespace agoge {

/**
 * @struct TimerData
 * @brief Stores total accumulated time and call count for a particular timer.
 */
struct TimerData {
    std::chrono::steady_clock::duration total{};
    long callCount = 0;
    int nestingLevel = 0;  // Track nesting level for indentation
    std::vector<std::string>
        children{};  // Child timers started within this timer
};

/**
 * @class PerformanceMonitor
 * @brief A singleton class that manages named timers for performance tracking.
 *        Supports nested timer hierarchies.
 *
 * Usage:
 *   PerformanceMonitor::instance().startTimer("computeL");
 *   // ... code ...
 *   PerformanceMonitor::instance().startTimer("computeSubStep"); // Nested
 * timer
 *   // ... nested code ...
 *   PerformanceMonitor::instance().stopTimer("computeSubStep");
 *   // ... more code ...
 *   PerformanceMonitor::instance().stopTimer("computeL");
 *
 *   At program end:
 *   PerformanceMonitor::instance().setSteps(steps);
 *   PerformanceMonitor::instance().setZones(zones);
 *   PerformanceMonitor::instance().compileReport();
 */
class PerformanceMonitor {
   public:
    /**
     * @brief Access the global singleton instance of PerformanceMonitor.
     */
    static PerformanceMonitor& instance() {
        static PerformanceMonitor s_instance;
        return s_instance;
    }

    /**
     * @brief Start a timer for the given name. If already running, it restarts.
     *        If another timer is active, this one becomes a child of it.
     *
     * @param name Unique string identifying this timer (e.g., "computeL").
     */
    void startTimer(const std::string& name) {
        auto now = std::chrono::steady_clock::now();
        startTimes_[name] = now;

        // Record the timer hierarchy
        const int nestingLevel =
            !_activeTimerStack.empty()
                ? timers_[_activeTimerStack.top()].nestingLevel + 1
                : 0;

        // Initialize the timer if it's the first time we're seeing it
        if (timers_.find(name) == timers_.end()) {
            timers_[name].nestingLevel = nestingLevel;
        }

        // Record parent-child relationship
        if (!_activeTimerStack.empty()) {
            const std::string& parent = _activeTimerStack.top();
            // Add this timer as child to parent (if not already there)
            auto& parentChildren = timers_[parent].children;
            if (std::find(parentChildren.begin(), parentChildren.end(), name) ==
                parentChildren.end()) {
                parentChildren.push_back(name);
            }
        }

        // Push this timer onto active stack
        _activeTimerStack.push(name);
    }

    /**
     * @brief Stop a timer for the given name, accumulate the elapsed time.
     *
     * @param name Must match a previously started timer.
     */
    void stopTimer(const std::string& name) {
        auto now = std::chrono::steady_clock::now();
        auto it = startTimes_.find(name);
        if (it != startTimes_.end()) {
            auto elapsed = now - it->second;
            timers_[name].total += elapsed;
            timers_[name].callCount += 1;

            // Check if we're stopping the timer at the top of the stack
            // This helps handle cases where timers are stopped out of order
            if (!_activeTimerStack.empty() && _activeTimerStack.top() == name) {
                _activeTimerStack.pop();
            } else if (!_activeTimerStack.empty()) {
                // This is unusual - stopping a timer that's not at the top of
                // the stack
                std::cerr << "[PerformanceMonitor] Warning: Timer '" << name
                          << "' stopped out of order.\n";
            }
        } else {
            std::cerr << "[PerformanceMonitor] Warning: Timer '" << name
                      << "' was not started.\n";
        }
    }

    /**
     * @brief Set the total number of simulation steps.
     *
     * This value is used to compute zone updates per second.
     *
     * @param steps Total number of simulation steps taken.
     */
    void setSteps(long steps) { steps_ = steps; }

    /**
     * @brief Set the total number of zones in the simulation domain.
     *
     * This value is used to compute zone updates per second.
     *
     * @param zones Total number of zones in the domain.
     */
    void setZones(long zones) { zones_ = zones; }

    /**
     * @brief Set the processor rank.
     *
     * @param rank Processor rank.
     */
    void setRank(int rank) { rank_ = rank; }

    /**
     * @brief Set the total number of ranks.
     *
     * @param size Total number of ranks.
     */
    void setCommSize(int size) { size_ = size; }

    /**
     * @brief Compile and print a summary of all timers.
     *
     * All MPI collectives within this method are called by every rank.
     * However, the detailed report is printed only by rank 0.
     */
    void compileReport() const {
        // Gather local timer names and data. All ranks participate.
        std::vector<std::string> timerNames;
        for (const auto& pair : timers_) {
            timerNames.push_back(pair.first);
        }
        int count = timerNames.size();
        std::vector<double> localTotals(count);
        std::vector<double> localCalls(count);
        int idx = 0;
        for (const auto& pair : timers_) {
            localTotals[idx] =
                std::chrono::duration<double>(pair.second.total).count();
            localCalls[idx] = static_cast<double>(pair.second.callCount);
            idx++;
        }

        std::vector<double> minTotals(count);
        std::vector<double> maxTotals(count);
        std::vector<double> sumTotals(count);

        MPI_Allreduce(localTotals.data(), minTotals.data(), count, MPI_DOUBLE,
                      MPI_MIN, MPI_COMM_WORLD);
        MPI_Allreduce(localTotals.data(), maxTotals.data(), count, MPI_DOUBLE,
                      MPI_MAX, MPI_COMM_WORLD);
        MPI_Allreduce(localTotals.data(), sumTotals.data(), count, MPI_DOUBLE,
                      MPI_SUM, MPI_COMM_WORLD);

        if (rank_ == 0) {
            printHierarchicalReport(timerNames, minTotals, maxTotals, sumTotals,
                                    localCalls);
        }
    }

   private:
    PerformanceMonitor() = default;  // Private constructor for singleton
    PerformanceMonitor(const PerformanceMonitor&) = delete;
    PerformanceMonitor& operator=(const PerformanceMonitor&) = delete;

    // Print report showing the timer hierarchy
    void printHierarchicalReport(const std::vector<std::string>& timerNames,
                                 const std::vector<double>& minTotals,
                                 const std::vector<double>& maxTotals,
                                 const std::vector<double>& sumTotals,
                                 const std::vector<double>& localCalls) const {
        const int headerWidth = 80;
        const int nameWidth = 28;
        const int colWidth = 12;

        std::cout
            << "\n====== HIERARCHICAL TIMINGS ACROSS ALL RANKS ======\n\n";
        std::cout << std::left << std::setw(nameWidth) << "Timer" << std::right
                  << std::setw(colWidth) << "Min(s)" << std::right
                  << std::setw(colWidth) << "Max(s)" << std::right
                  << std::setw(colWidth) << "Avg(s)" << std::right
                  << std::setw(colWidth) << "% of Parent" << std::right
                  << std::setw(colWidth) << "Calls"
                  << "\n";
        std::cout << std::string(headerWidth, '-') << "\n";

        int size = 0;
        MPI_Comm_size(MPI_COMM_WORLD, &size);
        double zoneUpdates = 0.0;

        // Find root timers (those with no parents or with nesting level 0)
        std::vector<std::string> rootTimers;
        for (const auto& pair : timers_) {
            if (pair.second.nestingLevel == 0) {
                rootTimers.push_back(pair.first);
            }
        }

        // For each root timer, recursively print its hierarchy
        for (const auto& rootTimer : rootTimers) {
            printTimerAndChildren(rootTimer, timerNames, minTotals, maxTotals,
                                  sumTotals, localCalls, 0, 0.0, size);
        }

        // Find the time loop timer for zone updates calculation
        for (size_t i = 0; i < timerNames.size(); i++) {
            if (timerNames[i] == "timeLoop") {
                zoneUpdates = (maxTotals[i] > 0.0)
                                  ? (static_cast<double>(steps_) * zones_ /
                                     maxTotals[i] / 1e6)
                                  : 0.0;
                break;
            }
        }

        std::cout << std::string(headerWidth, '-') << "\n\n";
        std::cout << "Zone Updates per Second (M): " << std::fixed
                  << std::setprecision(2) << zoneUpdates << "\n\n";
    }

    // Recursive helper to print timer and its children with proper indentation
    void printTimerAndChildren(const std::string& timerName,
                               const std::vector<std::string>& timerNames,
                               const std::vector<double>& minTotals,
                               const std::vector<double>& maxTotals,
                               const std::vector<double>& sumTotals,
                               const std::vector<double>& localCalls,
                               int indent, double parentTime, int size) const {
        // Find this timer's index
        auto it = std::find(timerNames.begin(), timerNames.end(), timerName);
        if (it == timerNames.end()) return;

        size_t idx = std::distance(timerNames.begin(), it);
        double avgTotal = sumTotals[idx] / static_cast<double>(size);
        double percentOfParent =
            (parentTime > 0.0) ? (maxTotals[idx] / parentTime) * 100.0 : 0.0;

        // Indentation for hierarchical display
        std::string indentStr(indent * 2, ' ');
        std::string displayName = indentStr + timerName;

        const int nameWidth = 28;
        const int colWidth = 12;

        std::cout << std::left << std::setw(nameWidth) << displayName
                  << std::right << std::setw(colWidth) << std::fixed
                  << std::setprecision(6) << minTotals[idx] << std::right
                  << std::setw(colWidth) << maxTotals[idx] << std::right
                  << std::setw(colWidth) << avgTotal;

        if (indent > 0) {
            std::cout << std::right << std::setw(colWidth) << std::fixed
                      << std::setprecision(2) << percentOfParent;
        } else {
            std::cout << std::right << std::setw(colWidth) << "-";
        }

        std::cout << std::right << std::setw(colWidth)
                  << static_cast<long>(localCalls[idx]) << "\n";

        // Process children
        const auto& timer = timers_.find(timerName);
        if (timer != timers_.end()) {
            for (const auto& childName : timer->second.children) {
                printTimerAndChildren(childName, timerNames, minTotals,
                                      maxTotals, sumTotals, localCalls,
                                      indent + 1, maxTotals[idx], size);
            }
        }
    }

    // Map from timer name -> accumulated data
    mutable std::unordered_map<std::string, TimerData> timers_;

    // Map from timer name -> start time
    std::unordered_map<std::string, std::chrono::steady_clock::time_point>
        startTimes_;

    // Stack of currently active timers to track nesting
    std::stack<std::string> _activeTimerStack;

    // Simulation parameters for zone updates
    long steps_ = 0;
    long zones_ = 0;
    int rank_ = 0;  // Store the processor rank
    int size_ = 1;  // Store total number of ranks
};

}  // namespace agoge