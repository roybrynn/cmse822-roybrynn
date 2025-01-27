#pragma once

#include <chrono>
#include <iomanip>
#include <iostream>
#include <string>
#include <unordered_map>

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
};

/**
 * @class PerformanceMonitor
 * @brief A singleton class that manages named timers for performance tracking.
 *
 * Usage:
 *   PerformanceMonitor::instance().startTimer("computeL");
 *   // ... code ...
 *   PerformanceMonitor::instance().stopTimer("computeL");
 *
 *   At program end:
 *   PerformanceMonitor::instance().setSteps(steps);
 *   PerformanceMonitor::instance().setZones(zones);
 *   PerformanceMonitor::instance().printReport();
 */
class PerformanceMonitor {
   public:
    /**
     * @brief Access the global singleton instance of PerformanceMonitor.
     */
    static PerformanceMonitor &instance() {
        static PerformanceMonitor s_instance;
        return s_instance;
    }

    /**
     * @brief Start a timer for the given name. If already running, it restarts.
     *
     * @param name Unique string identifying this timer (e.g., "computeL").
     */
    void startTimer(const std::string &name) {
        auto now = std::chrono::steady_clock::now();
        startTimes_[name] = now;
    }

    /**
     * @brief Stop a timer for the given name, accumulate the elapsed time.
     *
     * @param name Must match a previously started timer.
     */
    void stopTimer(const std::string &name) {
        auto now = std::chrono::steady_clock::now();
        auto it = startTimes_.find(name);
        if (it != startTimes_.end()) {
            auto elapsed = now - it->second;
            timers_[name].total += elapsed;
            timers_[name].callCount += 1;
            // Optionally erase the start time
            // startTimes_.erase(it);
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
     * @brief Print a summary of all timers to stdout.
     *
     * Shows total time in seconds, number of calls, average time in ms,
     * and zone updates per second.
     */
    void printReport() const {
        // Define column widths
        const int nameWidth = 20;  // Adjust as needed
        const int timeWidth = 18;
        const int callsWidth = 10;
        const int avgWidth = 13;
        const int zupsWidth = 16;  // Zone Updates per Second

        std::cout
            << "\n================= PERFORMANCE REPORT =================\n";

        // Set formatting for headers
        std::cout << std::left << std::setw(nameWidth) << "Timer Name"
                  << std::right << std::setw(timeWidth) << "Total Time (s)"
                  << std::right << std::setw(callsWidth) << "Calls"
                  << std::right << std::setw(avgWidth) << "Avg (ms)"
                  << std::right << std::setw(zupsWidth) << "Mega Zone Updates/s"
                  << "\n";

        std::cout << "---------------------------------------------------------"
                     "-----------------\n";

        // Set formatting for data rows
        for (const auto &pair : timers_) {
            const auto &name = pair.first;
            const auto &data = pair.second;

            double totalSec = std::chrono::duration<double>(data.total).count();
            double avgMs = 0.0;
            double MzoneUpdatesPerSec = 0.0;

            if (data.callCount > 0) {
                avgMs = (totalSec * 1000.0) / data.callCount;
            }

            if (totalSec > 0.0) {
                MzoneUpdatesPerSec = (static_cast<double>(steps_) *
                                     static_cast<double>(zones_)) /
                                    totalSec / 1e6;
            } else {
                MzoneUpdatesPerSec = 0.0;
            }

            std::cout << std::left << std::setw(nameWidth) << name << std::right
                      << std::setw(timeWidth) << std::fixed
                      << std::setprecision(6) << totalSec << std::right
                      << std::setw(callsWidth) << data.callCount << std::right
                      << std::setw(avgWidth) << std::fixed
                      << std::setprecision(6) << avgMs << std::right
                      << std::setw(zupsWidth) << std::fixed
                      << std::setprecision(2) << MzoneUpdatesPerSec << "\n";
        }

        std::cout << "---------------------------------------------------------"
                     "-----------------\n\n";
    }

   private:
    PerformanceMonitor() = default;  // Private constructor for singleton
    PerformanceMonitor(const PerformanceMonitor &) = delete;
    PerformanceMonitor &operator=(const PerformanceMonitor &) = delete;

    // Map from timer name -> accumulated data
    std::unordered_map<std::string, TimerData> timers_;

    // Map from timer name -> start time
    std::unordered_map<std::string, std::chrono::steady_clock::time_point>
        startTimes_;

    // Simulation parameters for zone updates
    long steps_ = 0;
    long zones_ = 0;
};

}  // namespace agoge
