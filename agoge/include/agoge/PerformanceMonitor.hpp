#pragma once

#include <string>
#include <unordered_map>
#include <chrono>
#include <iostream>

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
 *   PerformanceMonitor::instance().printReport();
 */
class PerformanceMonitor
{
public:
    /**
     * @brief Access the global singleton instance of PerformanceMonitor.
     */
    static PerformanceMonitor& instance()
    {
        static PerformanceMonitor s_instance;
        return s_instance;
    }

    /**
     * @brief Start a timer for the given name. If already running, it restarts.
     *
     * @param name Unique string identifying this timer (e.g., "computeL").
     */
    void startTimer(const std::string &name)
    {
        auto now = std::chrono::steady_clock::now();
        startTimes_[name] = now;
    }

    /**
     * @brief Stop a timer for the given name, accumulate the elapsed time.
     *
     * @param name Must match a previously started timer.
     */
    void stopTimer(const std::string &name)
    {
        auto now = std::chrono::steady_clock::now();
        auto it = startTimes_.find(name);
        if(it != startTimes_.end()) {
            auto elapsed = now - it->second;
            timers_[name].total += elapsed;
            timers_[name].callCount += 1;
            // Erase the start time or leave it for next iteration
            // startTimes_.erase(it);
        }
    }

    /**
     * @brief Print a summary of all timers to stdout.
     *
     * Shows total time in seconds, number of calls, and average time in ms.
     */
    void printReport() const
    {
        std::cout << "\n================= PERFORMANCE REPORT =================\n";
        std::cout << "Timer Name             Total Time (s)   Calls   Avg (ms)\n";
        std::cout << "-------------------------------------------------------\n";

        for(const auto &pair : timers_) {
            const auto &name = pair.first;
            const auto &data = pair.second;

            double totalSec = std::chrono::duration<double>(data.total).count();
            double avgMs = 0.0;
            if(data.callCount > 0) {
                avgMs = (totalSec * 1000.0) / data.callCount;
            }
            std::cout << name << "  "
                      << totalSec << "          "
                      << data.callCount << "       "
                      << avgMs << "\n";
        }
        std::cout << "=======================================================\n\n";
    }

private:
    PerformanceMonitor() = default; // Private constructor for singleton
    PerformanceMonitor(const PerformanceMonitor&) = delete;
    PerformanceMonitor& operator=(const PerformanceMonitor&) = delete;

    // Map from timer name -> accumulated data
    std::unordered_map<std::string, TimerData> timers_;

    // Map from timer name -> start time
    std::unordered_map<std::string, std::chrono::steady_clock::time_point> startTimes_;
};

} // namespace agoge
