#ifndef BENCHMARK_H
#define BENCHMARK_H

#include <iostream>
#include <iomanip>
#include <ctime>
#include <chrono>
#include "global.h"

using namespace std;

class MINIMUMSHARED_EXPORT Benchmark
{
public:
    Benchmark();

    auto tick() -> std::clock_t;
    auto tock() -> std::clock_t;

    auto CpuDurationClock() const -> std::clock_t;
    auto CpuDurationSecond() const -> double;
    auto printDuration() const -> void;
    auto printCpuClockDuration() const -> void;
    auto printWallClockDuration() const -> void;

private:
    std::clock_t m_tick;
    std::clock_t m_tock;
    //std::chrono::time_point< std::chrono::high_resolution_clock > m_chrone_tick;
    //std::chrono::time_point< std::chrono::high_resolution_clock > m_chrone_tock;
};

#endif // BENCHMARK_H
