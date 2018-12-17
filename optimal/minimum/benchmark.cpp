#include "benchmark.h"

Benchmark::Benchmark()
{
    m_tick = m_tock = 0;
}

auto Benchmark::tick() -> std::clock_t
{
    m_chrone_tick = std::chrono::high_resolution_clock::now();
    m_tick = clock();
    return m_tick;
}

auto Benchmark::tock() -> std::clock_t
{
    m_chrone_tock = std::chrono::high_resolution_clock::now();
    m_tock = clock();
    return m_tock;
}

auto Benchmark::CpuDurationClock() const -> std::clock_t
{
    return m_tock - m_tick;
}

auto Benchmark::CpuDurationSecond() const -> double
{
    return static_cast<double>(m_tock - m_tick)/CLOCKS_PER_SEC;
}

auto Benchmark::printDuration() const -> void
{
    std::cout << std::fixed << std::setprecision(2)
              << "CPU time used: "
              << CpuDurationSecond() << " ms" << " "
              << "Wall clock time passed: "
              << std::chrono::duration<double, std::milli>(m_chrone_tock-m_chrone_tick).count() << " ms"
              << std::endl;
}

auto Benchmark::printCpuClockDuration() const -> void
{
    std::cout << std::fixed << std::setprecision(2)
              << "CPU time used: "
              << CpuDurationSecond() << " ms" << " "
              << std::endl;
}

auto Benchmark::printWallClocDuration() const -> void
{
    std::cout << std::fixed << std::setprecision(2)
              << "Wall clock time passed: "
              << std::chrono::duration<double, std::milli>(m_chrone_tock-m_chrone_tick).count() << " ms"
              << std::endl;
}

