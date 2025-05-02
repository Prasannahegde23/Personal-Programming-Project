#pragma once
#include <chrono>
#include <iostream>
#include <mutex>
#include <string>
#include <unordered_map>

struct ScopedTimer 
{
  using clock = std::chrono::high_resolution_clock;
  const char* name;
  clock::time_point start;

  // shared accumulator (time in seconds)
  static std::mutex                                      mtx;
  static std::unordered_map<std::string,double> acc;

  ScopedTimer(const char* n)
    : name(n), start(clock::now()) {}

  ~ScopedTimer()
  {
    auto dur = std::chrono::duration<double>(clock::now() - start).count();
    std::lock_guard<std::mutex> l(mtx);
    acc[name] += dur;
  }

  /// Call once at very end to dump all timers:
  static void report()
  {
    std::lock_guard<std::mutex> l(mtx);
    std::cout << "\n=== Timing Summary ===\n";
    for(auto &p : acc)
      std::cout << std::setw(20) << std::left << p.first
                << ": " << p.second << " s\n";

    double h = acc["Hydrogen_subcycle"];
    double m = acc["Mechanical_PD"];
    double d = acc["Displacement"];
    std::cout << std::setw(20) << std::left << "Total_time"
              << ": " << (h + m + d) << " s\n";
  }
};

// static member definitions
std::mutex ScopedTimer::mtx;
std::unordered_map<std::string,double> ScopedTimer::acc;
