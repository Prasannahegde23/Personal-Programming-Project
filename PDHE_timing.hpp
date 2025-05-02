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
      std::cout << p.first
                << std::string(20 - p.first.length(), ' ')
                << ": " << p.second << " s\n";
  }
};

// static member definitions
std::mutex ScopedTimer::mtx;
std::unordered_map<std::string,double> ScopedTimer::acc;
