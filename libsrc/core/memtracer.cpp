#include <filesystem>
#include <vector>
#include <string>
#include <iostream>
#include <fstream>
#include <tuple>

#include "array.hpp"
#include "memtracer.hpp"


namespace ngcore
{

#if defined(NETGEN_TRACE_MEMORY) && !defined(__CUDA_ARCH__)
std::vector<std::string> MemoryTracer::names{"all"};
std::vector<int> MemoryTracer::parents{-1};
std::atomic<size_t> MemoryTracer::total_memory{0};
std::mutex MemoryTracer::create_id_mutex;
#endif

size_t MemoryTracer :: GetPageSize()
{
    if(!std::filesystem::exists("/proc/self/smaps") || !std::filesystem::exists("/proc/self/statm"))
        return 0;

    std::ifstream smaps("/proc/self/smaps");
    std::string line;
    while (std::getline(smaps, line))
    {
        if (line.rfind("KernelPageSize:", 0) == 0)
        {
            std::size_t value = 0;
            std::string unit;
            std::istringstream iss(line);
            iss >> unit >> value >> unit;
            if(unit == "kB" || unit == "KB")
                value = value * 1024;
            return value;
        }
    }
    return 0;
}

size_t MemoryTracer :: GetRSSMemory()
{
    static size_t page_size = GetPageSize();
    if(page_size == 0)
        return 0;
    std::ifstream statm("/proc/self/statm");
    long size = 0, resident = 0;
    statm >> size >> resident;
    return resident * page_size;
}

void MemoryTracer :: PrintMemoryUsage(const char  * file, int line, std::string msg, std::ostream & out)
{
    using std::setw;
    using std::fixed;
    using std::showpos;
    using std::noshowpos;
    using std::setprecision;

    static double last_memory = 0;
    static double last_rss_memory = 0;

    double to_mb = 1.0/(1024.0 * 1024.0);

    double memory = GetTotalMemory() * to_mb;
    double diff = memory-last_memory;

    double mem_rss = GetRSSMemory() * to_mb;
    double diff_rss = mem_rss-last_rss_memory;
    last_rss_memory = mem_rss;

    out << "mem: " << setw(9) << fixed << setprecision(1) << memory << " MB"
        << "  diff: " << showpos << setw(9) << fixed << setprecision(1) << diff << " MB  ";
    if(mem_rss > 0)
    {
        out  << "rss: " << noshowpos << setw(9) << fixed << setprecision(1) << mem_rss << " MB"
             << "  diff_rss: " << showpos << setw(9) << fixed << setprecision(1) << diff_rss << " MB";
    }
    out << noshowpos << "  " << file << ":" << line << " " << msg << std::endl;

    last_memory = memory;
}

} // namespace ngcore

