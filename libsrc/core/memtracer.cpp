#include <filesystem>
#include <vector>
#include <string>
#include <iostream>
#include <fstream>
#include <tuple>

#include <iomanip>
#include <sstream>

#include "memtracer.hpp"


namespace ngcore
{

size_t GetPageSize()
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

size_t GetRSSMemory()
{
    static size_t page_size = GetPageSize();
    if(page_size == 0)
        return 0;
    std::ifstream statm("/proc/self/statm");
    long size = 0, resident = 0;
    statm >> size >> resident;
    return resident * page_size;
}

void PrintMemoryUsage(const char  * file, int line, std::string msg, std::ostream & out)
{
    using std::setw;
    using std::fixed;
    using std::showpos;
    using std::noshowpos;
    using std::setprecision;

    static double last_rss_memory = 0;

    double to_mb = 1.0/(1024.0 * 1024.0);
    double mem_rss = GetRSSMemory() * to_mb;
    double diff_rss = mem_rss-last_rss_memory;
    last_rss_memory = mem_rss;

    if(mem_rss > 0)
    {
        out  << "rss: " << setw(9) << fixed << setprecision(1) << mem_rss << " MB"
             << "  diff_rss: " << showpos << setw(9) << fixed << setprecision(1) << diff_rss << " MB";
    }
    out << noshowpos << "  " << file << ":" << line << " " << msg << std::endl;
}

} // namespace ngcore

