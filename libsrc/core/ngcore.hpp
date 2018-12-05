#ifndef NG_CORE_HPP
#define NG_CORE_HPP

// std includes
#include <memory>
#include <string>
#include <vector>
#include <map>
#include <iterator>
#include <type_traits>
#include <functional>
#include <stdexcept>
#include <iostream>
#include <fstream>
#include <cstring>
#include <complex>


namespace ngcore
{
#if defined(__GNUC__)
  inline bool likely (bool x) { return __builtin_expect((x), true); }
  inline bool unlikely (bool x) { return __builtin_expect((x), false); }
#else
  inline bool likely (bool x) { return x; }
  inline bool unlikely (bool x) { return x; }
#endif
}

// own includes
#include "basearchive.hpp"
#include "version.hpp"
#include "archive.hpp"

#endif // NG_CORE_HPP
