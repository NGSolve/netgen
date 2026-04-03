#include <meshing.hpp>
#include "writeuser.hpp"

namespace netgen
{
DLL_HEADER void ReadMeditFormat (Mesh & mesh, const filesystem::path & filename, map<tuple<int, int>, int> & index_map);
DLL_HEADER void ReadMeditFormat (Mesh & mesh, const filesystem::path & filename);

DLL_HEADER void WriteMeditFormat (const Mesh & mesh, const filesystem::path & filename, map<tuple<int,int>, int> & index_map);
DLL_HEADER void WriteMeditFormat (const Mesh & mesh, const filesystem::path & filename);
} // namespace netgen
