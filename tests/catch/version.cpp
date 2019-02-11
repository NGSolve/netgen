
#include "catch.hpp"
#include <../core/ngcore.hpp>
using namespace ngcore;
using namespace std;


TEST_CASE("Version")
{
  VersionInfo v("v6.2.1811-3-asdf");
  CHECK(v.to_string() == "v6.2.1811-3-asdf");
  VersionInfo v2("6.2");
  CHECK(v2.to_string() == "v6.2");
  CHECK(v < "v7");
  CHECK(v >= "6.2");
  CHECK(v > "6.2.1811");
  CHECK(v < "6.2.1811-5");
  CHECK(v == "v6.2.1811-3");
}
