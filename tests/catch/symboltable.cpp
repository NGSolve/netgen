
#include "catch.hpp"
#include <../core/ngcore.hpp>
using namespace ngcore;
using namespace std;

TEST_CASE("Symboltable")
{
  SymbolTable<int> table;
  CHECK_THROWS_AS(table["test"], RangeException);
  table.Set("first", 1);
  CHECK(table["first"] == 1);
  table["first"] = 2;
  CHECK(table["first"] == 2);
  auto index = table.Index("first");
  CHECK(index == 0);
  CHECK(table[index] == 2);
  table[index] = 3;
  CHECK(table["first"] == 3);
#ifndef NDEBUG
  int a;
  CHECK_THROWS_AS(a = table[5], RangeException);
  CHECK_THROWS_AS(table.GetName(5), RangeException);
#endif
  CHECK(table.Used("first"));
  CHECK(!table.Used("second"));
  SymbolTable<int> another;
  another.Set("first", 1);
  another.Set("second", 2);
  table.Update(another);
  CHECK(table["first"] == 1);
  CHECK(table["second"] == 2);
  std::stringstream s;
  s << table;
  CHECK(s.str() == "first : 1\nsecond : 2\n");
  auto ss = std::make_shared<std::stringstream>();
  BinaryOutArchive ao(ss);
  ao & table;
  ao.FlushBuffer();
  BinaryInArchive ai(ss);
  SymbolTable<int> read;
  ai & read;
  for(size_t i = 0; i<table.Size(); i++)
    {
      CHECK(read[i] == table[i]);
      CHECK(read.GetName(i) == table.GetName(i));
    }
  table.DeleteAll();
  CHECK(table.Size() == 0);
  // SymbolTable<bool> is special because of vector<bool> is special...
  SymbolTable<bool> btable;
  btable.Set("true", true);
  btable.Set("false", false);
  CHECK(btable[0]);
  CHECK(!btable[1]);
  CHECK(btable["true"]);
  CHECK(!btable["false"]);
  ao & btable;
  ao.FlushBuffer();
  SymbolTable<bool> bread;
  ai & bread;
  CHECK(bread["true"]);
  CHECK(!bread["false"]);
}
