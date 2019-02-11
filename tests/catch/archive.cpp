
#include "catch.hpp"
#include <../core/ngcore.hpp>
using namespace ngcore;
using namespace std;

class CommonBase
{
public:
  int a;
  virtual ~CommonBase() {}

  virtual void DoArchive(Archive& archive) { archive & a; }
};

// pure abstract base class
class SharedPtrHolder : virtual public CommonBase
{
public:
  vector<shared_ptr<string>> names;
  virtual ~SharedPtrHolder()
  { }

  virtual void abstract() = 0;
  virtual void DoArchive(Archive& archive)
  {
    CommonBase::DoArchive(archive);
    archive & names;
  }
};

class PtrHolder : virtual public CommonBase
{
public:
  vector<int*> numbers;
  virtual ~PtrHolder() {}

  virtual void DoArchive(Archive& archive)
  {
    CommonBase::DoArchive(archive);
    archive & numbers;
  }
};

class SharedPtrAndPtrHolder : public SharedPtrHolder, public PtrHolder
{
public:
  virtual ~SharedPtrAndPtrHolder() {}
  virtual void DoArchive(Archive& archive)
  {
    SharedPtrHolder::DoArchive(archive);
    PtrHolder::DoArchive(archive);
  }
  virtual void abstract() {}
};

// Classes without virt. or multiple inheritance do not need to be registered
class SimpleClass : public CommonBase
{
public:
  double d;
  virtual void DoArchive(Archive& ar)
  {
    CommonBase::DoArchive(ar);
    ar & d;
  }
};

class NotRegisteredForArchive : public SharedPtrAndPtrHolder {};

class ClassWithConstPtr
{
private:
  const int* ptr;
public:
  ClassWithConstPtr(const int* aptr) : ptr(aptr) { }
  // constructor only for archive
  ClassWithConstPtr() {}
  void DoArchive(Archive& ar)
  {
    ar & ptr;
  }
  const int* getPtr() { return ptr; }
};

class OneMoreDerivedClass : public SharedPtrAndPtrHolder {};

static RegisterClassForArchive<CommonBase> regb;
static RegisterClassForArchive<SharedPtrHolder, CommonBase> regsp;
static RegisterClassForArchive<PtrHolder, CommonBase> regp;
static RegisterClassForArchive<SharedPtrAndPtrHolder, SharedPtrHolder, PtrHolder> regspp;
static RegisterClassForArchive<OneMoreDerivedClass, SharedPtrAndPtrHolder> regom;

void testNullPtr(Archive& in, Archive& out)
{
  SharedPtrHolder* p = nullptr;
  shared_ptr<string> sp = nullptr;
  out & p & sp;
  out.FlushBuffer();
  SharedPtrHolder* pin = nullptr;
  shared_ptr<string> spin = nullptr;
  in & pin & spin;
  CHECK(pin == nullptr);
  CHECK(spin == nullptr);
}

void testSharedPointer(Archive& in, Archive& out)
{
  SECTION("Same shared ptr")
    {
      static_assert(detail::has_DoArchive<SharedPtrHolder>::value, "");
      SharedPtrAndPtrHolder holder, holder2;
      holder.names.push_back(make_shared<string>("name"));
      holder2.names = holder.names; // same shared ptr
      out & holder & holder2;
      out.FlushBuffer();
      SharedPtrAndPtrHolder inholder, inholder2;
      in & inholder & inholder2;
      CHECK(inholder.names.size() == 1);
      CHECK(inholder.names[0] == inholder2.names[0]);
      CHECK(inholder.names[0].use_count() == 3); // one shared ptr is still kept in the archive
      CHECK(*inholder.names[0] == "name");
    }
}

void testPointer(Archive& in, Archive& out)
{
  PtrHolder holder, holder2;
  holder.numbers.push_back(new int(3));
  holder2.numbers = holder.numbers; // same shared ptr
  out & holder & holder2;
  out.FlushBuffer();
  PtrHolder inholder, inholder2;
  in & inholder & inholder2;
  CHECK(inholder.numbers.size() == 1);
  CHECK(inholder.numbers[0] == inholder2.numbers[0]);
  CHECK(*inholder.numbers[0] == 3);
}

void testConstPointer(Archive& in, Archive& out)
{
  SECTION("Const pointer")
    {
      int* iptr = new int(4);
      double d = 0.1;
      ClassWithConstPtr cls(iptr);
      out & cls & iptr & d;
      out.FlushBuffer();
      ClassWithConstPtr incls;
      int* iniptr;
      double ind;
      in & incls & iniptr & ind;
      CHECK(*incls.getPtr() == 4);
      CHECK(incls.getPtr() == iniptr);
      CHECK(ind == 0.1);
      delete iptr;
      delete iniptr;
    }
}

void testMultipleInheritance(Archive& in, Archive& out)
{
  PtrHolder* p = new OneMoreDerivedClass;
  p->numbers.push_back(new int(2));
  p->a = 5;
  auto p2 = dynamic_cast<SharedPtrHolder*>(p);
  p2->names.push_back(make_shared<string>("test"));
  auto sp1 = shared_ptr<PtrHolder>(p);
  auto sp2 = dynamic_pointer_cast<SharedPtrHolder>(sp1);
  auto checkPtr = [] (auto pin, auto pin2)
  {
      CHECK(typeid(*pin) == typeid(*pin2));
      CHECK(typeid(*pin) == typeid(OneMoreDerivedClass));
      CHECK(*pin2->names[0] == "test");
      CHECK(*pin->numbers[0] == 2);
      CHECK(dynamic_cast<SharedPtrAndPtrHolder*>(pin) == dynamic_cast<SharedPtrAndPtrHolder*>(pin2));
      CHECK(pin->a == pin2->a);
      CHECK(pin->a == 5);
      REQUIRE(dynamic_cast<SharedPtrAndPtrHolder*>(pin2) != nullptr);
      CHECK(*dynamic_cast<SharedPtrAndPtrHolder*>(pin2)->numbers[0] == 2);
      CHECK(*pin->numbers[0] == *dynamic_cast<SharedPtrAndPtrHolder*>(pin2)->numbers[0]);
      REQUIRE(dynamic_cast<SharedPtrAndPtrHolder*>(pin) != nullptr);
      CHECK(dynamic_cast<SharedPtrAndPtrHolder*>(pin)->names[0] == pin2->names[0]);
  };
  SECTION("Archive ptrs to leaves of mult. inh.")
    {
      out & p & p2;
      out.FlushBuffer();
      PtrHolder* pin = nullptr;
      SharedPtrHolder* pin2 = nullptr;
      in & pin & pin2;
      checkPtr(pin, pin2);
    }
  SECTION("Archive shared ptrs to leaves of mult. inh.")
    {
      out & sp1 & sp2;
      out.FlushBuffer();
      shared_ptr<PtrHolder> pin;
      shared_ptr<SharedPtrHolder> pin2;
      in & pin & pin2;
      checkPtr(pin.get(), pin2.get());
    }
  SECTION("Virtual base class")
    {
      CommonBase* b = dynamic_cast<CommonBase*>(p);
      out & b & p;
      PtrHolder* pin;
      CommonBase* bin;
      in & bin & pin;
      checkPtr(pin, dynamic_cast<SharedPtrHolder*>(bin));
    }
  SECTION("Simple class without register")
    {
      auto a = new SimpleClass;
      a->a = 5;
      a->d = 2.3;
      SECTION("check pointer")
        {
          out & a;
          out.FlushBuffer();
          SimpleClass* ain;
          in & ain;
          CHECK(ain->a == 5);
          CHECK(ain->d == 2.3);
        }
      SECTION("check shared pointer")
        {
          auto spa = shared_ptr<SimpleClass>(a);
          out & spa;
          out.FlushBuffer();
          shared_ptr<SimpleClass> spain;
          in & spain;
          CHECK(spain->a == 5);
          CHECK(spain->d == 2.3);
        }
    }
}

void testMap(Archive& in, Archive& out)
{
  map<string, VersionInfo> map1;
  map1["netgen"] = "v6.2.1901";
  out & map1;
  out.FlushBuffer();
  map<string, VersionInfo> map2;
  in & map2;
  CHECK(map2.size() == 1);
  CHECK(map2["netgen"] == "v6.2.1901");
}

enum MyEnum
  {
   CASE1,
   CASE2
  };

void testEnum(Archive& in, Archive& out)
  {
   MyEnum en = CASE2;
   out & en;
   out.FlushBuffer();
   MyEnum enin;
   in & enin;
   CHECK(enin == CASE2);
  }

void testArchive(Archive& in, Archive& out)
{
  SECTION("Empty String")
    {
      char* cstr = nullptr;
      char* empty = new char[1];
      char* simple = new char[7] {'s','i','m','p','l','e','\0'};
      empty[0] = '\0';
      out << string("") << cstr << empty << simple << string("simple") << long(1);
      out.FlushBuffer();
      string str; long i; char* readempty; char* readsimple;
      string simplestr;
      in & str & cstr & readempty & readsimple & simplestr & i;
      CHECK(str == "");
      CHECK(cstr == nullptr);
      CHECK(strcmp(readempty,"") == 0);
      CHECK(strcmp(readsimple,"simple") == 0);
      CHECK(i == 1);
      CHECK(simplestr == "simple");
      delete[] readempty;
      delete[] empty;
      delete[] simple;
      delete[] readsimple;
    }
  SECTION("SharedPtr")
    {
      testSharedPointer(in, out);
    }
  SECTION("Pointer")
    {
      testPointer(in, out);
    }
  SECTION("Const Pointer")
    {
      testConstPointer(in, out);
    }
  SECTION("Multiple inheritance")
    {
      testMultipleInheritance(in, out);
    }
  SECTION("Not registered")
    {
      SharedPtrAndPtrHolder* p = new NotRegisteredForArchive;
      REQUIRE_THROWS(out & p, Catch::Contains("not registered for archive"));
    }
  SECTION("nullptr")
    {
      testNullPtr(in, out);
    }
  SECTION("map")
    {
      testMap(in, out);
    }
  SECTION("enum")
    {
      testEnum(in, out);
    }
}

TEST_CASE("BinaryArchive")
{
  auto stream = make_shared<stringstream>();
  BinaryOutArchive out(stream);
  BinaryInArchive in(stream);
  testArchive(in, out);
}

TEST_CASE("TextArchive")
{
  auto stream = make_shared<stringstream>();
  TextOutArchive out(stream);
  TextInArchive in(stream);
  testArchive(in, out);
}
