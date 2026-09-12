
#include "logging.hpp"
#include "python_ngcore.hpp"

namespace py = pybind11;
using std::string;

namespace ngcore
{
  bool ngcore_have_numpy = false;
  bool parallel_pickling = true;

  namespace detail
  {
    // keyed by the mangled name: type_info objects are not unique across
    // shared libraries (hidden visibility, dlls)
    static std::map<std::string, PyArchiveCasters>& PyArchiveCasterRegistry()
    {
      static std::map<std::string, PyArchiveCasters> reg;
      return reg;
    }
    void SetPyArchiveCasters(const std::type_info& ti, const PyArchiveCasters& casters)
    { PyArchiveCasterRegistry()[ti.name()] = casters; }
    const PyArchiveCasters* FindPyArchiveCasters(const std::type_info& ti)
    {
      auto& reg = PyArchiveCasterRegistry();
      auto it = reg.find(ti.name());
      return it == reg.end() ? nullptr : &it->second;
    }
    const PyArchiveCasters& GetPyArchiveCasters(const std::type_info& ti)
    {
      if (auto c = FindPyArchiveCasters(ti))
        return *c;
      throw Exception("Type " + Demangle(ti.name()) +
                      " is not registered for python archiving, call RegisterPyArchiveCaster<T>()");
    }
  } // namespace detail

  py::object CastAnyToPy(const std::any& a)
  {
    return detail::GetPyArchiveCasters(a.type()).any_to_py(a);
  }

  // registered type of the object, or of its nearest registered base class
  std::any CastPyToAny(py::object& obj)
  {
    for (auto cls : py::type::of(obj).attr("__mro__"))
      {
        auto tinfo = py::detail::get_type_info((PyTypeObject*)cls.ptr());
        if (!tinfo) continue;
        if (auto c = detail::FindPyArchiveCasters(*tinfo->cpptype))
          return c->py_to_any(obj);
      }
    throw Exception("Class " + std::string(py::str(py::type::of(obj))) +
                    " is not registered for std::any conversion, call RegisterPyArchiveCaster<T>()");
  }
  
  void SetFlag(Flags &flags, string s, py::object value) 
  {
    if (py::isinstance<py::dict>(value))
      {             
        Flags nested_flags;
        for(auto item : value.cast<py::dict>())
          SetFlag(nested_flags, item.first.cast<string>(),
                  item.second.cast<py::object>());
        flags.SetFlag(s, nested_flags);
        return;
      }
    else if (py::isinstance<py::bool_>(value))
      flags.SetFlag(s, value.cast<bool>());
    else if (py::isinstance<py::float_>(value))
      flags.SetFlag(s, value.cast<double>());
    else if (py::isinstance<py::int_>(value))
      flags.SetFlag(s, double(value.cast<int>()));
    else if (py::isinstance<py::str>(value))
      flags.SetFlag(s, value.cast<string>());
    else if (py::isinstance<py::list>(value))
      {             
        auto vdl = py::cast<py::list>(value);
        if (py::len(vdl) > 0)
          {
            if(py::isinstance<py::float_>(vdl[0]) || py::isinstance<py::int_>(vdl[0]))
              flags.SetFlag(s, makeCArray<double>(vdl));
            else if(py::isinstance<py::str>(vdl[0]))
              flags.SetFlag(s, makeCArray<string>(vdl));
            else
              {
                /*
                Array<std::any> sta;
                for (auto el : vdl)
                  // sta.Append(CastPyToAny(dynamic_cast<py::object&>(el)));
                  {
                    auto obj = py::reinterpret_borrow<py::object>(el);
                    sta.Append(CastPyToAny(obj));
                  }
                */
                std::vector<std::any> sta;
                for (auto el : vdl)
                  // sta.Append(CastPyToAny(dynamic_cast<py::object&>(el)));
                  {
                    auto obj = py::reinterpret_borrow<py::object>(el);
                    sta.push_back(CastPyToAny(obj));
                  }
                
                flags.SetFlag(s, sta);
              }
          }
        else
          {
            Array<string> dummystr;
            Array<double> dummydbl;
            Array<std::any> dummyany;
            flags.SetFlag(s,dummystr);
            flags.SetFlag(s,dummydbl);
            flags.SetFlag(s,dummyany);            
          }
      }
    else if (py::isinstance<py::tuple>(value))
      {
        auto vdt = py::cast<py::tuple>(value);
        if (py::isinstance<py::float_>(value))
          flags.SetFlag(s, makeCArray<double>(vdt));
        if (py::isinstance<py::int_>(value))
          flags.SetFlag(s, makeCArray<double>(vdt));
        if (py::isinstance<py::str>(value))
          flags.SetFlag(s, makeCArray<string>(vdt));
      }
    else
      {
        flags.SetFlag(s, CastPyToAny(value));
      }
  }

  Flags CreateFlagsFromKwArgs(const py::kwargs& kwargs, py::object pyclass, py::list info)
  {
    static std::shared_ptr<Logger> logger = GetLogger("Flags");
    py::dict flags_dict;

    if (kwargs.contains("flags"))
      {
        logger->warn("WARNING: using flags as kwarg is deprecated{}, use the flag arguments as kwargs instead!",
                     pyclass.is_none() ? "" : std::string(" in ") + std::string(py::str(pyclass)));
        auto addflags = py::cast<py::dict>(kwargs["flags"]);
        for (auto item : addflags)
          flags_dict[item.first.cast<string>().c_str()] = item.second;
      }
    py::dict special;
    if(!pyclass.is_none())
      {
        auto flags_doc = pyclass.attr("__flags_doc__")();
        for (auto item : kwargs)
          if (!flags_doc.contains(item.first.cast<string>().c_str()) &&
              !(item.first.cast<string>() == "flags"))
            logger->warn("WARNING: kwarg '{}' is an undocumented flags option for class {}, maybe there is a typo?",
                         item.first.cast<string>(), std::string(py::str(pyclass)));
      
        if(py::hasattr(pyclass,"__special_treated_flags__"))
          special = pyclass.attr("__special_treated_flags__")();
      }
    for (auto item : kwargs)
      {
        auto name = item.first.cast<string>();
        if (name != "flags")
          {
            if(!special.contains(name.c_str()))
              flags_dict[name.c_str()] = item.second;
          }
      }

    Flags flags;
    for(auto item : flags_dict)
      SetFlag(flags, item.first.cast<string>(), item.second.cast<py::object>());

    for (auto item : kwargs)
      {
        auto name = item.first.cast<string>();
        if (name != "flags")
          {
            if(special.contains(name.c_str()))
              special[name.c_str()](item.second, &flags, info);
          }
      }
    return flags;
  }

  py::dict CreateDictFromFlags(const Flags& flags)
  {
    py::dict d;
    std::string key;
    for(auto i : Range(flags.GetNFlagsFlags()))
      {
        auto& f = flags.GetFlagsFlag(i, key);
        d[key.c_str()] = CreateDictFromFlags(f);
      }
    for(auto i : Range(flags.GetNStringListFlags()))
      {
        auto strlistflag = flags.GetStringListFlag(i, key);
        py::list lst;
        for(auto& val : *strlistflag)
          lst.append(val);
        d[key.c_str()] = lst;
      }
    for(auto i : Range(flags.GetNNumListFlags()))
      {
        auto numlistflag = flags.GetNumListFlag(i, key);
        py::list lst;
        for(auto& val : *numlistflag)
          lst.append(val);
        d[key.c_str()] = lst;
      }
    for(auto i : Range(flags.GetNStringFlags()))
      {
        auto val = flags.GetStringFlag(i, key);
        d[key.c_str()] = val;
      }
    for(auto i : Range(flags.GetNNumFlags()))
      {
        auto val = flags.GetNumFlag(i, key);
        d[key.c_str()] = val;
      }
    for(auto i : Range(flags.GetNDefineFlags()))
      {
        auto val = flags.GetDefineFlag(i, key);
        d[key.c_str()] = val;
      }
    for(auto i : Range(flags.GetNAnyFlags()))
      {
        auto& a = flags.GetAnyFlag(i, key);
        d[key.c_str()] = CastAnyToPy(a);
      }
    return d;
  }

} // namespace ngcore
