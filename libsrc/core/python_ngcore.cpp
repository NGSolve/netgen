
#include "logging.hpp"
#include "python_ngcore.hpp"

namespace py = pybind11;
using std::string;

namespace ngcore
{
  bool ngcore_have_numpy = false;
  bool parallel_pickling = true;
  
  void SetFlag(Flags &flags, string s, py::object value) 
  {
    if (py::isinstance<py::dict>(value))
      {             
        py::dict vdd(value);
        // call recursively to set dictionary
        for (auto item : vdd) {
          string name = item.first.cast<string>();
          py::object val = py::reinterpret_borrow<py::object>(item.second);
          SetFlag(flags, name, val);
        }
        return;
      }

    if (py::isinstance<py::bool_>(value))
      flags.SetFlag(s, value.cast<bool>());

    if (py::isinstance<py::float_>(value))
      flags.SetFlag(s, value.cast<double>());

    if (py::isinstance<py::int_>(value))
      flags.SetFlag(s, double(value.cast<int>()));

    if (py::isinstance<py::str>(value))
      flags.SetFlag(s, value.cast<string>());

    if (py::isinstance<py::list>(value))
      {             
        auto vdl = py::cast<py::list>(value);
        if (py::len(vdl) > 0)
          {
            if(py::isinstance<py::float_>(vdl[0]) || py::isinstance<py::int_>(vdl[0]))
              flags.SetFlag(s, makeCArray<double>(vdl));
            if(py::isinstance<py::str>(vdl[0]))
              flags.SetFlag(s, makeCArray<string>(vdl));
          }
        else
          {
            Array<string> dummystr;
            Array<double> dummydbl;
            flags.SetFlag(s,dummystr);
            flags.SetFlag(s,dummydbl);
          }
      }

    if (py::isinstance<py::tuple>(value))
      {
        auto vdt = py::cast<py::tuple>(value);
        if (py::isinstance<py::float_>(value))
          flags.SetFlag(s, makeCArray<double>(vdt));
        if (py::isinstance<py::int_>(value))
          flags.SetFlag(s, makeCArray<double>(vdt));
        if (py::isinstance<py::str>(value))
          flags.SetFlag(s, makeCArray<string>(vdt));
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

    auto flags = py::cast<Flags>(flags_dict);

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
    return d;
  }

} // namespace ngcore
