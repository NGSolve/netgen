#ifdef NG_PYTHON

// BEGIN EVIL HACK: Patch PyThread_get_key_value/PyThread_tss_get inside pybind11 to avoid deadlocks
// see https://github.com/pybind/pybind11/pull/1211 (please merge!)
#if defined(__GNUG__) && !defined(__clang__)
#  pragma GCC diagnostic ignored "-Wattributes"
#endif
#include <Python.h>
#include <pythread.h>
#include <pybind11/cast.h>
#undef PYBIND11_TLS_GET_VALUE
#if PY_VERSION_HEX >= 0x03070000
    inline void * PYBIND11_TLS_GET_VALUE(Py_tss_t *state) {
        PyThreadState *tstate = (PyThreadState *) PyThread_tss_get(state);
        if (!tstate) tstate = PyGILState_GetThisThreadState();
        return tstate;
    }
#else
    inline void * PYBIND11_TLS_GET_VALUE(int state) {
        PyThreadState *tstate = (PyThreadState *) PyThread_get_key_value(state);
        if (!tstate) tstate = PyGILState_GetThisThreadState();
        return tstate;
    }
#endif
// END EVIL HACK
#include <pybind11/pybind11.h>
#include <pybind11/operators.h>
#include <pybind11/numpy.h>
namespace py = pybind11;
#include <iostream>
#include <sstream>


template <typename T>
py::array MoveToNumpy(std::vector<T>& vec)
{
  auto newvec = new std::vector<T>();
  std::swap(*newvec, vec);
  auto capsule = py::capsule(newvec, [](void *v) { delete reinterpret_cast<std::vector<T>*>(v); });
  return py::array(newvec->size(), newvec->data(), capsule);
}

namespace PYBIND11_NAMESPACE {
template<typename T>
bool CheckCast( py::handle obj ) {
  try{
    obj.cast<T>();
    return true;
  }
  catch (py::cast_error &e) {
    return false;
  }
}


template <typename T>
struct extract
{
  py::handle obj;
  extract( py::handle aobj ) : obj(aobj) {}

  bool check() { return CheckCast<T>(obj); }
  T operator()() { return obj.cast<T>(); }
};
}

struct NGDummyArgument {};

inline void NOOP_Deleter(void *) { ; }

namespace netgen
{

  //////////////////////////////////////////////////////////////////////
  // Lambda to function pointer conversion
  template <typename Function>
  struct function_traits
    : public function_traits<decltype(&Function::operator())> {};

  template <typename ClassType, typename ReturnType, typename... Args>
  struct function_traits<ReturnType(ClassType::*)(Args...) const> {
    typedef ReturnType (*pointer)(Args...);
    typedef ReturnType return_type;
  };

  template <typename Function>
  typename function_traits<Function>::pointer
  FunctionPointer (const Function& lambda) {
    return static_cast<typename function_traits<Function>::pointer>(lambda);
  }


  template <class T>
  inline std::string ToString (const T& t)
  {
    std::stringstream ss;
    ss << t;
    return ss.str();
  }

}

#endif

