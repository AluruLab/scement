#ifndef PY_TYPES_HPP
#define PY_TYPES_HPP 

#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>

typedef pybind11::array_t<float,
                          pybind11::array::f_style | pybind11::array::forcecast>
    pyfarray_t;
typedef pybind11::array_t<int32_t,
                          pybind11::array::f_style | pybind11::array::forcecast>
    pyiarray_t;
typedef pybind11::dict pydict;
typedef pybind11::options pyoptions;
typedef pybind11::buffer_info pybuffer_info;
typedef pybind11::arg pyarg;

#endif // !PY_TYPES_HPP
