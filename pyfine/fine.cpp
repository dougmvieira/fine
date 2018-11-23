#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/eigen.h>

#include "fine.hpp"

namespace py = pybind11;
using namespace fine;


PYBIND11_PLUGIN(fine) {
    py::module m("fine", "The affine stochastic process models pricing module");

    m.def("heston_formula", &heston_formula, "Vectorised Heston formula",
          py::arg("strikes"), py::arg("maturities"),
          py::arg("underlying_price"), py::arg("volatility"),
          py::arg("interest_rate"), py::arg("kappa"), py::arg("theta"),
          py::arg("nu"), py::arg("rho"));

    m.def("heston_delta", &heston_delta, "Vectorised Heston delta",
          py::arg("strikes"), py::arg("maturities"),
          py::arg("underlying_price"), py::arg("volatility"),
          py::arg("interest_rate"), py::arg("kappa"), py::arg("theta"),
          py::arg("nu"), py::arg("rho"));

    m.def("heston_vega", &heston_vega, "Vectorised Heston vega",
          py::arg("strikes"), py::arg("maturities"),
          py::arg("underlying_price"), py::arg("volatility"),
          py::arg("interest_rate"), py::arg("kappa"), py::arg("theta"),
          py::arg("nu"), py::arg("rho"));

    m.def("heston_calibration", &heston_calibration,
          "Calibrate Heston from market prices",
          py::arg("prices"), py::arg("strikes"), py::arg("maturities"),
          py::arg("underlying_price"), py::arg("interest_rate"));

    m.def("heston_calibration_ts", &heston_calibration_ts,
          "Calibrate Heston from market price time series",
          py::arg("prices"), py::arg("strikes"), py::arg("maturities"),
          py::arg("underlying_prices"), py::arg("interest_rate"));

    m.def("heston_calibrate_vol", &heston_calibrate_vol,
          "Calibrate Heston vol from market prices",
          py::arg("prices"), py::arg("strikes"), py::arg("maturities"),
          py::arg("underlying_prices"), py::arg("interest_rate"),
          py::arg("kappa"), py::arg("theta"), py::arg("sigma"), py::arg("rho"));

    return m.ptr();
}
