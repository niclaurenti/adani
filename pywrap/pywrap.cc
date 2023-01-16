#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/operators.h>
#include <pybind11/functional.h>

#include <adani/adani.h>

namespace py = pybind11;

PYBIND11_MODULE(adanipy, m) {
    // Documentation
    m.doc() = "Python wrapper of libadani";

    m.def("C2m_g3_approximation", &C2m_g3_approximation, "", py::arg("x"), py::arg("mQ"), py::arg("mMu"), py::arg("nf"), py::arg("method_flag")=1, py::arg("calls")=500000);
    m.def("CLm_g3_approximation", &CLm_g3_approximation, "", py::arg("x"), py::arg("mQ"), py::arg("mMu"), py::arg("nf"), py::arg("method_flag")=1, py::arg("calls")=500000);
    m.def("C2m_ps3_approximation", &C2m_ps3_approximation, "", py::arg("x"), py::arg("mQ"), py::arg("mMu"), py::arg("nf"));
    m.def("CLm_ps3_approximation", &CLm_ps3_approximation, "", py::arg("x"), py::arg("mQ"), py::arg("mMu"), py::arg("nf"));

    m.def("C2m_g3_approximationA_klmv", &C2m_g3_approximationA_klmv, "", py::arg("x"), py::arg("mQ"), py::arg("mMu"), py::arg("nf"));
    m.def("C2m_g3_approximationB_klmv", &C2m_g3_approximationB_klmv, "", py::arg("x"), py::arg("mQ"), py::arg("mMu"), py::arg("nf"));
    m.def("C2m_g3_approximationBlowxi_klmv", &C2m_g3_approximationBlowxi_klmv, "", py::arg("x"), py::arg("mQ"), py::arg("mMu"), py::arg("nf"));

    m.def("C2m_ps3_approximationA_klmv", &C2m_ps3_approximationA_klmv, "", py::arg("x"), py::arg("mQ"), py::arg("mMu"), py::arg("nf"));
    m.def("C2m_ps3_approximationB_klmv", &C2m_ps3_approximationB_klmv, "", py::arg("x"), py::arg("mQ"), py::arg("mMu"), py::arg("nf"));

}
