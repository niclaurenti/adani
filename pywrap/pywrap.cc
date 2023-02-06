#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/operators.h>
#include <pybind11/functional.h>

#include <adani/adani.h>

namespace py = pybind11;

PYBIND11_MODULE(_core, m) {

    // Documentation

    m.doc() = "Python wrapper of libadani";

    // ApproximateCoefficientFunctions

    m.def("C2m_g3_approximation", &C2m_g3_approximation, py::arg("x"), py::arg("mQ"), py::arg("mMu"), py::arg("nf"), py::arg("method_flag")=1, py::arg("calls")=25000);
    m.def("CLm_g3_approximation", &CLm_g3_approximation, py::arg("x"), py::arg("mQ"), py::arg("mMu"), py::arg("nf"), py::arg("method_flag")=1, py::arg("calls")=25000);
    m.def("C2m_ps3_approximation", &C2m_ps3_approximation, py::arg("x"), py::arg("mQ"), py::arg("mMu"), py::arg("nf"));
    m.def("CLm_ps3_approximation", &CLm_ps3_approximation, py::arg("x"), py::arg("mQ"), py::arg("mMu"), py::arg("nf"));

    m.def("C2m_g2_approximation", &C2m_g2_approximation, py::arg("x"), py::arg("mQ"), py::arg("mMu")) ;
    m.def("CLm_g2_approximation", &CLm_g2_approximation, py::arg("x"), py::arg("mQ"), py::arg("mMu")) ;
    m.def("C2m_ps2_approximation", &C2m_ps2_approximation, py::arg("x"), py::arg("mQ"), py::arg("mMu")) ;
    m.def("CLm_ps2_approximation", &CLm_ps2_approximation, py::arg("x"), py::arg("mQ"), py::arg("mMu")) ;

    m.def("C2m_g3_approximationA_klmv", &C2m_g3_approximationA_klmv, py::arg("x"), py::arg("mQ"), py::arg("mMu"), py::arg("nf"));
    m.def("C2m_g3_approximationB_klmv", &C2m_g3_approximationB_klmv, py::arg("x"), py::arg("mQ"), py::arg("mMu"), py::arg("nf"));
    m.def("C2m_g3_approximationBlowxi_klmv", &C2m_g3_approximationBlowxi_klmv, py::arg("x"), py::arg("mQ"), py::arg("mMu"), py::arg("nf"));

    m.def("C2m_ps3_approximationA_klmv", &C2m_ps3_approximationA_klmv, py::arg("x"), py::arg("mQ"), py::arg("mMu"), py::arg("nf"));
    m.def("C2m_ps3_approximationB_klmv", &C2m_ps3_approximationB_klmv, py::arg("x"), py::arg("mQ"), py::arg("mMu"), py::arg("nf"));

    // AsymptoticCoefficientFunctions

    m.def("C2m_g3_asymptotic", &C2m_g3_asymptotic, py::arg("x"), py::arg("mQ"), py::arg("mMu"), py::arg("nf"), py::arg("v1")=0, py::arg("v2")=0);

    // HighscaleCoefficientFunctions

    m.def("CLm_g3_highscale", &CLm_g3_highscale, py::arg("x"), py::arg("mQ"), py::arg("mMu"), py::arg("nf"));
    m.def("DLm_g3_highscale", &DLm_g3_highscale, py::arg("x"), py::arg("mQ"), py::arg("mMu"), py::arg("nf"));
    m.def("C2m_g3_highscale", &C2m_g3_highscale, py::arg("x"), py::arg("mQ"), py::arg("mMu"), py::arg("nf"), py::arg("v")=0);
    m.def("CLm_ps3_highscale", &CLm_ps3_highscale, py::arg("x"), py::arg("mQ"), py::arg("mMu"), py::arg("nf"));
    m.def("C2m_ps3_highscale", &C2m_ps3_highscale, py::arg("x"), py::arg("mQ"), py::arg("mMu"), py::arg("nf"));
    m.def("D2m_ps3_highscale", &D2m_ps3_highscale, py::arg("x"), py::arg("mQ"), py::arg("mMu"), py::arg("nf"), py::arg("v")=0);

    m.def("C2m_g2_highscale", &C2m_g2_highscale, py::arg("x"), py::arg("mQ"), py::arg("mMu"));
    m.def("CLm_g2_highscale", &CLm_g2_highscale, py::arg("x"), py::arg("mQ"), py::arg("mMu"));
    m.def("C2m_ps2_highscale", &C2m_ps2_highscale, py::arg("x"), py::arg("mQ"), py::arg("mMu"));
    m.def("CLm_ps2_highscale", &CLm_ps2_highscale, py::arg("x"), py::arg("mQ"), py::arg("mMu"));
    m.def("D2m_g2_highscale", &D2m_g2_highscale, py::arg("x"), py::arg("mQ"), py::arg("mMu"));

    // MassiveCoefficientFunctions

    m.def("C2m_g2", &C2m_g2, py::arg("x"), py::arg("mQ"), py::arg("mMu"));
    m.def("C2m_ps2", &C2m_ps2, py::arg("x"), py::arg("mQ"), py::arg("mMu"));
    m.def("CLm_g2", &CLm_g2, py::arg("x"), py::arg("mQ"), py::arg("mMu"));
    m.def("CLm_ps2", &CLm_ps2, py::arg("x"), py::arg("mQ"), py::arg("mMu"));

    m.def("C2m_g21", &C2m_g21, py::arg("x"), py::arg("mQ"));
    m.def("C2m_ps21", &C2m_ps21, py::arg("x"), py::arg("mQ"));
    m.def("CLm_g21", &CLm_g21, py::arg("x"), py::arg("mQ"));
    m.def("CLm_ps21", &CLm_ps21, py::arg("x"), py::arg("mQ"));

    m.def("CLm_g32", &CLm_g32, py::arg("x"), py::arg("mQ"), py::arg("nf"), py::arg("method_flag")=1, py::arg("calls")=25000);
    m.def("C2m_g32", &CLm_g32, py::arg("x"), py::arg("mQ"), py::arg("nf"), py::arg("method_flag")=1, py::arg("calls")=25000);

    // MasslessCoefficientFunctions

    m.def("CL_g3", &CL_g3, py::arg("x"), py::arg("nf"));

    // SplittingFunctions

    m.def("Pgg1reg", &Pgg1reg, py::arg("x"), py::arg("nf"));
    m.def("Pgg1sing", &Pgg1sing, py::arg("x"), py::arg("nf"));
    m.def("Pgg1loc", &Pgg1loc, py::arg("nf"));

    // ThresholdCoefficientFunctions

    m.def("C2m_g3_threshold", &C2m_g3_threshold, py::arg("x"), py::arg("mQ"), py::arg("mMu"), py::arg("nf"));

    // HighEnergyCoefficientFunctions

    m.def("C2m_g3_power_terms", &C2m_g3_power_terms, py::arg("x"), py::arg("mQ"), py::arg("mMu"), py::arg("nf"), py::arg("v")=0);

}
