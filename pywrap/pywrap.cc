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

    m.def("C2_g3_approximation", &C2_g3_approximation, py::arg("x"), py::arg("mQ"), py::arg("mMu"), py::arg("nf"), py::arg("v")=0, py::arg("method_flag")=default_method, py::arg("calls")=default_calls);
    m.def("CL_g3_approximation", &CL_g3_approximation, py::arg("x"), py::arg("mQ"), py::arg("mMu"), py::arg("nf"), py::arg("v")=0, py::arg("method_flag")=default_method, py::arg("calls")=default_calls);
    m.def("C2_ps3_approximation", &C2_ps3_approximation, py::arg("x"), py::arg("mQ"), py::arg("mMu"), py::arg("nf"), py::arg("v")=0);
    m.def("CL_ps3_approximation", &CL_ps3_approximation, py::arg("x"), py::arg("mQ"), py::arg("mMu"), py::arg("nf"), py::arg("v")=0);

    m.def("C2_g2_approximation", &C2_g2_approximation, py::arg("x"), py::arg("mQ"), py::arg("mMu"), py::arg("v")=0) ;
    m.def("CL_g2_approximation", &CL_g2_approximation, py::arg("x"), py::arg("mQ"), py::arg("mMu"), py::arg("v")=0) ;
    m.def("C2_ps2_approximation", &C2_ps2_approximation, py::arg("x"), py::arg("mQ"), py::arg("mMu"), py::arg("v")=0) ;
    m.def("CL_ps2_approximation", &CL_ps2_approximation, py::arg("x"), py::arg("mQ"), py::arg("mMu"), py::arg("v")=0) ;

    m.def("C2_g3_approximationA_klmv", &C2_g3_approximationA_klmv, py::arg("x"), py::arg("mQ"), py::arg("mMu"), py::arg("nf"), py::arg("method_flag")=default_method, py::arg("calls")=default_calls);
    m.def("C2_g3_approximationB_klmv", &C2_g3_approximationB_klmv, py::arg("x"), py::arg("mQ"), py::arg("mMu"), py::arg("nf"), py::arg("method_flag")=default_method, py::arg("calls")=default_calls);
    m.def("C2_g3_approximationBlowxi_klmv", &C2_g3_approximationBlowxi_klmv, py::arg("x"), py::arg("mQ"), py::arg("mMu"), py::arg("nf"), py::arg("method_flag")=default_method, py::arg("calls")=default_calls);

    m.def("C2_ps3_approximationA_klmv", &C2_ps3_approximationA_klmv, py::arg("x"), py::arg("mQ"), py::arg("mMu"), py::arg("nf"));
    m.def("C2_ps3_approximationB_klmv", &C2_ps3_approximationB_klmv, py::arg("x"), py::arg("mQ"), py::arg("mMu"), py::arg("nf"));

    // AsymptoticCoefficientFunctions

    m.def("C2_g3_asymptotic", &C2_g3_asymptotic, py::arg("x"), py::arg("mQ"), py::arg("mMu"), py::arg("nf"), py::arg("v1")=0, py::arg("v2")=0);
    m.def("CL_g3_asymptotic", &CL_g3_asymptotic, py::arg("x"), py::arg("mQ"), py::arg("mMu"), py::arg("nf"), py::arg("v"));


    // HighscaleCoefficientFunctions

    m.def("CL_g3_highscale", &CL_g3_highscale, py::arg("x"), py::arg("mQ"), py::arg("mMu"), py::arg("nf"));
    m.def("DL_g3_highscale", &DL_g3_highscale, py::arg("x"), py::arg("mQ"), py::arg("mMu"), py::arg("nf"));
    m.def("C2_g3_highscale", &C2_g3_highscale, py::arg("x"), py::arg("mQ"), py::arg("mMu"), py::arg("nf"), py::arg("v")=0);
    m.def("CL_ps3_highscale", &CL_ps3_highscale, py::arg("x"), py::arg("mQ"), py::arg("mMu"), py::arg("nf"));
    m.def("C2_ps3_highscale", &C2_ps3_highscale, py::arg("x"), py::arg("mQ"), py::arg("mMu"), py::arg("nf"));
    m.def("D2_ps3_highscale", &D2_ps3_highscale, py::arg("x"), py::arg("mQ"), py::arg("mMu"), py::arg("nf"), py::arg("v")=0);

    m.def("C2_g2_highscale", &C2_g2_highscale, py::arg("x"), py::arg("mQ"), py::arg("mMu"));
    m.def("CL_g2_highscale", &CL_g2_highscale, py::arg("x"), py::arg("mQ"), py::arg("mMu"));
    m.def("C2_ps2_highscale", &C2_ps2_highscale, py::arg("x"), py::arg("mQ"), py::arg("mMu"));
    m.def("CL_ps2_highscale", &CL_ps2_highscale, py::arg("x"), py::arg("mQ"), py::arg("mMu"));
    m.def("D2_g2_highscale", &D2_g2_highscale, py::arg("x"), py::arg("mQ"), py::arg("mMu"));

    m.def("C2_g1_highscale", &C2_g1_highscale, py::arg("x"), py::arg("mQ"));

    // ExactCoefficientFunctions

    m.def("C2_g1", &C2_g1, py::arg("x"), py::arg("mQ"));

    m.def("C2_g2", &C2_g2, py::arg("x"), py::arg("mQ"), py::arg("mMu"));
    m.def("C2_ps2", &C2_ps2, py::arg("x"), py::arg("mQ"), py::arg("mMu"));
    m.def("CL_g2", &CL_g2, py::arg("x"), py::arg("mQ"), py::arg("mMu"));
    m.def("CL_ps2", &CL_ps2, py::arg("x"), py::arg("mQ"), py::arg("mMu"));

    m.def("C2_g21", &C2_g21, py::arg("x"), py::arg("mQ"));
    m.def("C2_ps21", &C2_ps21, py::arg("x"), py::arg("mQ"));
    m.def("CL_g21", &CL_g21, py::arg("x"), py::arg("mQ"));
    m.def("CL_ps21", &CL_ps21, py::arg("x"), py::arg("mQ"));

    m.def("CL_g32", &CL_g32, py::arg("x"), py::arg("mQ"), py::arg("nf"), py::arg("method_flag")=default_method, py::arg("calls")=default_calls);
    m.def("C2_g32", &CL_g32, py::arg("x"), py::arg("mQ"), py::arg("nf"), py::arg("method_flag")=default_method, py::arg("calls")=default_calls);

    // MasslessCoefficientFunctions

    m.def("CL_g3_massless", &CL_g3_massless, py::arg("x"), py::arg("nf"));

    // SplittingFunctions

    m.def("Pgg1reg", &Pgg1reg, py::arg("x"), py::arg("nf"));
    m.def("Pgg1sing", &Pgg1sing, py::arg("x"), py::arg("nf"));
    m.def("Pgg1loc", &Pgg1loc, py::arg("nf"));

    // ThresholdCoefficientFunctions

    m.def("C2_g1_threshold", &C2_g1_threshold, py::arg("x"), py::arg("mQ"));

    m.def("C2_g2_threshold", &C2_g2_threshold, py::arg("x"), py::arg("mQ"), py::arg("mMu"));

    m.def("C2_g3_threshold", &C2_g3_threshold, py::arg("x"), py::arg("mQ"), py::arg("mMu"), py::arg("nf"));
    m.def("CL_g3_threshold", &CL_g3_threshold, py::arg("x"), py::arg("mQ"), py::arg("mMu"), py::arg("nf"));

    // HighEnergyCoefficientFunctions

    m.def("C2_g3_power_terms", &C2_g3_power_terms, py::arg("x"), py::arg("mQ"), py::arg("mMu"), py::arg("nf"), py::arg("v"));
    m.def("CL_g2_highenergy", &CL_g2_highenergy, py::arg("x"), py::arg("mQ"), py::arg("mMu"));

    m.def("C2_g3_highenergy", &C2_g3_highenergy, py::arg("x"), py::arg("mQ"), py::arg("mMu"), py::arg("nf"), py::arg("v"));
    m.def("CL_g3_highenergy", &CL_g3_highenergy, py::arg("x"), py::arg("mQ"), py::arg("mMu"), py::arg("nf"), py::arg("v"));
    m.def("C2_g3_highenergyLL", &C2_g3_highenergyLL, py::arg("x"), py::arg("mQ"), py::arg("mMu"));
    m.def("CL_g3_highenergyLL", &CL_g3_highenergyLL, py::arg("x"), py::arg("mQ"), py::arg("mMu"));
    m.def("C2_g3_highenergy_highscale", &C2_g3_highenergy_highscale, py::arg("x"), py::arg("mQ"), py::arg("mMu"), py::arg("nf"), py::arg("v"));
    m.def("CL_g3_highenergy_highscale", &CL_g3_highenergy_highscale, py::arg("x"), py::arg("mQ"), py::arg("mMu"), py::arg("nf"), py::arg("v"));
    m.def("C2_g3_highenergy_highscaleLL", &C2_g3_highenergy_highscaleLL, py::arg("x"), py::arg("mQ"), py::arg("mMu"));
    m.def("CL_g3_highenergy_highscaleLL", &CL_g3_highenergy_highscaleLL, py::arg("x"), py::arg("mQ"), py::arg("mMu"));

}
