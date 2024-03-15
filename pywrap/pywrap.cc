#include <pybind11/functional.h>
#include <pybind11/operators.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include <adani/adani.h>

namespace py = pybind11;

PYBIND11_MODULE(_core, m) {

    // Documentation

    m.doc() = "Python wrapper of libadani";

    // ApproximateCoefficientFunction
    py::class_<ApproximateCoefficientFunction>(m, "ApproximateCoefficientFunction")
        .def(py::init<const int&, const char&, const char&, const bool&, const double&, const double&, const int&, const int&, const int&>())
        .def("MuIndependentTerms", &ApproximateCoefficientFunction::MuIndependentTerms)
        .def("MuDependentTerms", &ApproximateCoefficientFunction::MuDependentTerms)
        .def("fx", &ApproximateCoefficientFunction::fx);

    // AsymptoticCoefficientFunction
    py::class_<AsymptoticCoefficientFunction>(m, "AsymptoticCoefficientFunction")
        .def(py::init<const int&, const char&, const char&, const bool&>())
        .def("MuIndependentTerms", &AsymptoticCoefficientFunction::MuIndependentTerms)
        .def("MuDependentTerms", &AsymptoticCoefficientFunction::MuDependentTerms)
        .def("fx", &AsymptoticCoefficientFunction::fx);

    // ExactCoefficientFunctions
    py::class_<ExactCoefficientFunction>(m, "ExactCoefficientFunction")
        .def(py::init<const int&, const char&, const char&, const double&, const double&, const int&, const int&, const int&>())
        .def("MuIndependentTerms", &ExactCoefficientFunction::MuIndependentTerms)
        .def("MuDependentTerms", &ExactCoefficientFunction::MuDependentTerms)
        .def("fx", &ExactCoefficientFunction::fx);

    // HighEnergyCoefficientFunction
    py::class_<HighEnergyCoefficientFunction>(m, "HighEnergyCoefficientFunction")
        .def(py::init<const int&, const char&, const char&, const bool&>())
        .def("MuIndependentTerms", &HighEnergyCoefficientFunction::MuIndependentTerms)
        .def("MuDependentTerms", &HighEnergyCoefficientFunction::MuDependentTerms)
        .def("fx", &HighEnergyCoefficientFunction::fx);

    // HighEnergyHighScaleCoefficientFunction
    py::class_<HighEnergyHighScaleCoefficientFunction>(m, "HighEnergyHighScaleCoefficientFunction")
        .def(py::init<const int&, const char&, const char&, const bool&>())
        .def("MuIndependentTerms", &HighEnergyHighScaleCoefficientFunction::MuIndependentTerms)
        .def("MuDependentTerms", &HighEnergyHighScaleCoefficientFunction::MuDependentTerms)
        .def("fx", &HighEnergyHighScaleCoefficientFunction::fx);

    // PowerTermsCoefficientFunction
    py::class_<PowerTermsCoefficientFunction>(m, "PowerTermsCoefficientFunction")
        .def(py::init<const int&, const char&, const char&, const bool&>())
        .def("MuIndependentTerms", &PowerTermsCoefficientFunction::MuIndependentTerms)
        .def("MuDependentTerms", &PowerTermsCoefficientFunction::MuDependentTerms)
        .def("fx", &PowerTermsCoefficientFunction::fx);

    // HighScaleCoefficientFunction
    py::class_<HighScaleCoefficientFunction>(m, "HighScaleCoefficientFunction")
        .def(py::init<const int&, const char&, const char&>())
        .def("MuIndependentTerms", &HighScaleCoefficientFunction::MuIndependentTerms)
        .def("MuDependentTerms", &HighScaleCoefficientFunction::MuDependentTerms)
        .def("fx", &HighScaleCoefficientFunction::fx);

    // HighScaleSplitLogs
    py::class_<HighScaleSplitLogs>(m, "HighScaleSplitLogs")
        .def(py::init<const int&, const char&, const char&>())
        .def("MuIndependentTerms", (double (HighScaleSplitLogs::*)(double, double, int) const) &HighScaleSplitLogs::MuIndependentTerms)
        .def("MuDependentTerms", &HighScaleSplitLogs::MuDependentTerms)
        .def("fx", (double (HighScaleSplitLogs::*)(double, double, int, int) const) &HighScaleSplitLogs::fx)
        .def("LL", &HighScaleSplitLogs::LL)
        .def("NLL", &HighScaleSplitLogs::NLL)
        .def("N2LL", &HighScaleSplitLogs::N2LL)
        .def("N3LL", &HighScaleSplitLogs::N3LL);

    // MasslessCoefficientFunction
    py::class_<MasslessCoefficientFunction>(m, "MasslessCoefficientFunction")
        .def(py::init<const int&, const char&, const char&>())
        .def("MuIndependentTerms", (double (MasslessCoefficientFunction::*)(double, int) const) &MasslessCoefficientFunction::MuIndependentTerms)
        .def("MuDependentTerms", &MasslessCoefficientFunction::MuDependentTerms)
        .def("fx", &MasslessCoefficientFunction::fx);

    // SplittingFunction
    py::class_<SplittingFunction>(m, "SplittingFunction")
        .def(py::init<const int&, const char&, const char&>())
        .def("Regular", &SplittingFunction::Regular)
        .def("Singular", &SplittingFunction::Singular)
        .def("SingularIntegrated", &SplittingFunction::SingularIntegrated)
        .def("Local", &SplittingFunction::Local);

    // ConvolutedSplittingFunctions
    py::class_<ConvolutedSplittingFunctions>(m, "ConvolutedSplittingFunctions")
        .def(py::init<const int&, const char&, const char&, const int&, const char&, const char&>())
        .def("Regular", &ConvolutedSplittingFunctions::Regular)
        .def("Singular", &ConvolutedSplittingFunctions::Singular)
        .def("SingularIntegrated", &ConvolutedSplittingFunctions::SingularIntegrated)
        .def("Local", &ConvolutedSplittingFunctions::Local);

    // ThresholdCoefficientFunction
    py::class_<ThresholdCoefficientFunction>(m, "ThresholdCoefficientFunction")
        .def(py::init<const int&, const char&, const char&>())
        .def("MuIndependentTerms", &ThresholdCoefficientFunction::MuIndependentTerms)
        .def("MuDependentTerms", &ThresholdCoefficientFunction::MuDependentTerms)
        .def("fx", &ThresholdCoefficientFunction::fx);

    // for tests

    m.def("c2nlog_", (double (*)(double*, double*)) &c2nlog_);
    m.def("clnlog_", (double (*)(double*, double*)) &clnlog_);
    m.def("c2nloq_", (double (*)(double*, double*)) &c2nloq_);
    m.def("clnloq_", (double (*)(double*, double*)) &clnloq_);

    m.def("c2nlobarg_", (double (*)(double*, double*)) &c2nlobarg_);
    m.def("clnlobarg_", (double (*)(double*, double*)) &clnlobarg_);
    m.def("c2nlobarq_", (double (*)(double*, double*)) &c2nlobarq_);
    m.def("clnlobarq_", (double (*)(double*, double*)) &clnlobarq_);

}
