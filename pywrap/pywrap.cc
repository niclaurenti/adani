#include <pybind11/functional.h>
#include <pybind11/operators.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include <adani/adani.h>

#include <string>

using std::string;

namespace py = pybind11;

PYBIND11_MODULE(_core, m) {

    // Documentation

    m.doc() = "Python wrapper of libadani";

    // Value
    py::class_<Value>(m, "Value")
        .def(
            py::init<const double &, const double &, const double &>(),
            py::arg("central"), py::arg("higher"), py::arg("lower")
        )
        .def(
            py::init<const double &, const double &>(), py::arg("higher"),
            py::arg("lower")
        )
        .def(py::init<const double &>(), py::arg("central"))
        .def(py::self + py::self)
        .def(py::self + double())
        .def(double() + py::self)
        .def(py::self - double())
        .def(py::self * double())
        .def(double() * py::self)
        .def(py::self / double())
        .def(py::self *= double())
        .def(py::self /= double())
        .def("ToVect", &Value::ToVect)
        .def("GetCentral", &Value::GetCentral)
        .def("GetHigher", &Value::GetHigher)
        .def("GetLower", &Value::GetLower);

    // ApproximateCoefficientFunction
    py::class_<ApproximateCoefficientFunction>(
        m, "ApproximateCoefficientFunction"
    )
        .def(
            py::init<
                const int &, const char &, const char &, const bool &,
                const string &, const double &, const double &, const int &,
                const bool &, const int &>(),
            py::arg("order"), py::arg("kind"), py::arg("channel"),
            py::arg("NLL") = true, py::arg("highscale_version") = "klmv",
            py::arg("abserr") = 1e-3, py::arg("relerr") = 1e-3,
            py::arg("dim") = 1000, py::arg("MCintegral") = false,
            py::arg("MCcalls") = 25000
        )
        .def(
            "MuIndependentTerms",
            &ApproximateCoefficientFunction::MuIndependentTerms
        )
        .def(
            "MuDependentTerms",
            &ApproximateCoefficientFunction::MuDependentTerms
        )
        .def("fx", &ApproximateCoefficientFunction::fx)
        .def(
            "MuIndependentTermsBand",
            &ApproximateCoefficientFunction::MuIndependentTermsBand
        )
        .def(
            "MuDependentTermsBand",
            &ApproximateCoefficientFunction::MuDependentTermsBand
        )
        .def("fxBand", &ApproximateCoefficientFunction::fxBand);

    // ApproximateCoefficientFunctionKLMV
    py::class_<ApproximateCoefficientFunctionKLMV>(
        m, "ApproximateCoefficientFunctionKLMV"
    )
        .def(
            py::init<
                const int &, const char &, const char &, const string &,
                const bool &, const double &, const double &, const int &,
                const bool &, const int &>(),
            py::arg("order"), py::arg("kind"), py::arg("channel"),
            py::arg("highscale_version") = "klmv", py::arg("lowxi") = false,
            py::arg("abserr") = 1e-3, py::arg("relerr") = 1e-3,
            py::arg("dim") = 1000, py::arg("MCintegral") = false,
            py::arg("MCcalls") = 25000
        )
        .def(
            "MuIndependentTerms",
            &ApproximateCoefficientFunctionKLMV::MuIndependentTerms
        )
        .def(
            "MuDependentTerms",
            &ApproximateCoefficientFunctionKLMV::MuDependentTerms
        )
        .def("fx", &ApproximateCoefficientFunctionKLMV::fx)
        .def(
            "MuIndependentTermsBand",
            &ApproximateCoefficientFunctionKLMV::MuIndependentTermsBand
        )
        .def(
            "MuDependentTermsBand",
            &ApproximateCoefficientFunctionKLMV::MuDependentTermsBand
        )
        .def("fxBand", &ApproximateCoefficientFunctionKLMV::fxBand);

    // AsymptoticCoefficientFunction
    py::class_<AsymptoticCoefficientFunction>(
        m, "AsymptoticCoefficientFunction"
    )
        .def(
            py::init<
                const int &, const char &, const char &, const bool &,
                const string &>(),
            py::arg("order"), py::arg("kind"), py::arg("channel"),
            py::arg("NLL") = true, py::arg("highscale_version") = "klmv"
        )
        .def(
            "MuIndependentTerms",
            &AsymptoticCoefficientFunction::MuIndependentTerms
        )
        .def(
            "MuDependentTerms", &AsymptoticCoefficientFunction::MuDependentTerms
        )
        .def("fx", &AsymptoticCoefficientFunction::fx)
        .def(
            "MuIndependentTermsBand",
            &AsymptoticCoefficientFunction::MuIndependentTermsBand
        )
        .def(
            "MuDependentTermsBand",
            &AsymptoticCoefficientFunction::MuDependentTermsBand
        )
        .def("fxBand", &AsymptoticCoefficientFunction::fxBand);

    // ExactCoefficientFunction
    py::class_<ExactCoefficientFunction>(m, "ExactCoefficientFunction")
        .def(
            py::init<
                const int &, const char &, const char &, const double &,
                const double &, const int &, const bool &, const int &>(),
            py::arg("order"), py::arg("kind"), py::arg("channel"),
            py::arg("abserr") = 1e-3, py::arg("relerr") = 1e-3,
            py::arg("dim") = 1000, py::arg("MCintegral") = false,
            py::arg("MCcalls") = 25000
        )
        .def(
            "MuIndependentTerms", &ExactCoefficientFunction::MuIndependentTerms
        )
        .def("MuDependentTerms", &ExactCoefficientFunction::MuDependentTerms)
        .def("fx", &ExactCoefficientFunction::fx)
        .def(
            "MuIndependentTermsBand",
            &ExactCoefficientFunction::MuIndependentTermsBand
        )
        .def(
            "MuDependentTermsBand",
            &ExactCoefficientFunction::MuDependentTermsBand
        )
        .def("fxBand", &ExactCoefficientFunction::fxBand);

    // HighEnergyCoefficientFunction
    py::class_<HighEnergyCoefficientFunction>(
        m, "HighEnergyCoefficientFunction"
    )
        .def(
            py::init<const int &, const char &, const char &, const bool &>(),
            py::arg("order"), py::arg("kind"), py::arg("channel"),
            py::arg("NLL") = true
        )
        .def(
            "MuIndependentTerms",
            &HighEnergyCoefficientFunction::MuIndependentTerms
        )
        .def(
            "MuDependentTerms", &HighEnergyCoefficientFunction::MuDependentTerms
        )
        .def("fx", &HighEnergyCoefficientFunction::fx)
        .def(
            "MuIndependentTermsBand",
            &HighEnergyCoefficientFunction::MuIndependentTermsBand
        )
        .def(
            "MuDependentTermsBand",
            &HighEnergyCoefficientFunction::MuDependentTermsBand
        )
        .def("fxBand", &HighEnergyCoefficientFunction::fxBand);

    // HighEnergyHighScaleCoefficientFunction
    py::class_<HighEnergyHighScaleCoefficientFunction>(
        m, "HighEnergyHighScaleCoefficientFunction"
    )
        .def(
            py::init<const int &, const char &, const char &, const bool &>(),
            py::arg("order"), py::arg("kind"), py::arg("channel"),
            py::arg("NLL") = true
        )
        .def(
            "MuIndependentTerms",
            &HighEnergyHighScaleCoefficientFunction::MuIndependentTerms
        )
        .def(
            "MuDependentTerms",
            &HighEnergyHighScaleCoefficientFunction::MuDependentTerms
        )
        .def("fx", &HighEnergyHighScaleCoefficientFunction::fx)
        .def(
            "MuIndependentTermsBand",
            &HighEnergyHighScaleCoefficientFunction::MuIndependentTermsBand
        )
        .def(
            "MuDependentTermsBand",
            &HighEnergyHighScaleCoefficientFunction::MuDependentTermsBand
        )
        .def("fxBand", &HighEnergyHighScaleCoefficientFunction::fxBand);

    // PowerTermsCoefficientFunction
    py::class_<PowerTermsCoefficientFunction>(
        m, "PowerTermsCoefficientFunction"
    )
        .def(
            py::init<const int &, const char &, const char &, const bool &>(),
            py::arg("order"), py::arg("kind"), py::arg("channel"),
            py::arg("NLL") = true
        )
        .def(
            "MuIndependentTerms",
            &PowerTermsCoefficientFunction::MuIndependentTerms
        )
        .def(
            "MuDependentTerms", &PowerTermsCoefficientFunction::MuDependentTerms
        )
        .def("fx", &PowerTermsCoefficientFunction::fx)
        .def(
            "MuIndependentTermsBand",
            &PowerTermsCoefficientFunction::MuIndependentTermsBand
        )
        .def(
            "MuDependentTermsBand",
            &PowerTermsCoefficientFunction::MuDependentTermsBand
        )
        .def("fxBand", &PowerTermsCoefficientFunction::fxBand);

    // HighScaleCoefficientFunction
    py::class_<HighScaleCoefficientFunction>(m, "HighScaleCoefficientFunction")
        .def(
            py::init<const int &, const char &, const char &, const string &>(),
            py::arg("order"), py::arg("kind"), py::arg("channel"),
            py::arg("version") = "klmv"
        )
        .def(
            "MuIndependentTerms",
            &HighScaleCoefficientFunction::MuIndependentTerms
        )
        .def(
            "MuDependentTerms", &HighScaleCoefficientFunction::MuDependentTerms
        )
        .def("fx", &HighScaleCoefficientFunction::fx)
        .def(
            "MuIndependentTermsBand",
            &HighScaleCoefficientFunction::MuIndependentTermsBand
        )
        .def(
            "MuDependentTermsBand",
            &HighScaleCoefficientFunction::MuDependentTermsBand
        )
        .def("fxBand", &HighScaleCoefficientFunction::fxBand);

    // HighScaleSplitLogs
    py::class_<HighScaleSplitLogs>(m, "HighScaleSplitLogs")
        .def(
            py::init<const int &, const char &, const char &, const string &>(),
            py::arg("order"), py::arg("kind"), py::arg("channel"),
            py::arg("version") = "klmv"
        )
        .def(
            "fx", (double(HighScaleSplitLogs::*)(double, double, int) const)
                      & HighScaleSplitLogs::fx
        )
        .def("LL", &HighScaleSplitLogs::LL)
        .def("NLL", &HighScaleSplitLogs::NLL)
        .def("N2LL", &HighScaleSplitLogs::N2LL)
        .def("N3LL", &HighScaleSplitLogs::N3LL)
        .def(
            "fxBand", (Value(HighScaleSplitLogs::*)(double, double, int) const)
                          & HighScaleSplitLogs::fxBand
        );

    // MasslessCoefficientFunction
    py::class_<MasslessCoefficientFunction>(m, "MasslessCoefficientFunction")
        .def(
            py::init<const int &, const char &, const char &>(),
            py::arg("order"), py::arg("kind"), py::arg("channel")
        )
        .def(
            "MuIndependentTerms",
            (double(MasslessCoefficientFunction::*)(double, int) const)
                & MasslessCoefficientFunction::MuIndependentTerms
        )
        .def("MuDependentTerms", &MasslessCoefficientFunction::MuDependentTerms)
        .def("fx", &MasslessCoefficientFunction::fx);

    // SplittingFunction
    py::class_<SplittingFunction>(m, "SplittingFunction")
        .def(
            py::init<const int &, const char &, const char &>(),
            py::arg("order"), py::arg("entry1"), py::arg("entry2")
        )
        .def("Regular", &SplittingFunction::Regular)
        .def("Singular", &SplittingFunction::Singular)
        .def("SingularIntegrated", &SplittingFunction::SingularIntegrated)
        .def("Local", &SplittingFunction::Local);

    // ConvolutedSplittingFunctions
    py::class_<ConvolutedSplittingFunctions>(m, "ConvolutedSplittingFunctions")
        .def(
            py::init<
                const int &, const char &, const char &, const int &,
                const char &, const char &>(),
            py::arg("order1"), py::arg("entry1"), py::arg("entry2"),
            py::arg("order2"), py::arg("entry3"), py::arg("entry4")
        )
        .def("Regular", &ConvolutedSplittingFunctions::Regular)
        .def("Singular", &ConvolutedSplittingFunctions::Singular)
        .def(
            "SingularIntegrated",
            &ConvolutedSplittingFunctions::SingularIntegrated
        )
        .def("Local", &ConvolutedSplittingFunctions::Local);

    // ThresholdCoefficientFunction
    py::class_<ThresholdCoefficientFunction>(m, "ThresholdCoefficientFunction")
        .def(
            py::init<const int &, const char &, const char &>(),
            py::arg("order"), py::arg("kind"), py::arg("channel")
        )
        .def(
            "MuIndependentTerms",
            &ThresholdCoefficientFunction::MuIndependentTerms
        )
        .def(
            "MuDependentTerms", &ThresholdCoefficientFunction::MuDependentTerms
        )
        .def("fx", &ThresholdCoefficientFunction::fx)
        .def(
            "MuIndependentTermsBand",
            &ThresholdCoefficientFunction::MuIndependentTermsBand
        )
        .def(
            "BetaIndependentTerms",
            &ThresholdCoefficientFunction::BetaIndependentTerms
        )
        .def(
            "MuDependentTermsBand",
            &ThresholdCoefficientFunction::MuDependentTermsBand
        )
        .def("fxBand", &ThresholdCoefficientFunction::fxBand);
}
