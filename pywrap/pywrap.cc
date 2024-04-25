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

    // CoefficientFunction
    py::class_<CoefficientFunction>(m, "CoefficientFunction")
        .def("GetOrder", &CoefficientFunction::GetOrder)
        .def("GetKind", &CoefficientFunction::GetKind)
        .def("GetChannel", &CoefficientFunction::GetChannel);

    // ApproximateCoefficientFunction
    py::class_<ApproximateCoefficientFunction, CoefficientFunction>(
        m, "ApproximateCoefficientFunction"
    )
        .def(
            py::init<
                const int &, const char &, const char &, const bool &,
                const string &, const double &, const double &, const int &,
                const string &, const int &>(),
            py::arg("order"), py::arg("kind"), py::arg("channel"),
            py::arg("NLL") = true, py::arg("highscale_version") = "klmv",
            py::arg("abserr") = 1e-3, py::arg("relerr") = 1e-3,
            py::arg("dim") = 1000, py::arg("double_int_method") = "analytical",
            py::arg("MCcalls") = 25000
        )
        .def(
            "MuIndependentTerms",
            &ApproximateCoefficientFunction::MuIndependentTerms, py::arg("x"),
            py::arg("m2Q2"), py::arg("nf")
        )
        .def(
            "MuDependentTerms",
            &ApproximateCoefficientFunction::MuDependentTerms, py::arg("x"),
            py::arg("m2Q2"), py::arg("m2mu2"), py::arg("nf")
        )
        .def(
            "fx", &ApproximateCoefficientFunction::fx, py::arg("x"),
            py::arg("m2Q2"), py::arg("m2mu2"), py::arg("nf")
        )
        .def(
            "MuIndependentTermsBand",
            &ApproximateCoefficientFunction::MuIndependentTermsBand,
            py::arg("x"), py::arg("m2Q2"), py::arg("nf")
        )
        .def(
            "MuDependentTermsBand",
            &ApproximateCoefficientFunction::MuDependentTermsBand, py::arg("x"),
            py::arg("m2Q2"), py::arg("m2mu2"), py::arg("nf")
        )
        .def(
            "fxBand", &ApproximateCoefficientFunction::fxBand, py::arg("x"),
            py::arg("m2Q2"), py::arg("m2mu2"), py::arg("nf")
        );

    // ApproximateCoefficientFunctionKLMV
    py::class_<ApproximateCoefficientFunctionKLMV, CoefficientFunction>(
        m, "ApproximateCoefficientFunctionKLMV"
    )
        .def(
            py::init<
                const int &, const char &, const char &, const string &,
                const bool &, const double &, const double &, const int &,
                const string &, const int &>(),
            py::arg("order"), py::arg("kind"), py::arg("channel"),
            py::arg("highscale_version") = "klmv", py::arg("lowxi") = false,
            py::arg("abserr") = 1e-3, py::arg("relerr") = 1e-3,
            py::arg("dim") = 1000, py::arg("double_int_method") = "analytical",
            py::arg("MCcalls") = 25000
        )
        .def(
            "MuIndependentTerms",
            &ApproximateCoefficientFunctionKLMV::MuIndependentTerms,
            py::arg("x"), py::arg("m2Q2"), py::arg("nf")
        )
        .def(
            "MuDependentTerms",
            &ApproximateCoefficientFunctionKLMV::MuDependentTerms, py::arg("x"),
            py::arg("m2Q2"), py::arg("m2mu2"), py::arg("nf")
        )
        .def(
            "fx", &ApproximateCoefficientFunctionKLMV::fx, py::arg("x"),
            py::arg("m2Q2"), py::arg("m2mu2"), py::arg("nf")
        )
        .def(
            "MuIndependentTermsBand",
            &ApproximateCoefficientFunctionKLMV::MuIndependentTermsBand,
            py::arg("x"), py::arg("m2Q2"), py::arg("nf")
        )
        .def(
            "MuDependentTermsBand",
            &ApproximateCoefficientFunctionKLMV::MuDependentTermsBand,
            py::arg("x"), py::arg("m2Q2"), py::arg("m2mu2"), py::arg("nf")
        )
        .def(
            "fxBand", &ApproximateCoefficientFunctionKLMV::fxBand, py::arg("x"),
            py::arg("m2Q2"), py::arg("m2mu2"), py::arg("nf")
        );

    // AsymptoticCoefficientFunction
    py::class_<AsymptoticCoefficientFunction, CoefficientFunction>(
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
            &AsymptoticCoefficientFunction::MuIndependentTerms, py::arg("x"),
            py::arg("m2Q2"), py::arg("nf")
        )
        .def(
            "MuDependentTerms",
            &AsymptoticCoefficientFunction::MuDependentTerms, py::arg("x"),
            py::arg("m2Q2"), py::arg("m2mu2"), py::arg("nf")
        )
        .def(
            "fx", &AsymptoticCoefficientFunction::fx, py::arg("x"),
            py::arg("m2Q2"), py::arg("m2mu2"), py::arg("nf")
        )
        .def(
            "MuIndependentTermsBand",
            &AsymptoticCoefficientFunction::MuIndependentTermsBand,
            py::arg("x"), py::arg("m2Q2"), py::arg("nf")
        )
        .def(
            "MuDependentTermsBand",
            &AsymptoticCoefficientFunction::MuDependentTermsBand, py::arg("x"),
            py::arg("m2Q2"), py::arg("m2mu2"), py::arg("nf")
        )
        .def(
            "fxBand", &AsymptoticCoefficientFunction::fxBand, py::arg("x"),
            py::arg("m2Q2"), py::arg("m2mu2"), py::arg("nf")
        );

    // ExactCoefficientFunction
    py::class_<ExactCoefficientFunction, CoefficientFunction>(
        m, "ExactCoefficientFunction"
    )
        .def(
            py::init<
                const int &, const char &, const char &, const double &,
                const double &, const int &, const string &, const int &>(),
            py::arg("order"), py::arg("kind"), py::arg("channel"),
            py::arg("abserr") = 1e-3, py::arg("relerr") = 1e-3,
            py::arg("dim") = 1000, py::arg("double_int_method") = "analytical",
            py::arg("MCcalls") = 25000
        )
        .def(
            "MuIndependentTerms", &ExactCoefficientFunction::MuIndependentTerms,
            py::arg("x"), py::arg("m2Q2"), py::arg("nf")
        )
        .def(
            "MuDependentTerms", &ExactCoefficientFunction::MuDependentTerms,
            py::arg("x"), py::arg("m2Q2"), py::arg("m2mu2"), py::arg("nf")
        )
        .def(
            "fx", &ExactCoefficientFunction::fx, py::arg("x"), py::arg("m2Q2"),
            py::arg("m2mu2"), py::arg("nf")
        )
        .def(
            "MuIndependentTermsBand",
            &ExactCoefficientFunction::MuIndependentTermsBand, py::arg("x"),
            py::arg("m2Q2"), py::arg("nf")
        )
        .def(
            "MuDependentTermsBand",
            &ExactCoefficientFunction::MuDependentTermsBand, py::arg("x"),
            py::arg("m2Q2"), py::arg("m2mu2"), py::arg("nf")
        )
        .def(
            "fxBand", &ExactCoefficientFunction::fxBand, py::arg("x"),
            py::arg("m2Q2"), py::arg("m2mu2"), py::arg("nf")
        );

    // AbstractHighEnergyCoefficientFunction
    py::class_<AbstractHighEnergyCoefficientFunction, CoefficientFunction>(
        m, "AbstractHighEnergyCoefficientFunction"
    )
        .def("GetNLL", &AbstractHighEnergyCoefficientFunction::GetNLL);

    // HighEnergyCoefficientFunction
    py::class_<
        HighEnergyCoefficientFunction, AbstractHighEnergyCoefficientFunction>(
        m, "HighEnergyCoefficientFunction"
    )
        .def(
            py::init<const int &, const char &, const char &, const bool &>(),
            py::arg("order"), py::arg("kind"), py::arg("channel"),
            py::arg("NLL") = true
        )
        .def(
            "MuIndependentTerms",
            &HighEnergyCoefficientFunction::MuIndependentTerms, py::arg("x"),
            py::arg("m2Q2"), py::arg("nf")
        )
        .def(
            "MuDependentTerms",
            &HighEnergyCoefficientFunction::MuDependentTerms, py::arg("x"),
            py::arg("m2Q2"), py::arg("m2mu2"), py::arg("nf")
        )
        .def(
            "fx", &HighEnergyCoefficientFunction::fx, py::arg("x"),
            py::arg("m2Q2"), py::arg("m2mu2"), py::arg("nf")
        )
        .def(
            "MuIndependentTermsBand",
            &HighEnergyCoefficientFunction::MuIndependentTermsBand,
            py::arg("x"), py::arg("m2Q2"), py::arg("nf")
        )
        .def(
            "MuDependentTermsBand",
            &HighEnergyCoefficientFunction::MuDependentTermsBand, py::arg("x"),
            py::arg("m2Q2"), py::arg("m2mu2"), py::arg("nf")
        )
        .def(
            "fxBand", &HighEnergyCoefficientFunction::fxBand, py::arg("x"),
            py::arg("m2Q2"), py::arg("m2mu2"), py::arg("nf")
        );

    // HighEnergyHighScaleCoefficientFunction
    py::class_<
        HighEnergyHighScaleCoefficientFunction,
        AbstractHighEnergyCoefficientFunction>(
        m, "HighEnergyHighScaleCoefficientFunction"
    )
        .def(
            py::init<const int &, const char &, const char &, const bool &>(),
            py::arg("order"), py::arg("kind"), py::arg("channel"),
            py::arg("NLL") = true
        )
        .def(
            "MuIndependentTerms",
            &HighEnergyHighScaleCoefficientFunction::MuIndependentTerms,
            py::arg("x"), py::arg("m2Q2"), py::arg("nf")
        )
        .def(
            "MuDependentTerms",
            &HighEnergyHighScaleCoefficientFunction::MuDependentTerms,
            py::arg("x"), py::arg("m2Q2"), py::arg("m2mu2"), py::arg("nf")
        )
        .def(
            "fx", &HighEnergyHighScaleCoefficientFunction::fx, py::arg("x"),
            py::arg("m2Q2"), py::arg("m2mu2"), py::arg("nf")
        )
        .def(
            "MuIndependentTermsBand",
            &HighEnergyHighScaleCoefficientFunction::MuIndependentTermsBand,
            py::arg("x"), py::arg("m2Q2"), py::arg("nf")
        )
        .def(
            "MuDependentTermsBand",
            &HighEnergyHighScaleCoefficientFunction::MuDependentTermsBand,
            py::arg("x"), py::arg("m2Q2"), py::arg("m2mu2"), py::arg("nf")
        )
        .def(
            "fxBand", &HighEnergyHighScaleCoefficientFunction::fxBand,
            py::arg("x"), py::arg("m2Q2"), py::arg("m2mu2"), py::arg("nf")
        );

    // PowerTermsCoefficientFunction
    py::class_<
        PowerTermsCoefficientFunction, AbstractHighEnergyCoefficientFunction>(
        m, "PowerTermsCoefficientFunction"
    )
        .def(
            py::init<const int &, const char &, const char &, const bool &>(),
            py::arg("order"), py::arg("kind"), py::arg("channel"),
            py::arg("NLL") = true
        )
        .def(
            "MuIndependentTerms",
            &PowerTermsCoefficientFunction::MuIndependentTerms, py::arg("x"),
            py::arg("m2Q2"), py::arg("nf")
        )
        .def(
            "MuDependentTerms",
            &PowerTermsCoefficientFunction::MuDependentTerms, py::arg("x"),
            py::arg("m2Q2"), py::arg("m2mu2"), py::arg("nf")
        )
        .def(
            "fx", &PowerTermsCoefficientFunction::fx, py::arg("x"),
            py::arg("m2Q2"), py::arg("m2mu2"), py::arg("nf")
        )
        .def(
            "MuIndependentTermsBand",
            &PowerTermsCoefficientFunction::MuIndependentTermsBand,
            py::arg("x"), py::arg("m2Q2"), py::arg("nf")
        )
        .def(
            "MuDependentTermsBand",
            &PowerTermsCoefficientFunction::MuDependentTermsBand, py::arg("x"),
            py::arg("m2Q2"), py::arg("m2mu2"), py::arg("nf")
        )
        .def(
            "fxBand", &PowerTermsCoefficientFunction::fxBand, py::arg("x"),
            py::arg("m2Q2"), py::arg("m2mu2"), py::arg("nf")
        );

    // HighScaleCoefficientFunction
    py::class_<HighScaleCoefficientFunction, CoefficientFunction>(
        m, "HighScaleCoefficientFunction"
    )
        .def(
            py::init<const int &, const char &, const char &, const string &>(),
            py::arg("order"), py::arg("kind"), py::arg("channel"),
            py::arg("version") = "klmv"
        )
        .def(
            "MuIndependentTerms",
            &HighScaleCoefficientFunction::MuIndependentTerms, py::arg("x"),
            py::arg("m2Q2"), py::arg("nf")
        )
        .def(
            "MuDependentTerms", &HighScaleCoefficientFunction::MuDependentTerms,
            py::arg("x"), py::arg("m2Q2"), py::arg("m2mu2"), py::arg("nf")
        )
        .def(
            "fx", &HighScaleCoefficientFunction::fx, py::arg("x"),
            py::arg("m2Q2"), py::arg("m2mu2"), py::arg("nf")
        )
        .def(
            "MuIndependentTermsBand",
            &HighScaleCoefficientFunction::MuIndependentTermsBand, py::arg("x"),
            py::arg("m2Q2"), py::arg("nf")
        )
        .def(
            "MuDependentTermsBand",
            &HighScaleCoefficientFunction::MuDependentTermsBand, py::arg("x"),
            py::arg("m2Q2"), py::arg("m2mu2"), py::arg("nf")
        )
        .def(
            "fxBand", &HighScaleCoefficientFunction::fxBand, py::arg("x"),
            py::arg("m2Q2"), py::arg("m2mu2"), py::arg("nf")
        );

    // HighScaleSplitLogs
    py::class_<HighScaleSplitLogs, CoefficientFunction>(m, "HighScaleSplitLogs")
        .def(
            py::init<const int &, const char &, const char &, const string &>(),
            py::arg("order"), py::arg("kind"), py::arg("channel"),
            py::arg("version") = "klmv"
        )
        .def(
            "fx",
            (double(HighScaleSplitLogs::*)(double, double, int) const)
                & HighScaleSplitLogs::fx,
            py::arg("x"), py::arg("m2Q2"), py::arg("nf")
        )
        .def("LL", &HighScaleSplitLogs::LL, py::arg("x"), py::arg("nf"))
        .def("NLL", &HighScaleSplitLogs::NLL, py::arg("x"), py::arg("nf"))
        .def("N2LL", &HighScaleSplitLogs::N2LL, py::arg("x"), py::arg("nf"))
        .def("N3LL", &HighScaleSplitLogs::N3LL, py::arg("x"), py::arg("nf"))
        .def(
            "fxBand",
            (Value(HighScaleSplitLogs::*)(double, double, int) const)
                & HighScaleSplitLogs::fxBand,
            py::arg("x"), py::arg("m2Q2"), py::arg("nf")
        );

    // MasslessCoefficientFunction
    py::class_<MasslessCoefficientFunction, CoefficientFunction>(
        m, "MasslessCoefficientFunction"
    )
        .def(
            py::init<const int &, const char &, const char &>(),
            py::arg("order"), py::arg("kind"), py::arg("channel")
        )
        .def(
            "MuIndependentTerms",
            (double(MasslessCoefficientFunction::*)(double, int) const)
                & MasslessCoefficientFunction::MuIndependentTerms,
            py::arg("x"), py::arg("nf")
        )
        .def("MuDependentTerms", &MasslessCoefficientFunction::MuDependentTerms)
        .def("fx", &MasslessCoefficientFunction::fx);

    // SplittingFunction
    py::class_<SplittingFunction>(m, "SplittingFunction")
        .def(
            py::init<const int &, const char &, const char &>(),
            py::arg("order"), py::arg("entry1"), py::arg("entry2")
        )
        .def("GetOrder", &SplittingFunction::GetOrder)
        .def("GetEntry1", &SplittingFunction::GetEntry1)
        .def("GetEntry2", &SplittingFunction::GetEntry2)
        .def(
            "Regular", &SplittingFunction::Regular, py::arg("x"), py::arg("nf")
        )
        .def(
            "Singular", &SplittingFunction::Singular, py::arg("x"),
            py::arg("nf")
        )
        .def(
            "SingularIntegrated", &SplittingFunction::SingularIntegrated,
            py::arg("x"), py::arg("nf")
        )
        .def("Local", &SplittingFunction::Local, py::arg("nf"));

    // ConvolutedSplittingFunctions
    py::class_<ConvolutedSplittingFunctions>(m, "ConvolutedSplittingFunctions")
        .def(
            py::init<
                const int &, const char &, const char &, const int &,
                const char &, const char &>(),
            py::arg("order1"), py::arg("entry1"), py::arg("entry2"),
            py::arg("order2"), py::arg("entry3"), py::arg("entry4")
        )
        .def("GetOrder1", &ConvolutedSplittingFunctions::GetOrder1)
        .def("GetEntry1", &ConvolutedSplittingFunctions::GetEntry1)
        .def("GetEntry2", &ConvolutedSplittingFunctions::GetEntry2)
        .def("GetOrder2", &ConvolutedSplittingFunctions::GetOrder2)
        .def("GetEntry3", &ConvolutedSplittingFunctions::GetEntry3)
        .def("GetEntry4", &ConvolutedSplittingFunctions::GetEntry4)
        .def(
            "Regular", &ConvolutedSplittingFunctions::Regular, py::arg("x"),
            py::arg("nf")
        )
        .def(
            "Singular", &ConvolutedSplittingFunctions::Singular, py::arg("x"),
            py::arg("nf")
        )
        .def(
            "SingularIntegrated",
            &ConvolutedSplittingFunctions::SingularIntegrated, py::arg("x"),
            py::arg("nf")
        )
        .def("Local", &ConvolutedSplittingFunctions::Local, py::arg("nf"));

    // ThresholdCoefficientFunction
    py::class_<ThresholdCoefficientFunction, CoefficientFunction>(
        m, "ThresholdCoefficientFunction"
    )
        .def(
            py::init<const int &, const char &, const char &>(),
            py::arg("order"), py::arg("kind"), py::arg("channel")
        )
        .def(
            "MuIndependentTerms",
            &ThresholdCoefficientFunction::MuIndependentTerms, py::arg("x"),
            py::arg("m2Q2"), py::arg("nf")
        )
        .def(
            "MuDependentTerms", &ThresholdCoefficientFunction::MuDependentTerms,
            py::arg("x"), py::arg("m2Q2"), py::arg("m2mu2"), py::arg("nf")
        )
        .def(
            "fx", &ThresholdCoefficientFunction::fx, py::arg("x"),
            py::arg("m2Q2"), py::arg("m2mu2"), py::arg("nf")
        )
        .def(
            "MuIndependentTermsBand",
            &ThresholdCoefficientFunction::MuIndependentTermsBand, py::arg("x"),
            py::arg("m2Q2"), py::arg("nf")
        )
        .def(
            "BetaIndependentTerms",
            &ThresholdCoefficientFunction::BetaIndependentTerms, py::arg("x"),
            py::arg("m2Q2"), py::arg("m2mu2")
        )
        .def(
            "MuDependentTermsBand",
            &ThresholdCoefficientFunction::MuDependentTermsBand, py::arg("x"),
            py::arg("m2Q2"), py::arg("m2mu2"), py::arg("nf")
        )
        .def(
            "fxBand", &ThresholdCoefficientFunction::fxBand, py::arg("x"),
            py::arg("m2Q2"), py::arg("m2mu2"), py::arg("nf")
        );
}
