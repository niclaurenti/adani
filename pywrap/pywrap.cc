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
        .def(py::self - py::self)
        .def(py::self * py::self)
        .def(py::self / py::self)
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
        .def("GetChannel", &CoefficientFunction::GetChannel)
        .def(
            "MuIndependentTerms", &CoefficientFunction::MuIndependentTerms,
            py::arg("x"), py::arg("m2Q2"), py::arg("nf")
        )
        .def(
            "MuDependentTerms", &CoefficientFunction::MuDependentTerms,
            py::arg("x"), py::arg("m2Q2"), py::arg("m2mu2"), py::arg("nf")
        )
        .def(
            "fx", &CoefficientFunction::fx, py::arg("x"), py::arg("m2Q2"),
            py::arg("m2mu2"), py::arg("nf")
        )
        .def(
            "MuIndependentTermsBand",
            &CoefficientFunction::MuIndependentTermsBand, py::arg("x"),
            py::arg("m2Q2"), py::arg("nf")
        )
        .def(
            "MuDependentTermsBand", &CoefficientFunction::MuDependentTermsBand,
            py::arg("x"), py::arg("m2Q2"), py::arg("m2mu2"), py::arg("nf")
        )
        .def(
            "fxBand", &CoefficientFunction::fxBand, py::arg("x"),
            py::arg("m2Q2"), py::arg("m2mu2"), py::arg("nf")
        );

    // AbstractApproximate
    py::class_<AbstractApproximate, CoefficientFunction>(
        m, "AbstractApproximate"
    )
        .def(
            py::init<
                const int &, const char &, const char &, const double &,
                const double &, const int &>(),
            py::arg("order"), py::arg("kind"), py::arg("channel"),
            py::arg("abserr") = 1e-3, py::arg("relerr") = 1e-3,
            py::arg("dim") = 1000
        )
        .def(
            "SetDoubleIntegralMethod",
            &AbstractApproximate::SetDoubleIntegralMethod,
            py::arg("double_int_method"), py::arg("abserr") = 1e-3,
            py::arg("relerr") = 1e-3, py::arg("dim") = 1000,
            py::arg("MCcalls") = 25000
        );

    // ApproximateCoefficientFunction
    py::class_<ApproximateCoefficientFunction, AbstractApproximate>(
        m, "ApproximateCoefficientFunction"
    )
        .def(
            py::init<
                const int &, const char &, const char &, const bool &,
                const string &, const double &, const double &, const int &>(),
            py::arg("order"), py::arg("kind"), py::arg("channel"),
            py::arg("NLL") = true, py::arg("highscale_version") = "exact",
            py::arg("abserr") = 1e-3, py::arg("relerr") = 1e-3,
            py::arg("dim") = 1000
        )
        .def(
            "SetLegacyThreshold", &ApproximateCoefficientFunction::SetLegacyThreshold,
            py::arg("legacy_threshold")
        )
        .def(
            "SetLegacyVariation", &ApproximateCoefficientFunction::SetLegacyVariation,
            py::arg("legacy_var")
        )
        .def(
            "SetLegacyPowerTerms", &ApproximateCoefficientFunction::SetLegacyPowerTerms,
            py::arg("legacy_pt")
        );

    // ApproximateCoefficientFunctionKLMV
    py::class_<ApproximateCoefficientFunctionKLMV, AbstractApproximate>(
        m, "ApproximateCoefficientFunctionKLMV"
    )
        .def(
            py::init<
                const int &, const char &, const char &, const string &,
                const bool &, const double &, const double &, const int &>(),
            py::arg("order"), py::arg("kind"), py::arg("channel"),
            py::arg("highscale_version") = "exact", py::arg("lowxi") = false,
            py::arg("abserr") = 1e-3, py::arg("relerr") = 1e-3,
            py::arg("dim") = 1000
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
            py::arg("NLL") = true, py::arg("highscale_version") = "exact"
        )
        .def(
            "SetLegacyPowerTerms", &AsymptoticCoefficientFunction::SetLegacyPowerTerms
        );

    // ExactCoefficientFunction
    py::class_<ExactCoefficientFunction, CoefficientFunction>(
        m, "ExactCoefficientFunction"
    )
        .def(
            py::init<
                const int &, const char &, const char &, const double &,
                const double &, const int &>(),
            py::arg("order"), py::arg("kind"), py::arg("channel"),
            py::arg("abserr") = 1e-3, py::arg("relerr") = 1e-3,
            py::arg("dim") = 1000
        )
        .def(
            "SetDoubleIntegralMethod",
            &ExactCoefficientFunction::SetDoubleIntegralMethod,
            py::arg("double_int_method"), py::arg("abserr") = 1e-3,
            py::arg("relerr") = 1e-3, py::arg("dim") = 1000,
            py::arg("MCcalls") = 25000
        );

    // AbstractHighEnergyCoefficientFunction
    py::class_<AbstractHighEnergyCoefficientFunction, CoefficientFunction>(
        m, "AbstractHighEnergyCoefficientFunction"
    )
        .def("GetNLL", &AbstractHighEnergyCoefficientFunction::GetNLL)
        .def(
            "LL", &AbstractHighEnergyCoefficientFunction::LL,
            py::arg("m2Q2"), py::arg("m2mu2")
        )
        .def(
            "NLL", &AbstractHighEnergyCoefficientFunction::NLL,
            py::arg("m2Q2"), py::arg("m2mu2"), py::arg("nf")
        );

    // HighEnergyCoefficientFunction
    py::class_<
        HighEnergyCoefficientFunction, AbstractHighEnergyCoefficientFunction>(
        m, "HighEnergyCoefficientFunction"
    )
        .def(
            py::init<const int &, const char &, const char &, const bool &>(),
            py::arg("order"), py::arg("kind"), py::arg("channel"),
            py::arg("NLL") = true
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
        );

    // AbstractPowerTerms
    py::class_<AbstractPowerTerms, CoefficientFunction>(
        m, "AbstractPowerTerms"
    );

    // PowerTermsCoefficientFunction
    py::class_<PowerTermsCoefficientFunction, AbstractPowerTerms>(
        m, "PowerTermsCoefficientFunction"
    )
        .def(
            py::init<const int &, const char &, const char &, const bool &>(),
            py::arg("order"), py::arg("kind"), py::arg("channel"),
            py::arg("NLL") = true
        );

    // MultiplicativeAsymptotic
    py::class_<MultiplicativeAsymptotic, AbstractPowerTerms>(
        m, "MultiplicativeAsymptotic"
    )
        .def(
            py::init<const int &, const char &, const char &, const bool &>(),
            py::arg("order"), py::arg("kind"), py::arg("channel"),
            py::arg("NLL") = true
        );

    // HighScaleCoefficientFunction
    py::class_<HighScaleCoefficientFunction, CoefficientFunction>(
        m, "HighScaleCoefficientFunction"
    )
        .def(
            py::init<const int &, const char &, const char &, const string &>(),
            py::arg("order"), py::arg("kind"), py::arg("channel"),
            py::arg("version") = "exact"
        );

    // HighScaleSplitLogs
    py::class_<HighScaleSplitLogs, CoefficientFunction>(m, "HighScaleSplitLogs")
        .def(
            py::init<const int &, const char &, const char &, const string &>(),
            py::arg("order"), py::arg("kind"), py::arg("channel"),
            py::arg("version") = "exact"
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
        );
    ;

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
            "BetaIndependentTerms",
            &ThresholdCoefficientFunction::BetaIndependentTerms, py::arg("x"),
            py::arg("m2Q2"), py::arg("m2mu2")
        )
        .def(
            "SetLegacyThreshold", &ThresholdCoefficientFunction::SetLegacyThreshold,
            py::arg("legacy_threshold")
        );

    // MatchingCondition
    py::class_<MatchingCondition>(m, "MatchingCondition")
        .def(
            py::init<const int &, const char &, const char &, const string &>(),
            py::arg("order"), py::arg("entry1"), py::arg("entry2"),
            py::arg("version") = "exact"
        )
        .def(
            "MuIndependentNfIndependentTerm",
            &MatchingCondition::MuIndependentNfIndependentTerm, py::arg("x")
        )
        .def(
            "MuIndependentNfDependentTerm",
            &MatchingCondition::MuIndependentNfDependentTerm, py::arg("x")
        )
        .def(
            "MuIndependentTerm", &MatchingCondition::MuIndependentTerm,
            py::arg("x"), py::arg("nf")
        );

    // Exceptions
    //  m.def("throw_exception", &throw_exception);
}
