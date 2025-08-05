#include "adani/ThresholdCoefficientFunction.h"
#include "adani/Constants.h"
#include "adani/ExactCoefficientFunction.h"
#include "adani/SpecialFunctions.h"

#include <cmath>

//==========================================================================================//
//  ThresholdCoefficientFunction: constructor
//------------------------------------------------------------------------------------------//

ThresholdCoefficientFunction::ThresholdCoefficientFunction(
    const int &order, const char &kind, const char &channel
)
    : CoefficientFunction(order, kind, channel) {

    if (GetChannel() == 'g') {
        exact_as1_ = new ExactCoefficientFunction(1, GetKind(), GetChannel());
    } else
        exact_as1_ = nullptr;

    try {
        SetFunctions();
        legacy_threshold_ = false;
    } catch (UnexpectedException &e) {
        e.runtime_error();
    }
}

//==========================================================================================//
//  ThresholdCoefficientFunction: destructor
//------------------------------------------------------------------------------------------//

ThresholdCoefficientFunction::~ThresholdCoefficientFunction() {
    delete exact_as1_;
}

//==========================================================================================//
//  ThresholdCoefficientFunction: function that sets the pointer for fx
//------------------------------------------------------------------------------------------//

void ThresholdCoefficientFunction::SetFunctions() {
    switch (GetChannel()) {
        case 'q':
            expansion_beta_ = nullptr;
            expansion_no_beta_ = nullptr;
            fx_ = &ThresholdCoefficientFunction::ZeroFunction;
            break;
        case 'g':
            switch (GetOrder()) {
                case 1:
                    expansion_beta_ = nullptr;
                    expansion_no_beta_ = nullptr;
                    fx_ = &ThresholdCoefficientFunction::Order1;
                    break;
                case 2:
                    switch (GetKind()) {
                        case '2':
                            expansion_beta_ =
                                &ThresholdCoefficientFunction::C2_g2_threshold_expansion;
                            expansion_no_beta_ = &ThresholdCoefficientFunction::
                                                    C2_g2_threshold_expansion_const;
                            break;
                        case 'L':
                            expansion_beta_ =
                                &ThresholdCoefficientFunction::CL_g2_threshold_expansion;
                            expansion_no_beta_ = &ThresholdCoefficientFunction::
                                                    CL_g2_threshold_expansion_const;
                            break;
                        default:
                            throw UnexpectedException(
                                "Unexpected exception!", __PRETTY_FUNCTION__, __LINE__
                            );
                    }
                    fx_ = &ThresholdCoefficientFunction::ModifiedThreshold2;
                    break;
                case 3:
                    switch (GetKind()) {
                        case '2':
                            expansion_beta_ =
                                &ThresholdCoefficientFunction::C2_g3_threshold_expansion;
                            expansion_no_beta_ = &ThresholdCoefficientFunction::
                                                    C2_g3_threshold_expansion_const;
                            break;
                        case 'L':
                            expansion_beta_ =
                                &ThresholdCoefficientFunction::CL_g3_threshold_expansion;
                            expansion_no_beta_ = &ThresholdCoefficientFunction::
                                                    CL_g3_threshold_expansion_const;
                            break;
                        default:
                            throw UnexpectedException(
                                "Unexpected exception!", __PRETTY_FUNCTION__, __LINE__
                            );
                    }
                    fx_ = &ThresholdCoefficientFunction::ModifiedThreshold3;
                    break;
                default:
                    throw UnexpectedException(
                        "Unexpected exception!", __PRETTY_FUNCTION__, __LINE__
                );
            }
            break;
        default:
            throw UnexpectedException(
                "Unexpected exception!", __PRETTY_FUNCTION__, __LINE__
            );
    }

    if (GetKind() == '2') {
        threshold_as1_ = &ThresholdCoefficientFunction::C2_g1_threshold;
    } else if (GetKind() == 'L') {
        threshold_as1_ = &ThresholdCoefficientFunction::CL_g1_threshold;
    } else {
        throw UnexpectedException(
            "Unexpected exception!", __PRETTY_FUNCTION__, __LINE__
        );
    }
}

//==========================================================================================//
//  ThresholdCoefficientFunction: function that sets the legacy behavior for the
//  threshold
//------------------------------------------------------------------------------------------//

void ThresholdCoefficientFunction::SetLegacyThreshold(
    const bool &legacy_threshold
) {
    try {
        if (legacy_threshold == legacy_threshold_) {
            throw NotValidException(
                "Setting legacy threshold identical to its previous value!",
                __PRETTY_FUNCTION__, __LINE__
            );
        }

        legacy_threshold_ = legacy_threshold;

        if (GetOrder() > 1 && GetChannel() == 'g') {
            if (legacy_threshold) {
                fx_ = &ThresholdCoefficientFunction::PlainThreshold;
                if (GetKind() == 'L') {
                    if (GetOrder() == 2) {
                        expansion_beta_ = &ThresholdCoefficientFunction::
                                            C2_g2_threshold_expansion;
                        expansion_no_beta_ = &ThresholdCoefficientFunction::
                                                C2_g2_threshold_expansion_const;
                    } else if (GetOrder() == 3) {
                        expansion_beta_ = &ThresholdCoefficientFunction::
                                            C2_g3_threshold_expansion;
                        expansion_no_beta_ = &ThresholdCoefficientFunction::
                                                C2_g3_threshold_expansion_const;
                    } else {
                        throw UnexpectedException(
                            "Unexpected exception!", __PRETTY_FUNCTION__, __LINE__
                        );
                    }
                }
            } else {
                switch (GetOrder()) {
                    case 2:
                        fx_ = &ThresholdCoefficientFunction::ModifiedThreshold2;
                        break;
                    case 3:
                        fx_ = &ThresholdCoefficientFunction::ModifiedThreshold3;
                        break;
                    default:
                        throw UnexpectedException(
                            "Unexpected exception!", __PRETTY_FUNCTION__, __LINE__
                        );
                }

                if (GetKind() == 'L') {
                    if (GetOrder() == 2) {
                        expansion_beta_ = &ThresholdCoefficientFunction::
                                            CL_g2_threshold_expansion;
                        expansion_no_beta_ = &ThresholdCoefficientFunction::
                                                CL_g2_threshold_expansion_const;
                    } else if (GetOrder() == 3) {
                        expansion_beta_ = &ThresholdCoefficientFunction::
                                            CL_g3_threshold_expansion;
                        expansion_no_beta_ = &ThresholdCoefficientFunction::
                                                CL_g3_threshold_expansion_const;
                    } else {
                        throw UnexpectedException(
                            "Unexpected exception!", __PRETTY_FUNCTION__, __LINE__
                        );
                    }
                }
            }
        }
    } catch (NotValidException &e) {
        e.warning();
    } catch (UnexpectedException &e) {
        e.runtime_error();
    }
}

//==========================================================================================//
//  ThresholdCoefficientFunction: contral value of the full contribution
//------------------------------------------------------------------------------------------//

double ThresholdCoefficientFunction::fx(
    double x, double m2Q2, double m2mu2, int nf
) const {

    return fxBand(x, m2Q2, m2mu2, nf).GetCentral();
}

//==========================================================================================//
//  ThresholdCoefficientFunction: band of the full contribution
//------------------------------------------------------------------------------------------//

Value ThresholdCoefficientFunction::fxBand(
    double x, double m2Q2, double m2mu2, int nf
) const {

    return (this->*fx_)(x, m2Q2, m2mu2, nf);
}

//==========================================================================================//
//  ThresholdCoefficientFunction: first order
//------------------------------------------------------------------------------------------//

Value ThresholdCoefficientFunction::Order1(
    double x, double m2Q2, double /*m2mu2*/, int /*nf*/
) const {

    return Value((this->*threshold_as1_)(x, m2Q2));
}

//==========================================================================================//
//  ThresholdCoefficientFunction:
//------------------------------------------------------------------------------------------//

Value ThresholdCoefficientFunction::PlainThreshold(
    double x, double m2Q2, double m2mu2, int nf
) const {
    double exp = (this->*expansion_beta_)(x, m2Q2, m2mu2, nf)
                 + (this->*expansion_no_beta_)(m2Q2, m2mu2);
    double central = exact_as1_->fx(x, m2Q2, m2mu2, nf) * exp;
    return Value(central);
}

//==========================================================================================//
//  ThresholdCoefficientFunction:
//------------------------------------------------------------------------------------------//

Value ThresholdCoefficientFunction::ModifiedThreshold2(
    double x, double m2Q2, double m2mu2, int nf
) const {
    double exp = (this->*expansion_beta_)(x, m2Q2, m2mu2, nf)
                 + (this->*expansion_no_beta_)(m2Q2, m2mu2);
    double central = exact_as1_->fx(x, m2Q2, m2mu2, nf) * exp;
    double delta = std::abs(central - (this->*threshold_as1_)(x, m2Q2) * exp);

    return Value(central, central + delta, central - delta);
}

//==========================================================================================//
//  ThresholdCoefficientFunction:
//------------------------------------------------------------------------------------------//

Value ThresholdCoefficientFunction::ModifiedThreshold3(
    double x, double m2Q2, double m2mu2, int nf
) const {
    double exp = (this->*expansion_beta_)(x, m2Q2, m2mu2, nf)
                 + (this->*expansion_no_beta_)(m2Q2, m2mu2);
    double central = exact_as1_->fx(x, m2Q2, m2mu2, nf) * exp;
    double delta_prefactor = std::abs(central - (this->*threshold_as1_)(x, m2Q2) * exp);
    double delta_const = std::abs(BetaIndependentTerms(x, m2Q2, m2mu2));

    double delta = sqrt(delta_prefactor * delta_prefactor + delta_const * delta_const);

    return Value(central, central + delta, central - delta);
}

//==========================================================================================//
//  ThresholdCoefficientFunction: central value of the beta-independent terms
//------------------------------------------------------------------------------------------//

double ThresholdCoefficientFunction::BetaIndependentTerms(
    double x, double m2Q2, double m2mu2
) const {
    if (GetChannel() == 'q')
        return 0.;
    // exact_as1 is independent on nf so we call it with nf=0
    return exact_as1_->fx(x, m2Q2, m2mu2, 0)
           * (this->*expansion_no_beta_)(m2Q2, m2mu2);
}

//==========================================================================================//
//  Threshold limit (x->xmax) of the gluon coefficient function for F2 at
//  O(as). In order to pass to klmv normalization multiply
//  m2Q2*M_PI*x
//
//  Eq. (3.15) of Ref. [arXiv:1205.5727]
//------------------------------------------------------------------------------------------//

double
    ThresholdCoefficientFunction::C2_g1_threshold(double x, double m2Q2) const {

    double beta = sqrt(1. - 4. * m2Q2 * x / (1. - x));
    double xi = 1. / m2Q2;

    return xi * TR * beta / (1. + xi / 4.) / x;
}

//==========================================================================================//
//  Threshold limit (x->xmax) of the gluon coefficient function for FL at
//  O(as). In order to pass to klmv normalization multiply
//  m2Q2*M_PI*x
//
//  Eq. (3.15) of Ref. [arXiv:1205.5727]
//------------------------------------------------------------------------------------------//

double
    ThresholdCoefficientFunction::CL_g1_threshold(double x, double m2Q2) const {

    double beta = sqrt(1. - 4. * m2Q2 * x / (1. - x));
    double beta3 = beta * beta * beta;

    double xi = 1. / m2Q2;
    double xi_p_4 = (4. + xi);

    return 64. / 3 * xi * xi * beta3 / (xi_p_4 * xi_p_4 * xi_p_4) / x;
}

//==========================================================================================//
//
//------------------------------------------------------------------------------------------//

double ThresholdCoefficientFunction::C2_g2_threshold_expansion(
    double x, double m2Q2, double m2mu2, int /*nf*/
) const {

    double beta = sqrt(1. - 4. * m2Q2 * x / (1. - x));

    double logb = log(beta);
    double log2b = logb * logb;

    return 16. * CA * log2b + 16. * CA * (3. * ln2 - 5. / 2) * logb
           + (2 * CF - CA) * M_PI * M_PI / beta + 8. * CA * log(m2mu2) * logb;
}

//==========================================================================================//
//
//------------------------------------------------------------------------------------------//

double ThresholdCoefficientFunction::CL_g2_threshold_expansion(
    double x, double m2Q2, double m2mu2, int /*nf*/
) const {

    double beta = sqrt(1. - 4. * m2Q2 * x / (1. - x));

    double logb = log(beta);
    double log2b = logb * logb;

    return 16. * CA * log2b + 16. * CA * (3. * ln2 - 5. / 2 - 2. / 3) * logb
           + (2 * CF - CA) * M_PI * M_PI / beta + 8. * CA * log(m2mu2) * logb;
}

//==========================================================================================//
//
//------------------------------------------------------------------------------------------//

double ThresholdCoefficientFunction::C2_g2_threshold_expansion_const(
    double m2Q2, double m2mu2
) const {

    double xi = 1. / m2Q2;

    return c0(xi) + 36. * CA * ln2 * ln2 - 60 * CA * ln2
           + log(m2mu2) * (8. * CA * ln2 - c0_bar(xi));
}

//==========================================================================================//
//
//------------------------------------------------------------------------------------------//

double ThresholdCoefficientFunction::CL_g2_threshold_expansion_const(
    double m2Q2, double m2mu2
) const {

    double rhoq = -4. * m2Q2;
    double betaq = sqrt(1. - rhoq);
    double chiq = (betaq - 1) / (betaq + 1);

    double muterm =
        -3. / 4 * ln2 + 0.5 + log((1 + chiq) * (1 + chiq) / 2 / chiq) + 1. / 6;

    return 16 * (CA * aL_10_OK(m2Q2) + 2. * CF * aL_10_QED(m2Q2) - CA * log(m2mu2) * muterm);
}

//==========================================================================================//
//
//------------------------------------------------------------------------------------------//

double ThresholdCoefficientFunction::C2_g3_threshold_expansion(
    double x, double m2Q2, double m2mu2, int nf
) const {

    double xi = 1. / m2Q2;
    double beta = sqrt(1. - 4. * m2Q2 * x / (1. - x));

    double Lm = log(m2mu2);
    double Lm2 = Lm * Lm;
    double l = log(beta);
    double l2 = l * l;
    double l3 = l2 * l;
    double l4 = l3 * l;

    double ln2_2 = ln2 * ln2;
    double ln2_3 = ln2_2 * ln2;

    double pi2 = M_PI * M_PI;
    double pi4 = pi2 * pi2;

    double c_log4 = 128. * CA * CA;

    double c_log3 = (768. * ln2 - 6464. / 9.) * CA * CA + 128. / 9. * CA * nf
                    + 128. * CA * CA * Lm;

    double c_log2 =
        (1728. * ln2_2 - 3232. * ln2 - 208. / 3. * pi2 + 15520. / 9.) * CA * CA
        + (64. * ln2 - 640. / 9.) * CA * nf + 16. * CA * c0(xi)
        + 32. * CA * (CF - CA / 2) * pi2 / beta
        - ((-512 * ln2 + 1136. / 3.) * CA * CA - 32. / 3. * CA * nf
           + 16 * CA * c0_bar(xi))
              * Lm
        + 32 * CA * CA * Lm2;

    double c_log_const =
        (1728. * ln2_3 - 4848 * ln2_2 + 15520. / 3. * ln2 - 208 * pi2 * ln2
         + 936 * zeta3 + 608. / 3. * pi2 - 88856. / 27.)
            * CA * CA
        + (96. * ln2_2 - 640. / 3. * ln2 - 16. / 3. * pi2 + 4592. / 27.) * CA
              * nf
        - 32. * CF * (CF - CA / 2) * pi2 + (48. * ln2 - 40.) * CA * c0(xi);

    double c_log_fracbeta =
        ((-92. / 3. + 32. * ln2) * CA + 8. / 3. * nf) * (CF - CA / 2) * pi2;

    double c_log_Lm =
        -((-672. * ln2_2 + 976 * ln2 + 104. / 3. * pi2 - 4160. / 9.) * CA * CA
          + (-32. * ln2 + 320. / 9.) * CA * nf
          + (48. * ln2 - 40.) * CA * c0_bar(xi) - 8. * CA * c0(xi)
          - 16. * CA * (CF - CA / 2) * pi2 / beta);

    double c_log_Lm2 = (64. * ln2 - 44. / 3.) * CA * CA + 8. / 3. * CA * nf
                       - 8. * CA * c0_bar(xi);

    double c_log =
        c_log_const + c_log_fracbeta / beta + c_log_Lm * Lm + c_log_Lm2 * Lm2;

    double c_fracbeta =
        ((8. * ln2_2 - 68. / 3. * ln2 + 8. / 3. * pi2 - 658. / 9.) * CA
         + (8. / 3. * ln2 - 20. / 9.) * nf + 2 * c0(xi)
         + (26. / 3. * CA + 4. / 3. * nf - 2 * c0_bar(xi)) * Lm)
        * (CF - CA / 2) * pi2;

    double c_fracbeta2 = 4. / 3. * (CF - CA / 2) * (CF - CA / 2) * pi4;

    return c_log4 * l4 + c_log3 * l3 + c_log2 * l2 + c_log * l
           + c_fracbeta / beta + c_fracbeta2 / beta / beta;
}

//==========================================================================================//
//
//------------------------------------------------------------------------------------------//

double ThresholdCoefficientFunction::C2_g3_threshold_expansion_const(
    double m2Q2, double m2mu2
) const {

    double xi = 1. / m2Q2;
    double Lm = log(m2mu2);

    double c_const_sqrt = C2_g2_threshold_expansion_const(m2Q2, m2mu2);

    return c_const_sqrt * c_const_sqrt;
}

//==========================================================================================//
//
//------------------------------------------------------------------------------------------//

double ThresholdCoefficientFunction::CL_g3_threshold_expansion(
    double x, double m2Q2, double m2mu2, int nf
) const {
    double beta = sqrt(1. - 4. * m2Q2 * x / (1. - x));

    double rhoq = -4. * m2Q2;
    double betaq = sqrt(1. - rhoq);
    double chiq = (betaq - 1) / (betaq + 1);

    double Lmu = log(m2mu2);
    double Lmu2 = Lmu * Lmu;
    double logb = log(beta);
    double log2b = logb * logb;
    double log3b = log2b * logb;
    double log4b = log3b * logb;

    double pi2 = M_PI * M_PI;
    double pi4 = pi2 * pi2;

    double c_log4 = 128. * CA * CA;

    double c_log3 = (128 * CA * CA * Lmu + 128 * CA * nf / 9 + CA * CA * (-8000./9 + 768 * ln2));

    double c_log2 = (32 * CA * CA * Lmu2 + CA * CA * (21088./9 + 256 * aL_10_OK(m2Q2)
                    - 208 * pi2 / 3 - 2784 * ln2 + 1152 * ln2 * ln2)
                    + CA * (512 * aL_10_QED(m2Q2) * CF + nf * (-256./3 + 64 * ln2))
                    + Lmu * (32 * CA * nf/3 + CA * CA * (-1904./3 + 576 * ln2
                        - 64 * log((1 + chiq) * (1 + chiq)/(2 * chiq)))))
                    + (32 * CA * (-CA/2 + CF) * pi2 / beta) ;

    double c_log = (-32 * CF * CF * pi2
                          + CA * (CF * (-4864 * aL_10_QED(m2Q2) / 3 + 16 * pi2
                            + 1536 * aL_10_QED(m2Q2) * ln2)
                            + nf * (6128./27 - 16 * pi2 / 3 - 256 * ln2 + 96 * ln2*ln2))
                           + Lmu2 * (8 * CA * nf/3 + CA * CA * (-100 + 96 * ln2 - 32 * log((1 + chiq) * (1 + chiq)/(2 * chiq))))
                           + Lmu * (CA * (256 * aL_10_QED(m2Q2) * CF + nf * (-128./3 + 32 * ln2))
                           + CA * CA * (9632./9 + 128 * aL_10_OK(m2Q2) - (104 * pi2) / 3
                                 - 1296 * ln2 + 576 * ln2 * ln2 + 608./3 * log((1 + chiq) * (1 + chiq)/(2 * chiq))
                                 - 192 * ln2 * log((1 + chiq) * (1 + chiq)/(2 * chiq))))
                           + CA * CA * (-112792./27 - 2432 * aL_10_OK(m2Q2)/3
                           + 2240 * pi2/9 + 9536 * ln2 / 3 + 768 * aL_10_OK(m2Q2) * ln2
                           - 208 * pi2 * ln2 - 528 * ln2 * ln2 + 936 * zeta3)) ;
                        + (-CA/2 + CF) * (16 * CA * Lmu * pi2
                        + (8 * nf * pi2)/3
                        + CA * (-92./3 * pi2 + 32 * pi2 * ln2)) / beta;

    double c_fracbeta = (-CA/2 + CF) * (64 * aL_10_QED(m2Q2) * CF * pi2
                        + CA * (-962 * pi2 /9 + 32 * aL_10_OK(m2Q2) * pi2
                            + 8 * pi4/3 + 388./3 * pi2 * ln2 - 64 * pi2 * ln2 * ln2)
                        + nf * (-20 * pi2 / 9 + 4. /9 * pi2 * 6 * ln2)
                        + Lmu * (4 * nf * pi2/3
                            + CA * (-22 * pi2 /3 + 16 * pi2 * ln2 - 8 * pi2 * log(2 + 1/chiq + chiq))));

    double c_fracbeta2 = (CA - 2 * CF) * (CA - 2 * CF) * pi4 / 3.;

    return c_log4 * log4b + c_log3 * log3b + c_log2 * log2b + c_log * logb
           + c_fracbeta / beta + c_fracbeta2 / beta / beta;
}

//==========================================================================================//
//
//------------------------------------------------------------------------------------------//

double ThresholdCoefficientFunction::CL_g3_threshold_expansion_const(
    double m2Q2, double m2mu2
) const {
    double xi = 1. / m2Q2;
    double Lm = log(m2mu2);

    double c_const_sqrt = CL_g2_threshold_expansion_const(m2Q2, m2mu2);

    return c_const_sqrt * c_const_sqrt;
}

//==========================================================================================//
//  Function needed for the threshold limit.
//
//  Eq. (3.10) from Ref. [arXiv:1205.5727]
//------------------------------------------------------------------------------------------//

double ThresholdCoefficientFunction::c0(double xi) const {

    double y = sqrt(1. + 4. / xi);

    double L1 = log(1. + xi / 2);
    double L2 = log(2. + xi / 2);
    double L3 = log(sqrt(xi) * (y - 1.) / 2);

    double xip2 = 2. + xi;
    double xip4 = 4. + xi;

    double Li_2 = Li2(-2. / xip2);
    double pi2 = M_PI * M_PI;

    double c_CA = 50. - pi2 + 12. * L3 / y + 4. * L3 * L3 + L1 * L1 + 6. * L2
                  - 4. * L2 * L2 + 2 * Li_2 + 48. / xip2 - 4. * L2 / xip2
                  + 64. * L2 / xip2 / xip2 - 128. * L2 / (xip2 * xip2 * xip4)
                  - 160. / xip2 / xip4 - 64. * L2 / xip2 / xip4
                  + 128. / (xip2 * xip4 * xip4) - 12. * (4. + zeta2) / xip4
                  - 8. * L3 * L3 / xip4 + 64. / xip4 / xip4;

    double c_CF = -18. - 2. / 3. * pi2 - 24. * L3 / y - 8. * L3 * L3
                  + 2. * L1 * L1 - 6. * L2 + 4. * Li_2 - 48. / xip2
                  + 8. * L2 / xip2 + 360. / xip2 / xip4
                  + 128. * L2 / xip2 / xip4 - 544. / (xip2 * xip4 * xip4)
                  + 48. * L3 * L3 / xip4 - 8. * L1 * L1 / xip4
                  + (44. + 40. * zeta2) / xip4 - 120. * L2 / xip2 / xip2
                  + 256. * L2 / (xip2 * xip2 * xip4) - 16. * Li_2 / xip4
                  - 272. / xip4 / xip4;

    return CA * c_CA + CF * c_CF;
}

//==========================================================================================//
//  Function needed for the threshold limit.
//
//  Eq. (3.11) from Ref. [arXiv:1205.5727]
//------------------------------------------------------------------------------------------//

double ThresholdCoefficientFunction::c0_bar(double xi) const {

    return 4. * CA * (2. + log(1. + xi / 4.)) - 4. / 3. * TR;
}

double ThresholdCoefficientFunction::aL_10_QED(double m2Q2) const {
    double rhoq = -4. * m2Q2;
    double betaq = sqrt(1. - rhoq);
    double chiq = (betaq - 1) / (betaq + 1);
    double rhoq_p_1_2 = (-1 + rhoq) * (-1 + rhoq);
    double log_chi = log(chiq);
    double log_chi_2 = log_chi * log_chi;

    double g1 = (-M_PI * M_PI / 2. + Li2(-(rhoq / (-2 + rhoq)))
                 - (2 * (1 - rhoq) * log_chi) / betaq - (3 * log_chi_2) / 2.
                 + pow(log(rhoq / (-2 + rhoq)), 2) / 2.)
                / 8.;

    return (3 - 2 * rhoq) / (8. * (-2 + rhoq))
           - (g1 * (-1 + 6 * rhoq)) / rhoq_p_1_2
           - (M_PI * M_PI * (-1 + 6 * rhoq)) / (24. * rhoq_p_1_2)
           + ((-6 + rhoq + rhoq * rhoq) * log_chi) / (8. * betaq * (-2 + rhoq))
           - ((-1 + 6 * rhoq) * log_chi_2) / (8. * rhoq_p_1_2)
           + ((3 + 2 * rhoq * (5 + (-5 + rhoq) * rhoq))
              * log(rhoq / (2. * (-1 + rhoq))))
                 / (4. * (-2 + rhoq) * (-2 + rhoq) * (-1 + rhoq));
}

double ThresholdCoefficientFunction::aL_10_OK(double m2Q2) const {
    double rhoq = -4. * m2Q2;
    double betaq = sqrt(1. - rhoq);
    double chiq = (betaq - 1) / (betaq + 1);
    double rhoq_p_1_2 = (-1 + rhoq) * (-1 + rhoq);
    double log_chi = log(chiq);
    double log_chi_2 = log_chi * log_chi;

    double g1 = (-M_PI * M_PI / 2. + Li2(-(rhoq / (-2 + rhoq)))
                 - (2 * (1 - rhoq) * log_chi) / betaq - (3 * log_chi_2) / 2.
                 + pow(log(rhoq / (-2 + rhoq)), 2) / 2.)
                / 8.;

    double g2 = (12.5 - 15 * ln2 + 9 * ln2 * ln2 + log_chi_2
                 - pow(log(rhoq / (2. * (-1 + rhoq))), 2))
                / 4.;

    return 0.6805555555555556 + g2
           - (M_PI * M_PI * (-4 + rhoq) * rhoq) / (24. * rhoq_p_1_2)
           + (g1 * (1 + 2 * rhoq)) / rhoq_p_1_2 - ln2
           - ((-4 + rhoq) * rhoq * log_chi_2) / (8. * rhoq_p_1_2)
           + ((-1 + rhoq * (-3 - 2 * (-3 + rhoq) * rhoq))
              * log(rhoq / (2. * (-1 + rhoq))))
                 / (4. * (-2 + rhoq) * rhoq_p_1_2);
}
