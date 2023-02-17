#include "adani/HighEnergyCoefficientFunctions.h"
#include "adani/SpecialFunctions.h"
#include "adani/Constants.h"
#include "adani/SpecialFunctions.h"
#include <cmath>
#include <iostream>

using namespace std;

//==========================================================================================//
//  High energy limit of the gluon coefficient function for F2 at O(alpha_s^2).
//
//  Eq. (3.38) of Ref. [arXiv:1205.5727]
//------------------------------------------------------------------------------------------//

double C2_g2_highenergy(double x, double mQ, double mMu) {

    if(x>=1 || x<0) return 0;

    double z = sqrt(1. / (1. + 4. * mQ));

    double L = log((1. + z) / (1. - z));

    double Hmp = H(z, 1, 1) + H(z, 1, -1) - H(z, -1, 1) - H(z, -1, -1);

    double I = 4. * z * Hmp;
    double J = 4. * z * L ;

    double LMu = log(mMu);

    double c_const = 10./3. + (1. - mQ) * I + (13./6. - 5./3. * mQ) * J ;

    double c_LMu = 2. + (1. - mQ) * J ;

    return  4. / 3. * CA * (c_const + c_LMu * LMu) / x ;

}

//==========================================================================================//
//  High energy limit of the quark coefficient function for F2 at O(alpha_s^2).
//------------------------------------------------------------------------------------------//

double C2_ps2_highenergy(double x, double mQ, double mMu) {

    return CF / CA * C2_g2_highenergy(x, mQ, mMu) ;

}

//==========================================================================================//
//  High energy limit of the gluon coefficient function for FL at O(alpha_s^2).
//
//  Eq. (16, 21) of Ref. [arXiv:hep-ph/9411431]
//------------------------------------------------------------------------------------------//

double CL_g2_highenergy(double x, double mQ, double mMu) {

    if(x>=1 || x<0) return 0;

    double z = sqrt(1. / (1. + 4. * mQ));

    double L = log((1. + z) / (1. - z));

    double Hmp = H(z, 1, 1) + H(z, 1, -1) - H(z, -1, 1) - H(z, -1, -1);

    double I = 4. * z * Hmp;
    double J = 4. * z * L;

    double LMu = log(1./mMu);

    double c_const = (
        4. * (-1. + 12. * mQ) / (3. + 12. * mQ)
        + (5. - 12. * mQ + 1. / (1. + 4. * mQ)) * J / 6.
        - 4. * mQ * (1. + 3. * mQ) / (1. + 4. * mQ) * I
    ) / 3. ;

    double c_LMu = (
        - 4. * (1. + 6. * mQ) / (1. + 4. * mQ)
        + 4. * mQ * (1. + 3. * mQ) / (1. + 4. * mQ) * J
    ) / 3. ;

    return 8 * CA * TR * (c_const + c_LMu * LMu) / x ;

}

//==========================================================================================//
//  High energy limit of the quark coefficient function for FL at O(alpha_s^2).
//------------------------------------------------------------------------------------------//

double CL_ps2_highenergy(double x, double mQ, double mMu) {

    return CF / CA * CL_g2_highenergy(x, mQ, mMu) ;

}

//==========================================================================================//
//  Q^2>>m^2 limit of the high energy limit of the gluon coefficient function for F2 at O(alpha_s^2).
//
//  Eq. (3.40) of Ref. [arXiv:1205.5727]
//------------------------------------------------------------------------------------------//


double C2_g2_highenergy_highscale(double x, double mQ , double mMu) {

    if(x>=1 || x<0) return 0;

    double LQ = log(1. / mQ);
    double L2Q = LQ * LQ;

    double Lm = log(mMu);

    return CA * (
        8. / 3. * L2Q + 104. / 9. * LQ + 40. / 9. - 16. / 3. * zeta2
        + (16. / 3. * LQ + 8. / 3.) * Lm
    ) / x ;

}

//==========================================================================================//
//  Q^2>>m^2 limit of the high energy limit of the quark coefficient function for F2 at O(alpha_s^2).
//------------------------------------------------------------------------------------------//


double C2_ps2_highenergy_highscale(double x, double mQ, double mMu) {

    return CF / CA * C2_g2_highenergy_highscale(x, mQ, mMu) ;

}

//==========================================================================================//
//  Q^2>>m^2 limit of the high energy limit of the gluon coefficient function for FL at O(alpha_s^2).
//
//  Q^2>>m^2 limit of Eq. (16, 21) of Ref. [arXiv:hep-ph/9411431]
//------------------------------------------------------------------------------------------//

double CL_g2_highenergy_highscale(double x, double mQ , double mMu) {

    if(x>=1 || x<0) return 0;

    double LQ = log(1./mQ);

    double Lm = log(1./mMu);

    double c_const = - 2. / 9. * (1. - 3. * LQ);

    double c_log = - 2. / 3.;

    return 16 * CA * TR / x * (c_const + c_log * Lm);

}

//==========================================================================================//
//  Q^2>>m^2 limit of the high energy limit of the quark coefficient function for FL at O(alpha_s^2).
//------------------------------------------------------------------------------------------//

double CL_ps2_highenergy_highscale(double x, double mQ, double mMu) {

    return CF / CA * CL_g2_highenergy_highscale(x, mQ, mMu) ;

}

//==========================================================================================//
//  Power terms in the small x limit of the gluon coefficient function for F2 at O(alpha_s^2).
//------------------------------------------------------------------------------------------//

double C2_g2_power_terms(double x, double mQ , double mMu) {

    if (x<0 || x>=1) return 0;

    return (
        C2_g2_highenergy(x, mQ, mMu)
        - C2_g2_highenergy_highscale(x, mQ, mMu)
    ) ;

}

//==========================================================================================//
//  Power terms in the small x limit of the quark coefficient function for F2 at O(alpha_s^2).
//------------------------------------------------------------------------------------------//

double C2_ps2_power_terms(double x, double mQ , double mMu) {

    return CF / CA * C2_g2_power_terms(x, mQ, mMu);

}

//==========================================================================================//
//  Power terms in the small x limit the gluon coefficient function for FL at O(alpha_s^2).
//------------------------------------------------------------------------------------------//

double CL_g2_power_terms(double x, double mQ , double mMu) {

    if (x<0 || x>=1) return 0;

    return (
        CL_g2_highenergy(x, mQ, mMu)
        - CL_g2_highenergy_highscale(x, mQ, mMu)
    ) ;

}

//==========================================================================================//
//  Power terms in the small x limit the quark coefficient function for FL at O(alpha_s^2).
//------------------------------------------------------------------------------------------//

double CL_ps2_power_terms(double x, double mQ , double mMu) {

    return CF / CA * CL_g2_power_terms(x, mQ, mMu);

}

//==========================================================================================//
//  High energy limit of the gluon coefficient function for F2 at O(alpha_s^3) at leading log.
//
//  Eq. (3.39) of Ref. [arXiv:1205.5727]
//------------------------------------------------------------------------------------------//

double C2_g3_highenergyLL(double x, double mQ, double mMu) {

    if(x>=1 || x<0) return 0;

    double z = sqrt(1. / (1. + 4. * mQ));

    double L = log((1. + z) / (1. - z));

    double Hmp = H(z, 1, 1) + H(z, 1, -1) - H(z, -1, 1) - H(z, -1, -1);

    double Hmpm = (
        H(z, 1, 1, 1)- H(z, 1, 1, -1) + H(z, 1, -1, 1) - H(z, 1, -1, -1)
        - H(z, -1, 1, 1) + H(z, -1, 1, -1) - H(z, -1, -1, 1) + H(z, -1, -1, -1)
    );

    double I = 4 * z * Hmp;
    double J = 4 * z * L;
    double K = 4 * z * Hmpm;

    double LMu = log(mMu);
    double L2Mu = LMu * LMu;

    double c_const = (
        -184./27. - 1./3. * (1. - mQ) * I * log(1. + 1. / 4 / mQ) - 1./9. * (13. - 10. * mQ) * I
        - 1. / 27. * (71. - 92 * mQ) * J + 1./3. * (1. - mQ) * K
    );

    double c_LMu = -20./9. - 2./3. * (1. - mQ) * I - 1./9. * (13. - 10. * mQ) * J;

    double c_L2Mu = -2./3. - 1./3. * (1. - mQ) * J;

    return 8. * CA * CA * log(x) * (c_const + c_LMu * LMu + c_L2Mu * L2Mu) / x;

}

//==========================================================================================//
//  High energy limit of the gluon coefficient function for F2 at O(alpha_s^3).
//
//  Eq. (3.39) of Ref. [arXiv:1205.5727]
//------------------------------------------------------------------------------------------//

double C2_g3_highenergy(double x, double mQ, double mMu, int nf) {

    if(x>=1 || x<0) return 0;

    double z = sqrt(1. / (1. + 4. * mQ));

    double L = log((1. + z) / (1. - z));

    double Hmp = H(z, 1, 1) + H(z, 1, -1) - H(z, -1, 1) - H(z, -1, -1);

    double Hmpm = (
        H(z, 1, 1, 1) - H(z, 1, 1, -1) + H(z, 1, -1, 1) - H(z,1,-1,-1)
        - H(z, -1, 1, 1) + H(z, -1, 1, -1) - H(z, -1, -1, 1) + H(z, -1, -1, -1)
    );

    double II = 4 * z * Hmp;
    double J = 4 * z * L;
    double K = 4 * z * Hmpm;

    double a11 = CA ;
    double a21 = nf * (26. * CF - 23. * CA) / 36 ;
    double a10 = - (11. * CA + 2. * nf * (1. - 2. * CF / CA)) / 12. ;

    double beta0 = beta(0, nf) ;

    double Lmu = log(mMu);
    double Lmu2 = Lmu * Lmu;

    double Logxi = log(1. + 1. / (4. * mQ));

    double tmp = (

        a21 * (-II * (-1 + mQ) + 10./3 + J * (156./72 - 120./72 * mQ))

        + a10 * a11 * (
            184./9 + J * (568./72 - 736./72 * mQ)
            + K * (-1 + mQ)
            + II * (312./72 - 240./72 * mQ + (1 - mQ) * Logxi)
        )

        + a11 * a11 * (
            Lmu2 * (- 1. + J * (-36./72 + 36./72 * mQ))
            + Lmu * (-10./3 + II * (-1 + mQ) + J * (-156./72 + 120./72 * mQ))
            + (
                -92./9 + K * (36./72 - 36./72 * mQ)
                + J * (-284./72 + 368./72 * mQ)
                + II * (
                    -156. /72 + 120./72 * mQ - 36./72 * (1 - mQ) * Logxi
                )
            )
        ) * (log(x) + Logxi)

        + a11 * (
            - 92./9 + K * (36./72 - 36./72 * mQ)
            + J * (-284./72 + 368./72 * mQ)
            + II * (-156./72 + 120./72 * mQ - 36./72 * (1 - mQ) * Logxi)
        ) * beta0

        + Lmu2 * (
            a10 * a11 * (-J * (-1 + mQ) + 2)
            + a11 * (1./2 * J * (-1 + mQ) - 1.) * beta0
        )

        + Lmu * (
            a10 * a11 * (
                -2 * II * (-1 + mQ) + 20./3
                + J * (312./72 - 240./72 * mQ)
            )
            + a21 * (2 + J * (1 - mQ))
            + a11 * (
                II * (-1 + mQ) - 10./3 + J * (-156./72 + 120./72 * mQ)
            ) * beta0
        )
    ) ;

    return 16. / 3 * tmp / x ;

}

//==========================================================================================//
//  High energy limit of the quark coefficient function for F2 at O(alpha_s^3) at leading log.
//------------------------------------------------------------------------------------------//

double C2_ps3_highenergyLL(double x, double mQ, double mMu) {

    return CF / CA * C2_g3_highenergyLL(x, mQ, mMu);

}

//==========================================================================================//
//  High energy limit of the quark coefficient function for F2 at O(alpha_s^3).
//------------------------------------------------------------------------------------------//

double C2_ps3_highenergy(double x, double mQ, double mMu, int nf) {

    return CF / CA * C2_g3_highenergy(x, mQ, mMu, nf);

}
//==========================================================================================//
//  High energy limit of the gluon coefficient function for FL at O(alpha_s^3).
//------------------------------------------------------------------------------------------//

double CL_g3_highenergy(double x, double mQ, double mMu, int nf) {

        if(x>=1 || x<0) return 0;

    double mQ2 = mQ * mQ ;

    double z = sqrt(1. / (1. + 4. * mQ));

    double L = log((1. + z)/(1. - z));

    double Hmp = H(z, 1, 1) + H(z, 1, -1) - H(z, -1, 1) - H(z, -1, -1);

    double Hmpm = (
        H(z, 1, 1, 1) - H(z, 1, 1, -1) + H(z, 1, -1, 1) - H(z, 1, -1, -1)
        - H(z, -1, 1, 1) + H(z, -1, 1, -1) - H(z, -1, -1, 1) + H(z, -1, -1, -1)
    );

    double II = 4 * z * Hmp;
    double J = 4 * z * L ;
    double K = 4 * z * Hmpm;

    double a11 = CA ;
    double a21 = nf * (26 * CF - 23 * CA) / 36 ;
    double a10 = -(11 * CA + 2 * nf * (1 - 2 * CF / CA)) / 12 ;

    double beta0 = beta(0, nf) ;

    double Lmu = log(mMu);
    double Lmu2 = Lmu * Lmu;

    double Logxi = log(1 + 1./(4 * mQ));

    return (
        a21 * (
            - 64./3 * II * mQ * (1. + 3. * mQ)
            + 64./9 * (-1. + 12. * mQ)
            - 16./9 * J * (- 3. - 4. * mQ + 24. * mQ2)
        )
        + a10 * a11 * (
            64./3 * K * mQ * (1. + 3. * mQ) + 256./27 * (17. + 120. * mQ)
            - 32./27 * J * (3. + 136. * mQ + 480. * mQ2)
            + II * (
                - 64./3 * Logxi * mQ * (1. + 3. * mQ)
                - 32./9 * (-3. - 4. * mQ + 24. * mQ2)
            )
        )
        + a11 * (
            - 32./3 * K * mQ * (1. + 3. * mQ) - 128./27 * (17. + 120. * mQ)
            + 16./27 * J * (3. + 136. * mQ + 480. * mQ2)
            + II * (32./3 * Logxi * mQ * (1. + 3. * mQ)
            + 16./9 * (- 3. - 4. * mQ + 24 * mQ2))
        ) * beta0
        + (
            a21 * (- 64./3 * J * mQ * (1. + 3. * mQ) + 64./3 * (1. + 6. * mQ))
            + a10 * a11 * (
                - 128./3 * II * mQ * (1. + 3. * mQ)
                + 128./9 * (- 1. + 12. * mQ)
                - 32./9 * J * (-3. - 4. * mQ + 24. * mQ2)
            )
            + a11 * (
                64./3 * II * mQ * (1. + 3. * mQ) - 64./9 * (-1. + 12. * mQ)
                + 16./9 * J * (-3. - 4. * mQ + 24. * mQ2)
            ) * beta0
        ) * Lmu
        + (
            a10 * a11 * (-64./3 * J * mQ * (1. + 3. * mQ) + 64./3 * (1. + 6. * mQ))
            + a11 * (32./3 * J * mQ * (1. + 3. * mQ) - 32./3 * (1. + 6. * mQ)) * beta0
        ) * Lmu2
        + a11 * a11 * (
            (
                - 32./3 * K * mQ * (1. + 3. * mQ) - 128./27 * (17. + 120. * mQ)
                + 16./27 * J * (3. + 136. * mQ + 480. * mQ2)
                + II * (
                    32./3 * Logxi * mQ * (1. + 3. * mQ)
                    + 16./9 * (-3. - 4. * mQ + 24. * mQ2)
                )
            )
            + (
                64./3 * II * mQ * (1. + 3. * mQ) - 64./9 * (-1. + 12. * mQ)
                + 16./9 * J * (-3. - 4. * mQ + 24. * mQ2)
            ) * Lmu
            + (32./3 * J * mQ * (1. + 3. * mQ) - 32./3 * (1. + 6. * mQ)) * Lmu2
        ) * (log(x) + Logxi)
    ) / ((4. * mQ + 1.) * x) ;

}

//==========================================================================================//
//  High energy limit of the quark coefficient function for FL at O(alpha_s^3).
//------------------------------------------------------------------------------------------//

double CL_ps3_highenergy(double x, double mQ, double mMu, int nf) {

    return CF / CA * CL_g3_highenergy(x, mQ, mMu, nf);

}

//==========================================================================================//
//  High scale limit of the high energy limit of the gluon coefficient function for F2 at
//  O(alpha_s^3) at leading log.
//
//  Eq. (3.41) of Ref. [arXiv:1205.5727]
//------------------------------------------------------------------------------------------//

double C2_g3_highenergy_highscaleLL(double x, double mQ , double mMu) {

    if(x>=1 || x<0) return 0;

    double LQ = log(1. / mQ);
    double L2Q = LQ * LQ;
    double L3Q = LQ * L2Q;

    double LMu = log(mMu);
    double L2Mu = LMu * LMu;

    double c_L3Q = 32./9.;

    double c_L2Q = 208./9.;

    double c_LQ = 2272./27. - 64./3. * zeta2;

    double c_const = 1472./27. - 416./9.*zeta2 + 128./3.*zeta3;

    double c_LMu = 32./3. * L2Q + 416./9. * LQ + 160./9. - 64./3. * zeta2;

    double c_L2Mu = 32./3. * LQ + 16./3.;

    return - CA * CA * log(x) * (
        c_L3Q * L3Q + c_L2Q * L2Q + c_LQ * LQ
        + c_const + c_LMu * LMu + c_L2Mu * L2Mu
    ) / x;

}


//==========================================================================================//
//  High scale limit of the high energy limit of the gluon coefficient function for F2 at
//  O(alpha_s^3).
//------------------------------------------------------------------------------------------//

double C2_g3_highenergy_highscale(double x, double mQ, double mMu, int nf) {

    if(x>=1 || x<0) return 0;

    double Lmu = log(mMu);
    double Lmu2 = Lmu * Lmu;

    double pi2 = M_PI * M_PI ;

    double LQ = log(mQ);
    double LQ2 = LQ * LQ;
    double LQ3 = LQ2 * LQ;
    double LQ4 = LQ3 * LQ;

    double a11 = CA ;
    double a21 = nf * (26 * CF - 23 * CA) / 36. ;
    double a10 = - (11. * CA + 2 * nf * (1. - 2. * CF / CA)) / 12. ;

    double beta0 = beta(0, nf) ;

    return (
        -32./9 * a21 * (-5. + pi2)
        + (
            416. * a10 * a11 / 9. + 32. * a21 / 3
            - 208 * a11 * beta0 / 9
            + 32./27 * a11 * a11 * (-71. + 3 * pi2 + 39 * ln2)
        ) * Lmu2
        + (
            - 64. * a10 * a11 / 9 + 32. * a11 * beta0 / 9
            - 16./9 * a11 * a11 * (-13 + 4 * ln2)
        ) * LQ3
        - 32./9 * a11 * a11 * LQ4
        + Lmu2 * (
            32. * a10 * a11 / 3 - 16. * a11 * beta0 / 3 + 32./3 * a11 * a11 * ln2
            + (
                - 64. * a10 * a11 / 3 + 32. * a11 * beta0 / 3 - 16./3 * a11 * a11 * (-1. + 4. * ln2)
            ) * LQ
            - 32./3 * a11 * a11 * LQ2
        )
        + Lmu * (
            32. * a21 / 3 - 64./9 * a10 * a11 * (-5. + pi2)
            + 32./9 * a11 * beta0 * (-5. + pi2) - 64./9 * a11 * a11 * (-5. + pi2) * ln2
            + (
                - 832. * a10 * a11 / 9 - 64. * a21 / 3 + 416. * a11 * beta0 / 9
                - 32./9 * a11 * a11 * (-5. + pi2 + 26. * ln2)
            ) * LQ
            + (
                64. * a10 * a11 / 3 - 32. * a11 * beta0 / 3
                + 32. / 9 * a11 * a11 * (-13. + 6. * ln2)
            ) * LQ2
            + 32. / 3 * a11 * a11 * LQ3
        )
        + a11 * a11 * log(x) * (
            - 32./27 * (-71. + 3 * pi2) * LQ
            - 208./9 * LQ2 + 32./9 * LQ3
            + Lmu2 * (- 16./3 + 32./3 * LQ)
            + Lmu * (
                32./9 * (-5. + pi2) + 416./9 * LQ - 32./3 * LQ2
            )
            + 16./27 * (-92. + 13. * pi2 - 72. * zeta3)
        )
        - 32./27 * a10 * a11 * (-92. + 13. * pi2 - 72. * zeta3)
        + 16./27 * a11 * beta0 * (-92. + 13. * pi2 - 72. * zeta3)
        - 16./27 * a11*a11 * 2. * ln2 * (-92. + 13. * pi2 - 72. * zeta3)
        + LQ * (
            - 416. * a21 / 9 + 64./27 * a10 * a11 * (-71 + 3 * pi2)
            - 32./27 * a11 * beta0 * (-71. + 3. * pi2)
            + 16./27 * a11 * a11 * (92. - 13. * pi2 - 142. * 2 * ln2 + pi2 * 12. * ln2 + 72 * zeta3)
        )
    ) / x;

}

//==========================================================================================//
//  High scale limit of the high energy limit of the quark coefficient function for F2 at
//  O(alpha_s^3) at leading log.
//------------------------------------------------------------------------------------------//

double C2_ps3_highenergy_highscaleLL(double x, double mQ, double mMu) {

    return CF / CA * C2_g3_highenergy_highscaleLL(x, mQ, mMu);

}

//==========================================================================================//
//  High scale limit of the high energy limit of the quark coefficient function for F2 at
//  O(alpha_s^3).
//------------------------------------------------------------------------------------------//

double C2_ps3_highenergy_highscale(double x, double mQ, double mMu, int nf) {

    return CF / CA * C2_g3_highenergy_highscale(x, mQ, mMu, nf);

}

//==========================================================================================//
//  High scale limit of the high energy limit of the gluon coefficient function for FL at
//  O(alpha_s^3).
//------------------------------------------------------------------------------------------//

double CL_g3_highenergy_highscale(double x, double mQ, double mMu, int nf) {

    if(x>=1 || x<0) return 0;

    double Lmu = log(mMu);
    double Lmu2 = Lmu * Lmu;

    double pi2 = M_PI * M_PI ;

    double LQ = log(mQ);
    double LQ2 = LQ * LQ;
    double LQ3 = LQ2 * LQ;
    double LQ4 = LQ3 * LQ;

    double a11 = CA ;
    double a21 = nf * (26 * CF - 23 * CA) / 36. ;
    double a10 = - (11. * CA + 2 * nf * (1. - 2. * CF / CA)) / 12. ;

    double beta0 = beta(0, nf) ;

    return (
        - 64. * a21 / 9
        - 64. / 27 * a10 * a11 * (-68. + 3. * pi2)
        + 32./27 * a11 * beta0 * (-68. + 3. * pi2)
        - 64./27 * a11 * a11 * (-68. + 3. * pi2) * ln2
        + (
            128 * a10 * a11 / 9 - 64 * a21 / 3
            - 64 * a11 * beta0 / 9 - 32./27 * a11 * a11 * (-68. + 3. * pi2 - 12. * ln2)
        ) * LQ
        + (
            64 * a10 * a11 / 3 - 32 * a11 * beta0 / 3 + 32./9 * a11 * a11 * (2. + 6 * ln2)
        ) * LQ2
        + 32./3 * a11 * a11 * LQ3
        + Lmu2 * (
            64. * a10 * a11 / 3 - 32 * a11 * beta0 / 3
            + 32./3 * a11 * a11 * 2 * ln2 + 32./3 * a11 * a11 * LQ
        )
        + Lmu * (
            - 128 * a10 * a11 / 9 + 64. * a21 / 3
            + 64. * a11 * beta0 / 9 - 64. / 9 * a11 * a11 * 2 * ln2
            + (
                - 128. * a10 * a11 /3 + 64. * a11 * beta0 / 3 - 64./9 * a11 * a11 * (1. + 6 * ln2)
            ) * LQ
            - 64./3 * a11 * a11 * LQ2
        )
        + a11 * a11 * (
            32./27 * (-68. + 3. * pi2)
            - 32./3 * Lmu2 - 64./9 * LQ
            - 32./3 * LQ2
            + Lmu * (
                64. / 9 + 64. / 3 * LQ
            )
        ) * log(x)
    ) / x;

}

//==========================================================================================//
//  High scale limit of the high energy limit of the quark coefficient function for F2 at
//  O(alpha_s^3).
//------------------------------------------------------------------------------------------//

double CL_ps3_highenergy_highscale(double x, double mQ, double mMu, int nf) {

    return CF / CA * CL_g3_highenergy_highscale(x, mQ, mMu, nf);

}

//______________________________________________________________


double C2_g3_highenergy_ERR(double x, double mQ, double mMu, int nf) {

    if(x>=1 || x<0) return 0;

    double z = sqrt(1. / (1. + 4. * mQ));

    double L = log((1. + z) / (1. - z));

    double Hmp = H(z, 1, 1) + H(z, 1, -1) - H(z, -1, 1) - H(z, -1, -1);

    double Hmpm = (
        H(z, 1, 1, 1) - H(z, 1, 1, -1) + H(z, 1, -1, 1) - H(z, 1, -1, -1)
        - H(z, -1, 1, 1) + H(z, -1, 1, -1) - H(z, -1, -1, 1) + H(z, -1, -1, -1)
    );

    double II = 4 * z * Hmp;
    double J = 4 * z * L;
    double K = 4 * z * Hmpm;

    double a11 = CA ;
    //double a21=nf*(26*CF - 23*CA)/36/pi2;
    double a10 = -(11 * CA + 2 * nf * (1. - 2. * CF / CA)) / 12 ;

    double beta0 = beta(0, nf) ;

    double Lmu = log(mMu);
    double Lmu2 = Lmu * Lmu;

    double Logxi = log(1. + 1. / (4. * mQ));

    double tmp = (
        a10 * a11 * (
            184./9 + J * (1136./144 - 1472./144 * mQ) + K * (-1 + mQ)
            + II * (624./144 - 480./144 * mQ - 288./144 * (-1 + mQ) * 0.5 * Logxi)
        )
        + a11 * beta0 * (
            K * (72./144 - 72./144 * mQ)
            + II * (-312./144 + 240./144 * mQ + 576./144 * (-1 + mQ) * ln2 - 72./144 * (1 - mQ) * Logxi - 378./144 * (-1 + mQ) * zeta3)
            + J * (-568./144 + 736./144 * mQ - 1248./144 * ln2 + 960./144 * mQ * ln2 + 819./144 * zeta3 - 630./144 * mQ * zeta3)
            + (- 1472./144 - 1920./144 * ln2 + 1260./144 * zeta3)
        )
        + a11 * a11 * (
            Lmu2 * (1./2 * J * (-1 + mQ) - 1.)
            + Lmu * (II * (-1 + mQ) - 10./3 + J * (-312./144 + 240./144 * mQ))
            + (- 92./9 + K * (72./144 - 72./144 * mQ) + J * (-568./144 + 736./144 * mQ) + II * (-312./144 + 240./144 * mQ + (-1 + mQ) * 0.5 * Logxi))
        ) * (log(x) + Logxi)
        + Lmu2 * (
            a10 * a11 * (-J * (-1 + mQ) + 2)
            + a11 * (1./2 * J * (-1 + mQ) - 1.) * beta0
        )
        + Lmu * (
            a10 * a11 * (-2 * II * (-1 + mQ) + 20./3 + J * (624./144 - 480./144 * mQ))
            + a11 * beta0 * (
                II * (-1 + mQ) + J * (-312./144 + 240./144 * mQ - 576./144 * ln2 + 576./144 * mQ * ln2 + 378./144 * zeta3 - 378./144 * mQ * zeta3)
                + (-480./144 - 1152./144 * ln2 + 756./144 * zeta3)
            )
        )
    );

    return 16. / 3 * tmp / x ;

}

//____________________________________________________________________________________________

double C2_g3_highenergy_highscale_ERR(double x, double mQ, double mMu, int nf) {

    if(x>=1 || x<0) return 0;

    double pi2 = M_PI * M_PI;

    double Lmu = log(mMu);
    double Lmu2 = Lmu * Lmu;

    double LQ = log(mQ);
    double LQ2 = LQ * LQ;
    double LQ3 = LQ2 * LQ;

    double a11 = CA ;
    //double a21= nf*(26*CF - 23*CA)/36/pi2;
    double a10 = -(11 * CA + 2 * nf * (1. - 2. * CF / CA)) / 12 ;

    double beta0 = beta(0, nf) ;

    double tmp = (
        Lmu2 * (
            a11 * a10 * (1./6 - 1./3 * LQ)
            + a11 * a11 * (1./12 * (2 * ln2 + LQ) - 1./6 * LQ * (2 * ln2 + LQ))
            + a11 * (-1./12 + LQ/6) * beta0
        )
        + Lmu * (
            a11 * a10 * (-1./9 * (-5 + pi2) - 13./9 * LQ + 1./3 * LQ2)
            + a11 * a11 * (-1./18 * (-5 + pi2) * (2 * ln2 + LQ) - 13./18 * LQ * (2 * ln2 + LQ) + 1./6 * LQ2 * (2 * ln2 + LQ))
            + a11 * beta0 * (-1./6 * LQ2 + 1./432 * LQ * (312 + 576 * ln2 - 378 * zeta3) + 1./432 * (24 * (-5 + pi2 - 12 * ln2) + 189 * zeta3))
        )
        + a11 * a10 * (
            1./432 * (-1136 + 48 * pi2) * LQ
            + 13./18 * LQ2 - 1./9 * LQ3
            + 1./432 * (736 - 104 * pi2 + 576 * zeta3)
        )
        + a11 * beta0 * (
            LQ3/18 + 1./432 * LQ * (568 - 24 * pi2 + 1248 * ln2 - 819 * zeta3)
            + 1./432 * LQ2 * (-156 - 288 * ln2 + 189 * zeta3)
            + 1./432 * (-368 + 52 * pi2 + 96 * (-5 + pi2) * ln2 + 27 * zeta3 - 63 * pi2 * zeta3)
        )
        + a11 * a11 * (
            13./36 * LQ2 * (2 * ln2 + LQ)
            - 1./18 * LQ3 * (2 * ln2 + LQ)
            + 1./432 * LQ * (-568 * (2 * ln2 + LQ) + 24 * pi2 * (2 * ln2 + LQ))
            + 1./432 * (368 * (2 * ln2 + LQ) - 52 * pi2 * (2 * ln2 + LQ) + 288 * (2 * ln2 + LQ) * zeta3)
        )
        + log(x) * a11 * a11 * (
            Lmu2 * (-1./12 + LQ/6)
            + Lmu * (1./18 * (-5 + pi2) + (13 * LQ)/18 - LQ2/6)
            + (1./432 * (568 - 24 * pi2) * LQ - 13./36 * LQ2 + LQ3/18 + 1./432 * (-368 + 52 * pi2 - 288 * zeta3))
        )
    );

    return 64 * tmp / x;

}

//==========================================================================================//
//  Power terms in the small x limit of the gluon coefficient function for F2 at O(alpha_s^3) at leading log.
//------------------------------------------------------------------------------------------//

double C2_g3_power_termsLL(double x, double mQ , double mMu) {

    if (x<0 || x>=1) return 0;

    return (
        C2_g3_highenergyLL(x, mQ, mMu)
        - C2_g3_highenergy_highscaleLL(x, mQ, mMu)
    ) ;

}

//==========================================================================================//
//  Power terms in the small x limit of the gluon coefficient function for F2 at O(alpha_s^3).
//------------------------------------------------------------------------------------------//

double C2_g3_power_terms(double x, double mQ , double mMu, int nf, int v) {

    if (x<0 || x>=1) return 0;

    return (
        C2_g3_highenergy_BAND(x, mQ, mMu, nf, v)
        - C2_g3_highenergy_highscale_BAND(x, mQ, mMu, nf, v)
    ) ;

}

//==========================================================================================//
//  Power terms in the small x limit of the quark coefficient function for F2 at O(alpha_s^2)
//  at leading log.
//------------------------------------------------------------------------------------------//

double C2_ps3_power_termsLL(double x, double mQ , double mMu) {

    return CF / CA * C2_g3_power_termsLL(x, mQ, mMu);

}

//==========================================================================================//
//  Power terms in the small x limit of the quark coefficient function for F2 at O(alpha_s^3).
//------------------------------------------------------------------------------------------//

double C2_ps3_power_terms(double x, double mQ , double mMu, int nf) {

    return CF / CA * C2_g3_power_terms(x, mQ, mMu, nf, 0);

}

//==========================================================================================//
//  Power terms in the small x limit of the gluon coefficient function for FL at O(alpha_s^3).
//------------------------------------------------------------------------------------------//

double CL_g3_power_terms(double x, double mQ , double mMu, int nf) {

    if (x<0 || x>=1) return 0;

    return (
        CL_g3_highenergy(x, mQ, mMu, nf)
        - CL_g3_highenergy_highscale(x, mQ, mMu, nf)
    ) ;

}

//==========================================================================================//
//  Power terms in the small x limit of the gluon coefficient function for FL at O(alpha_s^3).
//------------------------------------------------------------------------------------------//

double CL_ps3_power_terms(double x, double mQ , double mMu, int nf) {

    return CF / CA * CL_g3_power_terms(x, mQ, mMu, nf);

}

//////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////


double C2_g3_highenergy_BAND(double x, double mQ, double mMu, int nf, int v) {

    double D=C2_g3_highenergy(x,mQ,mMu,nf);

    if(v==0) return D;

    double delta = fabs(D - C2_g3_highenergy_ERR(x,mQ,mMu,nf));

    if(v==1) return D + delta;
    if(v==2) return D - delta;

    else {
        cout<<"Choose either v=0 or v=1 or v=2!!\nExiting!!\n"<<endl;
        exit(-1);
    }

}

//________________________________________________________________

double C2_g3_highenergy_highscale_BAND(double x, double mQ, double mMu, int nf, int v) {

    double D=C2_g3_highenergy_highscale(x,mQ,mMu,nf);

    if(v==0) return D;

    double delta = fabs(D - C2_g3_highenergy_highscale_ERR(x,mQ,mMu,nf));

    if(v==1) return D + delta;
    if(v==2) return D - delta;

    else {
        cout<<"Choose either v=0 or v=1 or v=2!!\nExiting!!\n"<<endl;
        exit(-1);
    }

}

//________________________________________________________________
