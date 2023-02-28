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

    double Lmu = log(mMu);

    double c_const = 10./3. + (1. - mQ) * I + (13./6. - 5./3. * mQ) * J ;

    double c_Lmu = 2. + (1. - mQ) * J ;

    return  4. / 3. * CA * (c_const + c_Lmu * Lmu) / x ;

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

    double Lmu = log(1./mMu);

    double c_const = (
        4. * (-1. + 12. * mQ) / (3. + 12. * mQ)
        + (5. - 12. * mQ + 1. / (1. + 4. * mQ)) * J / 6.
        - 4. * mQ * (1. + 3. * mQ) / (1. + 4. * mQ) * I
    ) / 3. ;

    double c_Lmu = (
        - 4. * (1. + 6. * mQ) / (1. + 4. * mQ)
        + 4. * mQ * (1. + 3. * mQ) / (1. + 4. * mQ) * J
    ) / 3. ;

    return 8 * CA * TR * (c_const + c_Lmu * Lmu) / x ;

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

    double Logxi = log(1. + 1. / (4. * mQ));

    double Lmu = log(mMu);
    double Lmu2 = Lmu * Lmu;

    return CA * CA * (
        - 1472. / 27 - 8. / 3 * K * (-1. + mQ)
        + 8. / 27 * J * (-71. + 92. * mQ)
        + I * (
            8. / 3 * Logxi * (-1. + mQ)
            + 8. / 9 * (-13. + 10. * mQ)
        )
        + (
            - 160. / 9 + 16. / 3 * I * (-1. + mQ)
            + 8. / 9 * J * (-13. + 10. * mQ)
        ) * Lmu
        + (- 16. / 3 + 8. / 3 * J * (-1. + mQ)) * Lmu2
    ) * log(x) / x;

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

    double NLL = (
        a21 * (
            160. / 9 - 16. / 3 * II * (-1. + mQ)
            - 8. / 9 * J * (-13. + 10. * mQ)
        )
        + a10 * a11 * (
            2944. / 27 + 16. / 3 * K * (-1. + mQ)
            - 16. / 27 * J * (-71. + 92. * mQ)
            + II * (
                - 16. / 3 * Logxi * (-1. + mQ)
                - 16. / 9 * (-13. + 10. * mQ)
            )
        )
        + a11 * beta0 * (
            - 1472. / 27 - 8. / 3 * K * (-1. + mQ)
            + 8. / 27 * J * (-71. + 92. * mQ)
            + II * (
                8. / 3 * Logxi * (-1. + mQ)
                + 8. / 9 * (-13. + 10. * mQ)
            )
        )
        + (
            a21 * (32. / 3 - 16. / 3 * J * (-1. + mQ))
            + a10 * a11 * (
                320. / 9 - 32. / 3 * II * (-1. + mQ)
                - 16. / 9 * J * (-13. + 10. * mQ)
            )
            + a11 * beta0 * (
                - 160. / 9 + 16. / 3 * II * (-1. + mQ)
                + 8. / 9 * J * (-13. + 10. * mQ)
            )
        ) * Lmu
        + (
            a10 * a11 * (32. / 3 - 16. / 3 * J * (-1. + mQ))
            + a11 * beta0 * (- 16. / 3 + 8. / 3 * J * (-1. + mQ))
        ) * Lmu2
    ) / x;

    return C2_g3_highenergyLL(x, mQ, mMu) + NLL ;

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

double CL_g3_highenergyLL(double x, double mQ, double mMu) {

    if(x>=1 || x<0) return 0;

    double mQ_2 = mQ * mQ ;

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

    double Lmu = log(mMu);
    double Lmu2 = Lmu * Lmu;

    double Logxi = log(1 + 1./(4 * mQ));

    return a11 * a11 * (
        -32. / 3 * K * mQ * (1. + 3. * mQ)
        - 128. / 27 * (17. + 120. * mQ)
        + 16. / 27 * J * (3. + 136. * mQ + 480. * mQ_2)
        + II * (
            32. / 3 * Logxi * mQ * (1. + 3. * mQ)
            + 16. / 9 * (-3. - 4. * mQ + 24. * mQ_2)
        )
        + (
            64. / 3 * II * mQ * (1. + 3. * mQ)
            - 64. /9 * (-1. + 12 * mQ)
            + 16. / 9 * J * (-3. - 4. * mQ + 24. * mQ_2)
        ) * Lmu
        + (
            32. / 3 * J * mQ * (1. + 3. * mQ)
            - 32. / 3 * (1. + 6. * mQ)
        ) * Lmu2
    ) * log(x) / x / (1. + 4. * mQ) ;

}
//==========================================================================================//
//  High energy limit of the gluon coefficient function for FL at O(alpha_s^3).
//------------------------------------------------------------------------------------------//

double CL_g3_highenergy(double x, double mQ, double mMu, int nf) {

    if(x>=1 || x<0) return 0;

    double mQ_2 = mQ * mQ ;

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

    double NLL = (
        a21 * (
            -64. / 3 * II * mQ * (1. + 3. * mQ)
            + 64. / 9 * (-1. + 12. * mQ)
            - 16. / 9 * J * (-3. - 4. * mQ + 24. * mQ_2)
        )
        + a10 * a11 * (
            64. / 3 * K * mQ * (1. + 3. * mQ)
            + 256. / 27 * (17. + 120. * mQ)
            - 32. / 27 * J * (3. + 136. * mQ + 480. * mQ_2)
            + II * (
                -64. / 3 * Logxi * mQ * (1. + 3. * mQ)
                - 32. / 9 * (-3. - 4. * mQ + 24. * mQ_2)
            )
        )
        + a11 * beta0 * (
            - 32. / 3 * K * mQ * (1. + 3. * mQ)
            - 128. / 27 * (17. + 120. * mQ)
            + 16. / 27 * J * (3. + 136. * mQ + 480. * mQ_2)
            + II * (32. / 3 * Logxi * mQ * (1. + 3. * mQ)
                + 16. / 9 * (-3. - 4. * mQ + 24. * mQ_2)
            )
        )
        + (
            a21 * (
                - 64. / 3 * J * mQ * (1. + 3. * mQ)
                + 64. / 3 * (1. + 6. * mQ)
            )
            + a10 * a11 * (
                - 128. / 3 * II * mQ * (1. + 3. * mQ)
                + 128. / 9 * (-1. + 12. * mQ)
                - 32. / 9 * J * (-3. - 4. * mQ + 24. * mQ_2)
            )
            + a11 * beta0 * (
                64. / 3 * II * mQ * (1. + 3. * mQ)
                - 64./9 * (-1. + 12. * mQ)
                + 16. / 9 * J * (-3. - 4. * mQ + 24. * mQ_2)
            )
        ) * Lmu
        + (
            a11 * beta0 * (
                32. / 3 * J * mQ * (1. + 3. * mQ)
                - 32. / 3 * (1. + 6. * mQ)
            )
            + a10 * a11 * (
                - 64. / 3 * J * mQ * (1. + 3. * mQ)
                + 64. / 3 * (1. + 6. * mQ)
            )
        ) * Lmu2
    ) / x / (1. + 4. * mQ) ;

    return CL_g3_highenergyLL(x, mQ, mMu) + NLL ;

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

    double LQ = log(mQ);
    double LQ2 = LQ * LQ;
    double LQ3 = LQ * LQ2;

    double Lmu = log(mMu);
    double Lmu2 = Lmu * Lmu;

    return CA * CA * log(x) * (
        - 32. / 27 * (-71. + 18 * zeta2) * LQ
        - 208. / 9 * LQ2 + 32. / 9 * LQ3
        + Lmu2 * (
            - 16. / 3 + 32. / 3 * LQ
        )
        + Lmu * (
            32. / 9 * (-5. + 6 * zeta2)
            + 416. / 9 * LQ - 32. / 3 * LQ2
        )
        + 16. / 27 * (-92. + 78. * zeta2 - 72. * zeta3)
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

    double a11 = CA ;
    double a21 = nf * (26 * CF - 23 * CA) / 36. ;
    double a10 = - (11. * CA + 2 * nf * (1. - 2. * CF / CA)) / 12. ;

    double beta0 = beta(0, nf) ;

    double NLL = (
        - 32. / 9 * a21 * (-5. + pi2)
        + (
            - 416. * a21 / 9
            + 64. / 27 * a10 * a11 * (-71. + 3. * pi2)
            - 32. / 27 * a11 * beta0 * (-71. + 3. * pi2)
        ) * LQ
        + (
            416. * a10 * a11 / 9
            + 32. * a21 / 3
            - 208. * a11 * beta0 / 9
        ) * LQ2
        + (- 64. * a10 * a11 / 9 + 32. * a11 * beta0 / 9) * LQ3
        + Lmu2 * (
            32. * a10 * a11 / 3 - 16. * a11 * beta0 / 3
            + (- 64 * a10 * a11 / 3 + 32. * a11 * beta0 / 3) * LQ
        )
        + Lmu * (
            32. * a21 / 3
            - 64. / 9 * a10 * a11 * (-5. + pi2)
            + 32. / 9 * a11 * beta0 * (-5. + pi2)
            + (
                - 832. * a10 * a11 / 9 - 64 * a21 / 3
                + 416. * a11 * beta0 / 9
            ) * LQ
            + ( 64. * a10 * a11 / 3 - 32. * a11 * beta0 / 3) * LQ2
        )
        - 32. / 27 * a10 * a11 * (-92. + 13. * pi2 - 72. * zeta3)
        + 16. / 27 * a11 * beta0 * (-92. + 13. * pi2 - 72. * zeta3)
    ) / x;

    return C2_g3_highenergy_highscaleLL(x, mQ, mMu) + NLL ;

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
//  O(alpha_s^3) at leading log.
//------------------------------------------------------------------------------------------//

double CL_g3_highenergy_highscaleLL(double x, double mQ, double mMu) {

    if(x>=1 || x<0) return 0;

    double Lmu = log(mMu);
    double Lmu2 = Lmu * Lmu;

    double LQ = log(mQ);
    double LQ2 = LQ * LQ;

    return CA * CA * (
        32. / 27 * (-68. + 18. * zeta2)
        - 32. / 3 * Lmu2 - 64. / 9 * LQ
        - 32. / 3 * LQ2
        + Lmu * (64. / 9 + 64. / 3 * LQ)
    ) * log(x) / x ;

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

    double a11 = CA ;
    double a21 = nf * (26 * CF - 23 * CA) / 36. ;
    double a10 = - (11. * CA + 2 * nf * (1. - 2. * CF / CA)) / 12. ;

    double beta0 = beta(0, nf) ;

    double NLL = (
        - 64. * a21 / 9 - 64. / 27 * a10 * a11 * (-68. + 3 * pi2)
        + 32. / 27 * a11 * beta0 * (-68. + 3 * pi2)
        + (64. * a10 * a11 / 3 - 32. * a11 * beta0 / 3) * Lmu2
        + (128. * a10 * a11 / 9 - 64. * a21 / 3 - 64. * a11 * beta0 / 9) * LQ
        + (64. * a10 * a11 / 3 - 32. * a11 * beta0 / 3) * LQ2
        + Lmu * (
            - 128. * a10 * a11 / 9 + 64. * a21 / 3 + 64. * a11 * beta0 / 9
            + (- 128. * a10 * a11 / 3 + 64. * a11 * beta0 / 3) * LQ
        )
    ) / x;

    return CL_g3_highenergy_highscaleLL(x, mQ, mMu) + NLL ;

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
