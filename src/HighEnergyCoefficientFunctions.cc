#include "adani/HighEnergyCoefficientFunctions.h"
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

double C2_g2_highenergy(double x, double m2Q2, double m2mu2) {

    double z = sqrt(1. / (1. + 4. * m2Q2));

    double L = log((1. + z) / (1. - z));

    double Hmp = H(z, 1, 1) + H(z, 1, -1) - H(z, -1, 1) - H(z, -1, -1);

    double I = 4. * z * Hmp;
    double J = 4. * z * L;

    double Lmu = log(m2mu2);

    double c_const =
        10. / 3. + (1. - m2Q2) * I + (13. / 6. - 5. / 3. * m2Q2) * J;

    double c_Lmu = 2. + (1. - m2Q2) * J;

    return 4. / 3. * CA * (c_const + c_Lmu * Lmu) / x;
}

//==========================================================================================//
//  High energy limit of the quark coefficient function for F2 at O(alpha_s^2).
//------------------------------------------------------------------------------------------//

double C2_ps2_highenergy(double x, double m2Q2, double m2mu2) {

    return CF / CA * C2_g2_highenergy(x, m2Q2, m2mu2);
}

//==========================================================================================//
//  High energy limit of the gluon coefficient function for FL at O(alpha_s^2).
//
//  Eq. (16, 21) of Ref. [arXiv:hep-ph/9411431]
//------------------------------------------------------------------------------------------//

double CL_g2_highenergy(double x, double m2Q2, double m2mu2) {

    double z = sqrt(1. / (1. + 4. * m2Q2));

    double L = log((1. + z) / (1. - z));

    double Hmp = H(z, 1, 1) + H(z, 1, -1) - H(z, -1, 1) - H(z, -1, -1);

    double I = 4. * z * Hmp;
    double J = 4. * z * L;

    double Lmu = log(1. / m2mu2);

    double c_const = (4. * (-1. + 12. * m2Q2) / (3. + 12. * m2Q2)
                      + (5. - 12. * m2Q2 + 1. / (1. + 4. * m2Q2)) * J / 6.
                      - 4. * m2Q2 * (1. + 3. * m2Q2) / (1. + 4. * m2Q2) * I)
                     / 3.;

    double c_Lmu = (-4. * (1. + 6. * m2Q2) / (1. + 4. * m2Q2)
                    + 4. * m2Q2 * (1. + 3. * m2Q2) / (1. + 4. * m2Q2) * J)
                   / 3.;

    return 8 * CA * TR * (c_const + c_Lmu * Lmu) / x;
}

//==========================================================================================//
//  High energy limit of the quark coefficient function for FL at O(alpha_s^2).
//------------------------------------------------------------------------------------------//

double CL_ps2_highenergy(double x, double m2Q2, double m2mu2) {

    return CF / CA * CL_g2_highenergy(x, m2Q2, m2mu2);
}

//==========================================================================================//
//  Q^2>>m^2 limit of the high energy limit of the gluon coefficient function
//  for F2 at O(alpha_s^2).
//
//  Eq. (3.40) of Ref. [arXiv:1205.5727]
//------------------------------------------------------------------------------------------//

double C2_g2_highenergy_highscale(double x, double m2Q2, double m2mu2) {

    double LQ = log(1. / m2Q2);
    double L2Q = LQ * LQ;

    double Lm = log(m2mu2);

    return CA
           * (8. / 3. * L2Q + 104. / 9. * LQ + 40. / 9. - 16. / 3. * zeta2
              + (16. / 3. * LQ + 8. / 3.) * Lm)
           / x;
}

//==========================================================================================//
//  Q^2>>m^2 limit of the high energy limit of the quark coefficient function
//  for F2 at O(alpha_s^2).
//------------------------------------------------------------------------------------------//

double C2_ps2_highenergy_highscale(double x, double m2Q2, double m2mu2) {

    return CF / CA * C2_g2_highenergy_highscale(x, m2Q2, m2mu2);
}

//==========================================================================================//
//  Q^2>>m^2 limit of the high energy limit of the gluon coefficient function
//  for FL at O(alpha_s^2).
//
//  Q^2>>m^2 limit of Eq. (16, 21) of Ref. [arXiv:hep-ph/9411431]
//------------------------------------------------------------------------------------------//

double CL_g2_highenergy_highscale(double x, double m2Q2, double m2mu2) {

    double LQ = log(1. / m2Q2);

    double Lm = log(1. / m2mu2);

    double c_const = -2. / 9. * (1. - 3. * LQ);

    double c_log = -2. / 3.;

    return 16 * CA * TR / x * (c_const + c_log * Lm);
}

//==========================================================================================//
//  Q^2>>m^2 limit of the high energy limit of the quark coefficient function
//  for FL at O(alpha_s^2).
//------------------------------------------------------------------------------------------//

double CL_ps2_highenergy_highscale(double x, double m2Q2, double m2mu2) {

    return CF / CA * CL_g2_highenergy_highscale(x, m2Q2, m2mu2);
}

//==========================================================================================//
//  Power terms in the small x limit of the gluon coefficient function for F2 at
//  O(alpha_s^2).
//------------------------------------------------------------------------------------------//

double C2_g2_power_terms(double x, double m2Q2, double m2mu2) {

    return C2_g2_highenergy(x, m2Q2, m2mu2)
           - C2_g2_highenergy_highscale(x, m2Q2, m2mu2);
}

//==========================================================================================//
//  Power terms in the small x limit of the quark coefficient function for F2 at
//  O(alpha_s^2).
//------------------------------------------------------------------------------------------//

double C2_ps2_power_terms(double x, double m2Q2, double m2mu2) {

    return CF / CA * C2_g2_power_terms(x, m2Q2, m2mu2);
}

//==========================================================================================//
//  Power terms in the small x limit the gluon coefficient function for FL at
//  O(alpha_s^2).
//------------------------------------------------------------------------------------------//

double CL_g2_power_terms(double x, double m2Q2, double m2mu2) {

    return CL_g2_highenergy(x, m2Q2, m2mu2)
           - CL_g2_highenergy_highscale(x, m2Q2, m2mu2);
}

//==========================================================================================//
//  Power terms in the small x limit the quark coefficient function for FL at
//  O(alpha_s^2).
//------------------------------------------------------------------------------------------//

double CL_ps2_power_terms(double x, double m2Q2, double m2mu2) {

    return CF / CA * CL_g2_power_terms(x, m2Q2, m2mu2);
}

//==========================================================================================//
//  High energy limit of the gluon coefficient function for F2 at O(alpha_s^3)
//  at leading log.
//
//  Eq. (3.39) of Ref. [arXiv:1205.5727]
//------------------------------------------------------------------------------------------//

double C2_g3_highenergyLL(double x, double m2Q2, double m2mu2) {

    double z = sqrt(1. / (1. + 4. * m2Q2));

    double L = log((1. + z) / (1. - z));

    double Hmp = H(z, 1, 1) + H(z, 1, -1) - H(z, -1, 1) - H(z, -1, -1);

    double Hmpm = H(z, 1, 1, 1) - H(z, 1, 1, -1) + H(z, 1, -1, 1)
                  - H(z, 1, -1, -1) - H(z, -1, 1, 1) + H(z, -1, 1, -1)
                  - H(z, -1, -1, 1) + H(z, -1, -1, -1);

    double I = 4 * z * Hmp;
    double J = 4 * z * L;
    double K = 4 * z * Hmpm;

    double Logxi = log(1. + 1. / (4. * m2Q2));

    double Lmu = log(m2mu2);
    double Lmu2 = Lmu * Lmu;

    return CA * CA
           * (-1472. / 27 - 8. / 3 * K * (-1. + m2Q2)
              + 8. / 27 * J * (-71. + 92. * m2Q2)
              + I
                    * (8. / 3 * Logxi * (-1. + m2Q2)
                       + 8. / 9 * (-13. + 10. * m2Q2))
              + (-160. / 9 + 16. / 3 * I * (-1. + m2Q2)
                 + 8. / 9 * J * (-13. + 10. * m2Q2))
                    * Lmu
              + (-16. / 3 + 8. / 3 * J * (-1. + m2Q2)) * Lmu2)
           * log(x) / x;
}

//==========================================================================================//
//  High energy limit of the gluon coefficient function for F2 at O(alpha_s^3)
//  at next to leading log.
//
//  Eq. (3.39) of Ref. [arXiv:1205.5727]
//------------------------------------------------------------------------------------------//

double C2_g3_highenergyNLL(double x, double m2Q2, double m2mu2, int nf, int v) {

    double z = sqrt(1. / (1. + 4. * m2Q2));

    double L = log((1. + z) / (1. - z));

    double Hmp = H(z, 1, 1) + H(z, 1, -1) - H(z, -1, 1) - H(z, -1, -1);

    double Hmpm = H(z, 1, 1, 1) - H(z, 1, 1, -1) + H(z, 1, -1, 1)
                  - H(z, 1, -1, -1) - H(z, -1, 1, 1) + H(z, -1, 1, -1)
                  - H(z, -1, -1, 1) + H(z, -1, -1, -1);

    double II = 4 * z * Hmp;
    double J = 4 * z * L;
    double K = 4 * z * Hmpm;

    double a11 = a_11();
    double a21 = a_21(nf);
    double a10 = a_10(nf);

    double beta0 = beta(0, nf);

    double Lmu = log(m2mu2);
    double Lmu2 = Lmu * Lmu;

    double Logxi = log(1. + 1. / (4. * m2Q2));

    if (v == 0) {

        return (a21
                    * (160. / 9 - 16. / 3 * II * (-1. + m2Q2)
                       - 8. / 9 * J * (-13. + 10. * m2Q2))
                + a10 * a11
                      * (2944. / 27 + 16. / 3 * K * (-1. + m2Q2)
                         - 16. / 27 * J * (-71. + 92. * m2Q2)
                         + II
                               * (-16. / 3 * Logxi * (-1. + m2Q2)
                                  - 16. / 9 * (-13. + 10. * m2Q2)))
                + a11 * beta0
                      * (-1472. / 27 - 8. / 3 * K * (-1. + m2Q2)
                         + 8. / 27 * J * (-71. + 92. * m2Q2)
                         + II
                               * (8. / 3 * Logxi * (-1. + m2Q2)
                                  + 8. / 9 * (-13. + 10. * m2Q2)))
                + (a21 * (32. / 3 - 16. / 3 * J * (-1. + m2Q2))
                   + a10 * a11
                         * (320. / 9 - 32. / 3 * II * (-1. + m2Q2)
                            - 16. / 9 * J * (-13. + 10. * m2Q2))
                   + a11 * beta0
                         * (-160. / 9 + 16. / 3 * II * (-1. + m2Q2)
                            + 8. / 9 * J * (-13. + 10. * m2Q2)))
                      * Lmu
                + (a10 * a11 * (32. / 3 - 16. / 3 * J * (-1. + m2Q2))
                   + a11 * beta0 * (-16. / 3 + 8. / 3 * J * (-1. + m2Q2)))
                      * Lmu2)
               / x;
    }

    double central_value = C2_g3_highenergyNLL(x, m2Q2, m2mu2, nf, 0);

    double error = (a10 * a11
                        * (2944. / 27 + 16. / 3 * K * (-1. + m2Q2)
                           - 16. / 27 * J * (-71 + 92 * m2Q2)
                           + II
                                 * (-16. / 3 * Logxi * (-1. + m2Q2)
                                    - 16. / 9 * (-13 + 10 * m2Q2)))
                    + (a10 * a11 * (32. / 3 - 16. / 3 * J * (-1. + m2Q2))
                       + a11 * beta0 * (-16. / 3 + 8. / 3 * J * (-1 + m2Q2)))
                          * Lmu2
                    + a11 * beta0
                          * (-1472. / 27 - 8. / 3 * K * (-1 + m2Q2)
                             - 640. * ln2 / 9 + 140. * zeta3 / 3
                             + II
                                   * (8. / 3 * Logxi * (-1. + m2Q2)
                                      + 8. / 9 * (-13 + 10 * m2Q2)
                                      + 64. / 3 * (-1 + m2Q2) * ln2
                                      - 14. * (-1 + m2Q2) * zeta3)
                             + J
                                   * (8. / 27 * (-71 + 92 * m2Q2)
                                      + 32. / 9 * (-13 + 10 * m2Q2) * ln2
                                      - 7. / 3 * (-13 + 10 * m2Q2) * zeta3))
                    + Lmu
                          * (a10 * a11
                                 * (320. / 9 - 32. / 3 * II * (-1 + m2Q2)
                                    - 16. / 9 * J * (-13 + 10 * m2Q2))
                             + a11 * beta0
                                   * (-160. / 9 + 16. / 3 * II * (-1 + m2Q2)
                                      - 128. * ln2 / 3 + 28 * zeta3
                                      + J
                                            * (8. / 9 * (-13 + 10 * m2Q2)
                                               + 64. / 3 * (-1 + m2Q2) * ln2
                                               - 14 * (-1 + m2Q2) * zeta3))))
                   / x;

    double delta = fabs(central_value - error);

    if (v == 1)
        return central_value + delta;
    if (v == -1)
        return central_value - delta;
    else {
        std::cout << "C2_g3_highenergyNLL: Choose either v=0, v=1 or "
                     "v=-1!!\nExiting!!\n"
                  << std::endl;
        exit(-1);
    }
}

//==========================================================================================//
//  High energy limit of the gluon coefficient function for F2 at O(alpha_s^3).
//
//  Eq. (3.39) of Ref. [arXiv:1205.5727]
//------------------------------------------------------------------------------------------//

double C2_g3_highenergy(double x, double m2Q2, double m2mu2, int nf, int v) {

    return C2_g3_highenergyLL(x, m2Q2, m2mu2)
           + C2_g3_highenergyNLL(x, m2Q2, m2mu2, nf, v);
}

//==========================================================================================//
//  High energy limit of the quark coefficient function for F2 at O(alpha_s^3)
//  at leading log.
//------------------------------------------------------------------------------------------//

double C2_ps3_highenergyLL(double x, double m2Q2, double m2mu2) {

    return CF / CA * C2_g3_highenergyLL(x, m2Q2, m2mu2);
}

//==========================================================================================//
//  High energy limit of the quark coefficient function for F2 at O(alpha_s^3)
//  at next to leading log.
//------------------------------------------------------------------------------------------//

double
C2_ps3_highenergyNLL(double x, double m2Q2, double m2mu2, int nf, int v) {

    return CF / CA * C2_g3_highenergyNLL(x, m2Q2, m2mu2, nf, v);
}

//==========================================================================================//
//  High energy limit of the quark coefficient function for F2 at O(alpha_s^3).
//------------------------------------------------------------------------------------------//

double C2_ps3_highenergy(double x, double m2Q2, double m2mu2, int nf, int v) {

    return CF / CA * C2_g3_highenergy(x, m2Q2, m2mu2, nf, v);
}

//==========================================================================================//
//  High energy limit of the gluon coefficient function for FL at O(alpha_s^3)
//  at leading log.
//------------------------------------------------------------------------------------------//

double CL_g3_highenergyLL(double x, double m2Q2, double m2mu2) {

    double m4Q4 = m2Q2 * m2Q2;

    double z = sqrt(1. / (1. + 4. * m2Q2));

    double L = log((1. + z) / (1. - z));

    double Hmp = H(z, 1, 1) + H(z, 1, -1) - H(z, -1, 1) - H(z, -1, -1);

    double Hmpm = H(z, 1, 1, 1) - H(z, 1, 1, -1) + H(z, 1, -1, 1)
                  - H(z, 1, -1, -1) - H(z, -1, 1, 1) + H(z, -1, 1, -1)
                  - H(z, -1, -1, 1) + H(z, -1, -1, -1);

    double II = 4 * z * Hmp;
    double J = 4 * z * L;
    double K = 4 * z * Hmpm;

    double a11 = CA;

    double Lmu = log(m2mu2);
    double Lmu2 = Lmu * Lmu;

    double Logxi = log(1 + 1. / (4 * m2Q2));

    return a11 * a11
           * (-32. / 3 * K * m2Q2 * (1. + 3. * m2Q2)
              - 128. / 27 * (17. + 120. * m2Q2)
              + 16. / 27 * J * (3. + 136. * m2Q2 + 480. * m4Q4)
              + II
                    * (32. / 3 * Logxi * m2Q2 * (1. + 3. * m2Q2)
                       + 16. / 9 * (-3. - 4. * m2Q2 + 24. * m4Q4))
              + (64. / 3 * II * m2Q2 * (1. + 3. * m2Q2)
                 - 64. / 9 * (-1. + 12 * m2Q2)
                 + 16. / 9 * J * (-3. - 4. * m2Q2 + 24. * m4Q4))
                    * Lmu
              + (32. / 3 * J * m2Q2 * (1. + 3. * m2Q2)
                 - 32. / 3 * (1. + 6. * m2Q2))
                    * Lmu2)
           * log(x) / x / (1. + 4. * m2Q2);
}

//==========================================================================================//
//  High energy limit of the gluon coefficient function for FL at O(alpha_s^3)
//  at next to leading log.
//------------------------------------------------------------------------------------------//

double CL_g3_highenergyNLL(double x, double m2Q2, double m2mu2, int nf, int v) {

    double m4Q4 = m2Q2 * m2Q2;

    double z = sqrt(1. / (1. + 4. * m2Q2));

    double L = log((1. + z) / (1. - z));

    double Hmp = H(z, 1, 1) + H(z, 1, -1) - H(z, -1, 1) - H(z, -1, -1);

    double Hmpm = H(z, 1, 1, 1) - H(z, 1, 1, -1) + H(z, 1, -1, 1)
                  - H(z, 1, -1, -1) - H(z, -1, 1, 1) + H(z, -1, 1, -1)
                  - H(z, -1, -1, 1) + H(z, -1, -1, -1);

    double II = 4 * z * Hmp;
    double J = 4 * z * L;
    double K = 4 * z * Hmpm;

    double a11 = a_11();
    double a21 = a_21(nf);
    double a10 = a_10(nf);

    double beta0 = beta(0, nf);

    double Lmu = log(m2mu2);
    double Lmu2 = Lmu * Lmu;

    double Logxi = log(1. + 1. / (4. * m2Q2));

    if (v == 0) {

        return (a21
                    * (-64. / 3 * II * m2Q2 * (1. + 3. * m2Q2)
                       + 64. / 9 * (-1. + 12. * m2Q2)
                       - 16. / 9 * J * (-3. - 4. * m2Q2 + 24. * m4Q4))
                + a10 * a11
                      * (64. / 3 * K * m2Q2 * (1. + 3. * m2Q2)
                         + 256. / 27 * (17. + 120. * m2Q2)
                         - 32. / 27 * J * (3. + 136. * m2Q2 + 480. * m4Q4)
                         + II
                               * (-64. / 3 * Logxi * m2Q2 * (1. + 3. * m2Q2)
                                  - 32. / 9 * (-3. - 4. * m2Q2 + 24. * m4Q4)))
                + a11 * beta0
                      * (-32. / 3 * K * m2Q2 * (1. + 3. * m2Q2)
                         - 128. / 27 * (17. + 120. * m2Q2)
                         + 16. / 27 * J * (3. + 136. * m2Q2 + 480. * m4Q4)
                         + II
                               * (32. / 3 * Logxi * m2Q2 * (1. + 3. * m2Q2)
                                  + 16. / 9 * (-3. - 4. * m2Q2 + 24. * m4Q4)))
                + (a21
                       * (-64. / 3 * J * m2Q2 * (1. + 3. * m2Q2)
                          + 64. / 3 * (1. + 6. * m2Q2))
                   + a10 * a11
                         * (-128. / 3 * II * m2Q2 * (1. + 3. * m2Q2)
                            + 128. / 9 * (-1. + 12. * m2Q2)
                            - 32. / 9 * J * (-3. - 4. * m2Q2 + 24. * m4Q4))
                   + a11 * beta0
                         * (64. / 3 * II * m2Q2 * (1. + 3. * m2Q2)
                            - 64. / 9 * (-1. + 12. * m2Q2)
                            + 16. / 9 * J * (-3. - 4. * m2Q2 + 24. * m4Q4)))
                      * Lmu
                + (a11 * beta0
                       * (32. / 3 * J * m2Q2 * (1. + 3. * m2Q2)
                          - 32. / 3 * (1. + 6. * m2Q2))
                   + a10 * a11
                         * (-64. / 3 * J * m2Q2 * (1. + 3. * m2Q2)
                            + 64. / 3 * (1. + 6. * m2Q2)))
                      * Lmu2)
               / x / (1. + 4. * m2Q2);
    }

    double central_value = CL_g3_highenergyNLL(x, m2Q2, m2mu2, nf, 0);

    double error =
        (a10 * a11
             * (64. / 3 * K * m2Q2 * (1 + 3 * m2Q2)
                + 256. / 27 * (17 + 120 * m2Q2)
                - 32. / 27 * J * (3 + 136 * m2Q2 + 480 * m4Q4)
                + II
                      * (-64. / 3 * Logxi * m2Q2 * (1 + 3 * m2Q2)
                         - 32. / 9 * (-3 - 4 * m2Q2 + 24 * m4Q4)))
         + (a11 * beta0
                * (32. / 3 * J * m2Q2 * (1 + 3 * m2Q2)
                   - 32. / 3 * (1 + 6 * m2Q2))
            + a10 * a11
                  * (-64. / 3 * J * m2Q2 * (1 + 3 * m2Q2)
                     + 64. / 3 * (1 + 6 * m2Q2)))
               * Lmu2
         + a11 * beta0
               * (-32. / 3 * K * m2Q2 * (1 + 3 * m2Q2)
                  - 128. / 27 * (17 + 120 * m2Q2)
                  - 256. / 9 * (-1 + 12 * m2Q2) * ln2
                  + 56. / 3 * (-1 + 12 * m2Q2) * zeta3
                  + II
                        * (32. / 3 * Logxi * m2Q2 * (1 + 3 * m2Q2)
                           + 16. / 9 * (-3 - 4 * m2Q2 + 24 * m4Q4)
                           + 256. / 3 * m2Q2 * (1 + 3 * m2Q2) * ln2
                           - 56 * m2Q2 * (1 + 3 * m2Q2) * zeta3)
                  + J
                        * (16. / 27 * (3 + 136 * m2Q2 + 480 * m4Q4)
                           + 64. / 9 * (-3 - 4 * m2Q2 + 24 * m4Q4) * ln2
                           - 14. / 3 * (-3 - 4 * m2Q2 + 24 * m4Q4) * zeta3))
         + Lmu
               * (a10 * a11
                      * (-128. / 3 * II * m2Q2 * (1 + 3 * m2Q2)
                         + 128. / 9 * (-1 + 12 * m2Q2)
                         - 32. / 9 * J * (-3 - 4 * m2Q2 + 24 * m4Q4))
                  + a11 * beta0
                        * (64. / 3 * II * m2Q2 * (1 + 3 * m2Q2)
                           - 64. / 9 * (-1 + 12 * m2Q2)
                           - 256. / 3 * (1 + 6 * m2Q2) * ln2
                           + 56 * (1 + 6 * m2Q2) * zeta3
                           + J
                                 * (16. / 9 * (-3 - 4 * m2Q2 + 24 * m4Q4)
                                    + 256. / 3 * m2Q2 * (1 + 3 * m2Q2) * ln2
                                    - 56 * m2Q2 * (1 + 3 * m2Q2) * zeta3))))
        / x / (1. + 4. * m2Q2);

    double delta = fabs(central_value - error);

    if (v == 1)
        return central_value + delta;
    if (v == -1)
        return central_value - delta;
    else {
        std::cout << "CL_g3_highenergyNLL: Choose either v=0, v=1 or "
                     "v=-1!!\nExiting!!\n"
                  << std::endl;
        exit(-1);
    }
}

//==========================================================================================//
//  High energy limit of the gluon coefficient function for FL at O(alpha_s^3).
//------------------------------------------------------------------------------------------//

double CL_g3_highenergy(double x, double m2Q2, double m2mu2, int nf, int v) {

    return CL_g3_highenergyLL(x, m2Q2, m2mu2)
           + CL_g3_highenergyNLL(x, m2Q2, m2mu2, nf, v);
}

//==========================================================================================//
//  High energy limit of the quark coefficient function for FL at O(alpha_s^3)
//  at leading log.
//------------------------------------------------------------------------------------------//

double CL_ps3_highenergyLL(double x, double m2Q2, double m2mu2) {

    return CF / CA * CL_g3_highenergyLL(x, m2Q2, m2mu2);
}

//==========================================================================================//
//  High energy limit of the quark coefficient function for FL at O(alpha_s^3)
//  at next to leading log.
//------------------------------------------------------------------------------------------//

double
CL_ps3_highenergyNLL(double x, double m2Q2, double m2mu2, int nf, int v) {

    return CF / CA * CL_g3_highenergyNLL(x, m2Q2, m2mu2, nf, v);
}

//==========================================================================================//
//  High energy limit of the quark coefficient function for FL at O(alpha_s^3).
//------------------------------------------------------------------------------------------//

double CL_ps3_highenergy(double x, double m2Q2, double m2mu2, int nf, int v) {

    return CF / CA * CL_g3_highenergy(x, m2Q2, m2mu2, nf, v);
}

//==========================================================================================//
//  High scale limit of the high energy limit of the gluon coefficient function
//  for F2 at O(alpha_s^3) at leading log.
//
//  Eq. (3.41) of Ref. [arXiv:1205.5727]
//------------------------------------------------------------------------------------------//

double C2_g3_highenergy_highscaleLL(double x, double m2Q2, double m2mu2) {

    double LQ = log(m2Q2);
    double LQ2 = LQ * LQ;
    double LQ3 = LQ * LQ2;

    double Lmu = log(m2mu2);
    double Lmu2 = Lmu * Lmu;

    return CA * CA * log(x)
           * (-32. / 27 * (-71. + 18 * zeta2) * LQ - 208. / 9 * LQ2
              + 32. / 9 * LQ3 + Lmu2 * (-16. / 3 + 32. / 3 * LQ)
              + Lmu
                    * (32. / 9 * (-5. + 6 * zeta2) + 416. / 9 * LQ
                       - 32. / 3 * LQ2)
              + 16. / 27 * (-92. + 78. * zeta2 - 72. * zeta3))
           / x;
}

//==========================================================================================//
//  High scale limit of the high energy limit of the gluon coefficient function
//  for F2 at O(alpha_s^3) at next to leading log.
//------------------------------------------------------------------------------------------//

double C2_g3_highenergy_highscaleNLL(
    double x, double m2Q2, double m2mu2, int nf, int v
) {

    double Lmu = log(m2mu2);
    double Lmu2 = Lmu * Lmu;

    double pi2 = M_PI * M_PI;

    double LQ = log(m2Q2);
    double LQ2 = LQ * LQ;
    double LQ3 = LQ2 * LQ;

    double a11 = a_11();
    double a21 = a_21(nf);
    double a10 = a_10(nf);

    double beta0 = beta(0, nf);

    if (v == 0) {

        return (-32. / 9 * a21 * (-5. + pi2)
                + (-416. * a21 / 9 + 64. / 27 * a10 * a11 * (-71. + 3. * pi2)
                   - 32. / 27 * a11 * beta0 * (-71. + 3. * pi2))
                      * LQ
                + (416. * a10 * a11 / 9 + 32. * a21 / 3 - 208. * a11 * beta0 / 9
                  ) * LQ2
                + (-64. * a10 * a11 / 9 + 32. * a11 * beta0 / 9) * LQ3
                + Lmu2
                      * (32. * a10 * a11 / 3 - 16. * a11 * beta0 / 3
                         + (-64 * a10 * a11 / 3 + 32. * a11 * beta0 / 3) * LQ)
                + Lmu
                      * (32. * a21 / 3 - 64. / 9 * a10 * a11 * (-5. + pi2)
                         + 32. / 9 * a11 * beta0 * (-5. + pi2)
                         + (-832. * a10 * a11 / 9 - 64 * a21 / 3
                            + 416. * a11 * beta0 / 9)
                               * LQ
                         + (64. * a10 * a11 / 3 - 32. * a11 * beta0 / 3) * LQ2)
                - 32. / 27 * a10 * a11 * (-92. + 13. * pi2 - 72. * zeta3)
                + 16. / 27 * a11 * beta0 * (-92. + 13. * pi2 - 72. * zeta3))
               / x;
    }

    double central_value = C2_g3_highenergy_highscaleNLL(x, m2Q2, m2mu2, nf, 0);

    double error =
        ((-64 * a10 * a11 / 9 + 32 * a11 * beta0 / 9) * LQ3
         + Lmu2
               * (32 * a10 * a11 / 3 - 16 * a11 * beta0 / 3
                  + (-64 * a10 * a11 / 3 + 32 * a11 * beta0 / 3) * LQ)
         + LQ2
               * (416 * a10 * a11 / 9
                  - 4. / 9 * a11 * beta0 * (52 + 96 * ln2 - 63 * zeta3))
         - 32. / 27 * a10 * a11 * (-92 + 13 * pi2 - 72 * zeta3)
         + 4. / 27 * a11 * beta0
               * (-368 + 52 * pi2 - 480 * ln2 + 96 * pi2 * ln2 + 27 * zeta3
                  - 63 * pi2 * zeta3)
         + Lmu
               * (-64. / 9 * a10 * a11 * (-5 + pi2)
                  + (64 * a10 * a11 / 3 - 32 * a11 * beta0 / 3) * LQ2
                  + LQ
                        * (-832 * a10 * a11 / 9
                           + 8. / 9 * a11 * beta0 * (52 + 96 * ln2 - 63 * zeta3)
                        )
                  + 4. / 9 * a11 * beta0
                        * (-40 + 8 * pi2 - 96 * ln2 + 63 * zeta3))
         + LQ
               * (64. / 27 * a10 * a11 * (-71 + 3 * pi2)
                  - 4. / 27 * a11 * beta0
                        * (-568 + 24 * pi2 - 1248 * ln2 + 819 * zeta3)))
        / x;

    double delta = fabs(central_value - error);

    if (v == 1)
        return central_value + delta;
    if (v == -1)
        return central_value - delta;
    else {
        std::cout << "C2_g3_highenergy_highscaleNLL: Choose either v=0, v=1 or "
                     "v=-1!!\nExiting!!\n"
                  << std::endl;
        exit(-1);
    }
}

//==========================================================================================//
//  High scale limit of the high energy limit of the gluon coefficient function
//  for F2 at O(alpha_s^3).
//------------------------------------------------------------------------------------------//

double
C2_g3_highenergy_highscale(double x, double m2Q2, double m2mu2, int nf, int v) {

    return C2_g3_highenergy_highscaleLL(x, m2Q2, m2mu2)
           + C2_g3_highenergy_highscaleNLL(x, m2Q2, m2mu2, nf, v);
}

//==========================================================================================//
//  High scale limit of the high energy limit of the quark coefficient function
//  for F2 at O(alpha_s^3) at leading log.
//------------------------------------------------------------------------------------------//

double C2_ps3_highenergy_highscaleLL(double x, double m2Q2, double m2mu2) {

    return CF / CA * C2_g3_highenergy_highscaleLL(x, m2Q2, m2mu2);
}

//==========================================================================================//
//  High scale limit of the high energy limit of the quark coefficient function
//  for F2 at O(alpha_s^3) at next to leading log.
//------------------------------------------------------------------------------------------//

double C2_ps3_highenergy_highscaleNLL(
    double x, double m2Q2, double m2mu2, int nf, int v
) {

    return CF / CA * C2_g3_highenergy_highscaleNLL(x, m2Q2, m2mu2, nf, v);
}

//==========================================================================================//
//  High scale limit of the high energy limit of the quark coefficient function
//  for F2 at O(alpha_s^3).
//------------------------------------------------------------------------------------------//

double C2_ps3_highenergy_highscale(
    double x, double m2Q2, double m2mu2, int nf, int v
) {

    return CF / CA * C2_g3_highenergy_highscale(x, m2Q2, m2mu2, nf, v);
}

//==========================================================================================//
//  High scale limit of the high energy limit of the gluon coefficient function
//  for FL at O(alpha_s^3) at leading log.
//------------------------------------------------------------------------------------------//

double CL_g3_highenergy_highscaleLL(double x, double m2Q2, double m2mu2) {

    double Lmu = log(m2mu2);
    double Lmu2 = Lmu * Lmu;

    double LQ = log(m2Q2);
    double LQ2 = LQ * LQ;

    return CA * CA
           * (32. / 27 * (-68. + 18. * zeta2) - 32. / 3 * Lmu2 - 64. / 9 * LQ
              - 32. / 3 * LQ2 + Lmu * (64. / 9 + 64. / 3 * LQ))
           * log(x) / x;
}

//==========================================================================================//
//  High scale limit of the high energy limit of the gluon coefficient function
//  for FL at O(alpha_s^3) at next to leading log.
//------------------------------------------------------------------------------------------//

double CL_g3_highenergy_highscaleNLL(
    double x, double m2Q2, double m2mu2, int nf, int v
) {

    double Lmu = log(m2mu2);
    double Lmu2 = Lmu * Lmu;

    double pi2 = M_PI * M_PI;

    double LQ = log(m2Q2);
    double LQ2 = LQ * LQ;

    double a11 = a_11();
    double a21 = a_21(nf);
    double a10 = a_10(nf);

    double beta0 = beta(0, nf);

    if (v == 0) {

        return (-64. * a21 / 9 - 64. / 27 * a10 * a11 * (-68. + 3 * pi2)
                + 32. / 27 * a11 * beta0 * (-68. + 3 * pi2)
                + (64. * a10 * a11 / 3 - 32. * a11 * beta0 / 3) * Lmu2
                + (128. * a10 * a11 / 9 - 64. * a21 / 3 - 64. * a11 * beta0 / 9)
                      * LQ
                + (64. * a10 * a11 / 3 - 32. * a11 * beta0 / 3) * LQ2
                + Lmu
                      * (-128. * a10 * a11 / 9 + 64. * a21 / 3
                         + 64. * a11 * beta0 / 9
                         + (-128. * a10 * a11 / 3 + 64. * a11 * beta0 / 3) * LQ)
               )
               / x;
    }

    double central_value = CL_g3_highenergy_highscaleNLL(x, m2Q2, m2mu2, nf, 0);

    double error =
        (-64. / 27 * a10 * a11 * (-68 + 3 * pi2)
         + (64 * a10 * a11 / 3 - 32 * a11 * beta0 / 3) * Lmu2
         + (64 * a10 * a11 / 3 - 32 * a11 * beta0 / 3) * LQ2
         + Lmu
               * (-128 * a10 * a11 / 9
                  + (-128 * a10 * a11 / 3 + 64 * a11 * beta0 / 3) * LQ
                  - 8. / 9 * a11 * beta0 * (-8 + 96 * ln2 - 63 * zeta3))
         + LQ
               * (128 * a10 * a11 / 9
                  + 8. / 9 * a11 * beta0 * (-8 + 96 * ln2 - 63 * zeta3))
         + 8. / 27 * a11 * beta0 * (-272 + 12 * pi2 + 96 * ln2 - 63 * zeta3))
        / x;

    double delta = fabs(central_value - error);

    if (v == 1)
        return central_value + delta;
    if (v == -1)
        return central_value - delta;
    else {
        std::cout << "CL_g3_highenergy_highscaleNLL: choose either v=0, v=1 or "
                     "v=-1!!\nExiting!!\n"
                  << std::endl;
        exit(-1);
    }
}

//==========================================================================================//
//  High scale limit of the high energy limit of the gluon coefficient function
//  for FL at O(alpha_s^3).
//------------------------------------------------------------------------------------------//

double
CL_g3_highenergy_highscale(double x, double m2Q2, double m2mu2, int nf, int v) {

    return CL_g3_highenergy_highscaleLL(x, m2Q2, m2mu2)
           + CL_g3_highenergy_highscaleNLL(x, m2Q2, m2mu2, nf, v);
}

//==========================================================================================//
//  High scale limit of the high energy limit of the quark coefficient function
//  for F2 at O(alpha_s^3) at leading log.
//------------------------------------------------------------------------------------------//

double CL_ps3_highenergy_highscaleLL(double x, double m2Q2, double m2mu2) {

    return CF / CA * CL_g3_highenergy_highscaleLL(x, m2Q2, m2mu2);
}

//==========================================================================================//
//  High scale limit of the high energy limit of the quark coefficient function
//  for F2 at O(alpha_s^3) at next to leading log.
//------------------------------------------------------------------------------------------//

double CL_ps3_highenergy_highscaleNLL(
    double x, double m2Q2, double m2mu2, int nf, int v
) {

    return CF / CA * CL_g3_highenergy_highscaleNLL(x, m2Q2, m2mu2, nf, v);
}

//==========================================================================================//
//  High scale limit of the high energy limit of the quark coefficient function
//  for F2 at O(alpha_s^3).
//------------------------------------------------------------------------------------------//

double CL_ps3_highenergy_highscale(
    double x, double m2Q2, double m2mu2, int nf, int v
) {

    return CF / CA * CL_g3_highenergy_highscale(x, m2Q2, m2mu2, nf, v);
}

//==========================================================================================//
//  Power terms in the small x limit of the gluon coefficient function for F2 at
//  O(alpha_s^3) at leading log.
//------------------------------------------------------------------------------------------//

double C2_g3_power_termsLL(double x, double m2Q2, double m2mu2) {

    if (x < 0 || x >= 1)
        return 0;

    return C2_g3_highenergyLL(x, m2Q2, m2mu2)
           - C2_g3_highenergy_highscaleLL(x, m2Q2, m2mu2);
}

//==========================================================================================//
//  Power terms in the small x limit of the gluon coefficient function for F2 at
//  O(alpha_s^3).
//------------------------------------------------------------------------------------------//

double C2_g3_power_terms(double x, double m2Q2, double m2mu2, int nf, int v) {

    if (x < 0 || x >= 1)
        return 0;

    return C2_g3_highenergy(x, m2Q2, m2mu2, nf, v)
           - C2_g3_highenergy_highscale(x, m2Q2, m2mu2, nf, v);
}

//==========================================================================================//
//  Power terms in the small x limit of the quark coefficient function for F2 at
//  O(alpha_s^2) at leading log.
//------------------------------------------------------------------------------------------//

double C2_ps3_power_termsLL(double x, double m2Q2, double m2mu2) {

    return CF / CA * C2_g3_power_termsLL(x, m2Q2, m2mu2);
}

//==========================================================================================//
//  Power terms in the small x limit of the quark coefficient function for F2 at
//  O(alpha_s^3).
//------------------------------------------------------------------------------------------//

double C2_ps3_power_terms(double x, double m2Q2, double m2mu2, int nf, int v) {

    return CF / CA * C2_g3_power_terms(x, m2Q2, m2mu2, nf, v);
}

//==========================================================================================//
//  Power terms in the small x limit of the gluon coefficient function for FL at
//  O(alpha_s^3).
//------------------------------------------------------------------------------------------//

double CL_g3_power_terms(double x, double m2Q2, double m2mu2, int nf, int v) {

    if (x < 0 || x >= 1)
        return 0;

    return CL_g3_highenergy(x, m2Q2, m2mu2, nf, v)
           - CL_g3_highenergy_highscale(x, m2Q2, m2mu2, nf, v);
}

//==========================================================================================//
//  Power terms in the small x limit of the gluon coefficient function for FL at
//  O(alpha_s^3).
//------------------------------------------------------------------------------------------//

double CL_ps3_power_terms(double x, double m2Q2, double m2mu2, int nf, int v) {

    return CF / CA * CL_g3_power_terms(x, m2Q2, m2mu2, nf, v);
}

//==========================================================================================//
//                  Color factors O(alpha_s^3)
//------------------------------------------------------------------------------------------//

double a_10(int nf) {

    return -(11. * CA + 2. * nf * (1. - 2. * CF / CA)) / 12.;
}

double a_11() { return CA; }

double a_21(int nf) { return nf * (26. * CF - 23. * CA) / 36.; }
