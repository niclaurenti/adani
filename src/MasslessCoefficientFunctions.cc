#include "adani/MasslessCoefficientFunctions.h"
#include "adani/Constants.h"
#include "adani/SpecialFunctions.h"
#include <cmath>

//==========================================================================================//
//  Massless gluon coefficient functions for F2 at
//  O(alpha_s) for mu=Q.
//
//  Eq. (4.4) of Ref. [hep-ph/0504242v1]
//------------------------------------------------------------------------------------------//

double C2_g1_massless(double x, int nf) {

    return 4. * nf * TR
           * (-8. * x * x + 8. * x - 1.
              + log((1. - x) / x) * (2. * x * x - 2. * x + 1.));
}

//==========================================================================================//
//  Massless gluon coefficient functions for FL at
//  O(alpha_s) for mu=Q.
//
// Eq. (3) of Ref. [arXiv:hep-ph/0411112v2]
//------------------------------------------------------------------------------------------//

double CL_g1_massless(double x, int nf) { return 16. * nf * TR * x * (1. - x); }

//==========================================================================================//
//  Massless gluon coefficient functions for F2 at
//  O(alpha_s^2) for mu=Q.
//
//  Eq. (B.6) of Ref. [hep-ph/0504242v1]
//------------------------------------------------------------------------------------------//

double C2_g2_massless(double x, int nf) {

    double x2 = x * x;
    double x3 = x2 * x;

    double Hm1 = H(x, -1);
    double H0 = H(x, 0);
    double H1 = H(x, 1);

    double H00 = H(x, 0, 0);
    double H10 = H(x, 1, 0);
    double Hm10 = H(x, -1, 0);
    double H01 = H(x, 0, 1);
    double H11 = H(x, 1, 1);

    double Hm1m10 = H(x, -1, -1, 0);
    double H0m10 = H(x, 0, -1, 0);
    double Hm100 = H(x, -1, 0, 0);
    double H000 = H(x, 0, 0, 0);
    double H100 = H(x, 1, 0, 0);
    double H010 = H(x, 0, 1, 0);
    double H110 = H(x, 1, 1, 0);
    double Hm101 = H(x, -1, 0, 1);
    double H001 = H(x, 0, 0, 1);
    double H101 = H(x, 1, 0, 1);
    double H011 = H(x, 0, 1, 1);
    double H111 = H(x, 1, 1, 1);

    return nf * CF
               * (-647. / 15 - 104. / 3 * zeta2 * x + 72 * zeta2 * x2
                  + 96. / 5 * zeta2 * x3 + 16 * zeta2 + 72 * zeta3 * x2
                  + 32 * zeta3 + 8. / 15 / x + 239. / 5 * x - 36. / 5 * x2
                  + 32 * H0m10 + 32 * H0m10 * x2 - 32 * Hm1 * zeta2 * x
                  - 16 * Hm1 * zeta2 * x2 - 16 * Hm1 * zeta2 - 32 * Hm1m10
                  - 64 * Hm1m10 * x - 32 * Hm1m10 * x2 + 48 * Hm10
                  + 8. / 15 * Hm10 / x2 + 64. / 3 * Hm10 * x
                  + 96. / 5 * Hm10 * x3 + 16 * Hm100 + 32 * Hm100 * x
                  + 16 * Hm100 * x2 - 236. / 15 * H0 - 32 * H0 * zeta2 * x
                  + 48 * H0 * zeta2 * x2 + 16 * H0 * zeta2 - 8. / 15 * H0 / x
                  + 113. / 5 * H0 * x - 216. / 5 * H0 * x2 - 3 * H00
                  + 44. / 3 * H00 * x - 72 * H00 * x2 - 96. / 5 * H00 * x3
                  - 10 * H000 + 20 * H000 * x - 40 * H000 * x2 - 14 * H1
                  - 16 * H1 * zeta2 * x + 32 * H1 * zeta2 * x2 + 8 * H1 * zeta2
                  + 40 * H1 * x - 24 * H1 * x2 - 26 * H10 + 80 * H10 * x
                  - 72 * H10 * x2 - 4 * H100 + 8 * H100 * x - 24 * H100 * x2
                  - 26 * H11 + 80 * H11 * x - 72 * H11 * x2 - 16 * H110
                  + 32 * H110 * x - 32 * H110 * x2 - 20 * H111 + 40 * H111 * x
                  - 40 * H111 * x2 - 24 * H101 + 48 * H101 * x - 48 * H101 * x2
                  - 16 * H01 + 56 * H01 * x - 72 * H01 * x2 - 12 * H010
                  + 24 * H010 * x - 32 * H010 * x2 - 16 * H011 + 32 * H011 * x
                  - 40 * H011 * x2 - 16 * H001 + 32 * H001 * x - 48 * H001 * x2)
           + nf * CA
                 * (+239. / 9 - 16. / 3 * zeta2 / x - 144 * zeta2 * x
                    + 148 * zeta2 * x2 + 8 * zeta2 - 48 * zeta3 * x
                    + 24 * zeta3 * x2 + 4 * zeta3 + 344. / 27 / x
                    + 1072. / 9 * x - 4493. / 27 * x2 + 16 * H0m10 * x2
                    - 8 * Hm1 * zeta2 * x - 16 * Hm1 * zeta2 * x2
                    - 4 * Hm1 * zeta2 + 8 * Hm1m10 + 16 * Hm1m10 * x - 24 * Hm10
                    - 16. / 3 * Hm10 / x + 80. / 3 * Hm10 * x2 + 8 * Hm100
                    + 16 * Hm100 * x + 24 * Hm100 * x2 + 8 * Hm101
                    + 16 * Hm101 * x + 16 * Hm101 * x2 + 58 * H0
                    - 64 * H0 * zeta2 * x + 16 * H0 * zeta2 * x2
                    - 8 * H0 * zeta2 + 584. / 3 * H0 * x - 2090. / 9 * H0 * x2
                    - 2 * H00 + 176 * H00 * x - 388. / 3 * H00 * x2 + 20 * H000
                    + 56 * H000 * x + 62. / 3 * H1 - 16 * H1 * zeta2 * x
                    + 8 * H1 * zeta2 * x2 + 8 * H1 * zeta2 - 104. / 9 * H1 / x
                    + 454. / 3 * H1 * x - 1570. / 9 * H1 * x2 - 4 * H10
                    + 16. / 3 * H10 / x + 80 * H10 * x - 268. / 3 * H10 * x2
                    - 12 * H100 + 24 * H100 * x - 16 * H100 * x2 - 4 * H11
                    + 16. / 3 * H11 / x + 72 * H11 * x - 244. / 3 * H11 * x2
                    - 12 * H110 + 24 * H110 * x - 24 * H110 * x2 - 4 * H111
                    + 8 * H111 * x - 8 * H111 * x2 - 4 * H101 + 8 * H101 * x
                    - 8 * H101 * x2 - 8 * H01 + 144 * H01 * x - 148 * H01 * x2
                    + 48 * H010 * x - 16 * H010 * x2 + 48 * H011 * x
                    - 16 * H011 * x2 + 8 * H001 + 64 * H001 * x
                    - 16 * H001 * x2);
}

//==========================================================================================//
//  Massless quark coefficient functions for F2 at
//  O(alpha_s^2) for mu=Q.
//
//  Eq. (B.7) of Ref. [hep-ph/0504242v1]
//------------------------------------------------------------------------------------------//

double C2_ps2_massless(double x, int nf) {

    double x2 = x * x;

    double H0 = H(x, 0);
    double H1 = H(x, 1);

    double H00 = H(x, 0, 0);
    double H10 = H(x, 1, 0);
    double Hm10 = H(x, -1, 0);
    double H01 = H(x, 0, 1);
    double H11 = H(x, 1, 1);

    double H000 = H(x, 0, 0, 0);
    double H010 = H(x, 0, 1, 0);
    double H001 = H(x, 0, 0, 1);
    double H011 = H(x, 0, 1, 1);

    return nf * CF
           * (+158. / 9 - 16. / 3 * zeta2 / x - 16 * zeta2 * x + 16 * zeta2 * x2
              - 8 * zeta3 * x - 8 * zeta3 + 344. / 27 / x - 422. / 9 * x
              + 448. / 27 * x2 - 16 * Hm10 - 16. / 3 * Hm10 / x - 16 * Hm10 * x
              - 16. / 3 * Hm10 * x2 + 56 * H0 - 16 * H0 * zeta2 * x
              - 16 * H0 * zeta2 - 88. / 3 * H0 * x - 128. / 9 * H0 * x2
              - 2 * H00 + 30 * H00 * x - 64. / 3 * H00 * x2 + 20 * H000
              + 20 * H000 * x + 104. / 3 * H1 - 104. / 9 * H1 / x
              - 80. / 3 * H1 * x + 32. / 9 * H1 * x2 + 4 * H10
              + 16. / 3 * H10 / x - 4 * H10 * x - 16. / 3 * H10 * x2 + 4 * H11
              + 16. / 3 * H11 / x - 4 * H11 * x - 16. / 3 * H11 * x2
              - 16 * H01 * x2 + 8 * H010 + 8 * H010 * x + 8 * H011
              + 8 * H011 * x + 16 * H001 + 16 * H001 * x);
}

//==========================================================================================//
//  Massless gluon coefficient functions for FL at
//  O(alpha_s^2) for mu=Q.
//
//  Eq. (B.14) of Ref. [hep-ph/0504242v1]
//------------------------------------------------------------------------------------------//

double CL_g2_massless(double x, int nf) {

    double x2 = x * x;
    double x3 = x2 * x;

    double H0 = H(x, 0);
    double H1 = H(x, 1);

    double H00 = H(x, 0, 0);
    double H10 = H(x, 1, 0);
    double Hm10 = H(x, -1, 0);
    double H01 = H(x, 0, 1);
    double H11 = H(x, 1, 1);

    return nf * CF
               * (-128. / 15 + 16. / 3 * zeta2 * x + 64. / 5 * zeta2 * x3
                  + 32. / 15 / x - 304. / 5 * x + 336. / 5 * x2
                  + 32. / 15 * Hm10 / x2 - 32. / 3 * Hm10 * x
                  + 64. / 5 * Hm10 * x3 - 104. / 15 * H0 - 32. / 15 * H0 / x
                  - 208. / 5 * H0 * x + 96. / 5 * H0 * x2 - 64. / 3 * H00 * x
                  - 64. / 5 * H00 * x3 - 8 * H1 - 24 * H1 * x + 32 * H1 * x2
                  - 16 * H01 * x)
           + nf * CA
                 * (+16. / 3 - 64 * zeta2 * x + 32 * zeta2 * x2 - 16. / 9 / x
                    + 272. / 3 * x - 848. / 9 * x2 + 32 * Hm10 * x
                    + 32 * Hm10 * x2 + 16 * H0 + 128 * H0 * x - 208 * H0 * x2
                    + 96 * H00 * x + 16 * H1 - 16. / 3 * H1 / x + 144 * H1 * x
                    - 464. / 3 * H1 * x2 + 32 * H10 * x - 32 * H10 * x2
                    + 32 * H11 * x - 32 * H11 * x2 + 96 * H01 * x
                    - 32 * H01 * x2);
}

//==========================================================================================//
//  Massless gluon coefficient functions for FL at
//  O(alpha_s^2) for mu=Q.
//
//  Eq. (B.15) of Ref. [hep-ph/0504242v1]
//------------------------------------------------------------------------------------------//

double CL_ps2_massless(double x, int nf) {

    double x2 = x * x;

    double H0 = H(x, 0);
    double H1 = H(x, 1);

    double H00 = H(x, 0, 0);
    double H01 = H(x, 0, 1);

    return nf * CF
           * (+16. / 3 - 16 * zeta2 * x - 16. / 9 / x - 64. / 3 * x
              + 160. / 9 * x2 + 16 * H0 - 16 * H0 * x - 32 * H0 * x2
              + 32 * H00 * x + 16 * H1 - 16. / 3 * H1 / x - 32. / 3 * H1 * x2
              + 16 * H01 * x);
}

//==========================================================================================//
//  Massless gluon coefficient functions for F2 at
//  O(alpha_s^3) for mu=Q.
//  The term fl_g_11 is put to zero for the reason explained
//  in page 15 of arXiv:1205.5727
//
//  Eq. (B.9) of Ref. [hep-ph/0504242v1]
//------------------------------------------------------------------------------------------//

double C2_g3_massless(double x, int nf) {

    // the commented varibles are needed for the flav term

    double x2 = x * x;
    double x3 = x2 * x;

    // double xm = 1 - x;
    // double xm2 = xm * xm;
    // double xm3 = xm2 * xm;
    // double xm4 = xm3 * xm;

    // double xp = 1 + x;
    // double xp2 = xp * xp;
    // double xp3 = xp2 * xp;
    // double xp4 = xp3 * xp;

    // Allocate pointers for the harmonic polylogs
    double wx = x;
    int nw = 5;
    int n1 = -1;
    int n2 = 1;
    int sz = n2 - n1 + 1;
    double *Hr1 = new double[sz];
    double *Hr2 = new double[sz * sz];
    double *Hr3 = new double[sz * sz * sz];
    double *Hr4 = new double[sz * sz * sz * sz];
    double *Hr5 = new double[sz * sz * sz * sz * sz];

    // Call polylogs
    apf_hplog_(&wx, &nw, Hr1, Hr2, Hr3, Hr4, Hr5, &n1, &n2);

    // weight 1
    const double Hm1 = Hr1[0];
    const double H0 = Hr1[1];
    const double H1 = Hr1[2];

    // weight 2
    const double Hm1m1 = Hr2[0];
    const double H0m1 = Hr2[1];
    const double Hm10 = Hr2[3];
    const double H00 = Hr2[4];
    const double H10 = Hr2[5];
    const double H01 = Hr2[7];
    const double H11 = Hr2[8];

    // weight 3
    const double Hm1m1m1 = Hr3[0];
    const double H0m1m1 = Hr3[1];
    const double Hm10m1 = Hr3[3];
    const double H00m1 = Hr3[4];
    const double H10m1 = Hr3[5];
    const double Hm1m10 = Hr3[9];
    const double H0m10 = Hr3[10];
    const double Hm100 = Hr3[12];
    const double H000 = Hr3[13];
    const double H100 = Hr3[14];
    const double H010 = Hr3[16];
    const double H110 = Hr3[17];
    const double Hm101 = Hr3[21];
    const double H001 = Hr3[22];
    const double H101 = Hr3[23];
    const double H011 = Hr3[25];
    const double H111 = Hr3[26];

    // weight 4
    const double Hm1m1m10 = Hr4[27];
    const double H0m1m10 = Hr4[28];
    const double Hm10m10 = Hr4[30];
    const double H00m10 = Hr4[31];
    const double H10m10 = Hr4[32];
    const double Hm1m100 = Hr4[36];
    const double H0m100 = Hr4[37];
    const double Hm1000 = Hr4[39];
    const double H0000 = Hr4[40];
    const double H1000 = Hr4[41];
    const double H0100 = Hr4[43];
    const double H1100 = Hr4[44];
    const double Hm1010 = Hr4[48];
    const double H0010 = Hr4[49];
    const double H1010 = Hr4[50];
    const double H0110 = Hr4[52];
    const double H1110 = Hr4[53];
    const double Hm1m101 = Hr4[63];
    const double H0m101 = Hr4[64];
    const double Hm1001 = Hr4[66];
    const double H0001 = Hr4[67];
    const double H1001 = Hr4[68];
    const double H0101 = Hr4[70];
    const double H1101 = Hr4[71];
    const double Hm1011 = Hr4[75];
    const double H0011 = Hr4[76];
    const double H1011 = Hr4[77];
    const double H0111 = Hr4[79];
    const double H1111 = Hr4[80];

    //  weight 5
    const double Hm1m1m1m10 = Hr5[81];
    const double H0m1m1m10 = Hr5[82];
    const double Hm10m1m10 = Hr5[84];
    const double H00m1m10 = Hr5[85];
    const double H10m1m10 = Hr5[86];
    const double Hm1m10m10 = Hr5[90];
    const double H0m10m10 = Hr5[91];
    const double Hm100m10 = Hr5[93];
    const double H000m10 = Hr5[94];
    const double H100m10 = Hr5[95];
    const double H010m10 = Hr5[97];
    const double H110m10 = Hr5[98];
    const double Hm1m1m100 = Hr5[108];
    const double H0m1m100 = Hr5[109];
    const double Hm10m100 = Hr5[111];
    const double H00m100 = Hr5[112];
    const double H10m100 = Hr5[113];
    const double Hm1m1000 = Hr5[117];
    const double H0m1000 = Hr5[118];
    const double Hm10000 = Hr5[120];
    const double H00000 = Hr5[121];
    const double H10000 = Hr5[122];
    const double H01000 = Hr5[124];
    const double H11000 = Hr5[125];
    const double Hm10100 = Hr5[129];
    const double H00100 = Hr5[130];
    const double H10100 = Hr5[131];
    const double H01100 = Hr5[133];
    const double H11100 = Hr5[134];
    const double Hm1m1010 = Hr5[144];
    const double H0m1010 = Hr5[145];
    const double Hm10010 = Hr5[147];
    const double H00010 = Hr5[148];
    const double H10010 = Hr5[149];
    const double H01010 = Hr5[151];
    const double H11010 = Hr5[152];
    const double Hm10110 = Hr5[156];
    const double H00110 = Hr5[157];
    const double H10110 = Hr5[158];
    const double H01110 = Hr5[160];
    const double H11110 = Hr5[161];
    const double Hm1m1m101 = Hr5[189];
    const double H0m1m101 = Hr5[190];
    const double Hm10m101 = Hr5[192];
    const double H00m101 = Hr5[193];
    const double H10m101 = Hr5[194];
    const double Hm1m1001 = Hr5[198];
    const double H0m1001 = Hr5[199];
    const double Hm10001 = Hr5[201];
    const double H00001 = Hr5[202];
    const double H10001 = Hr5[203];
    const double H01001 = Hr5[205];
    const double H11001 = Hr5[206];
    const double Hm10101 = Hr5[210];
    const double H00101 = Hr5[211];
    const double H10101 = Hr5[212];
    const double H01101 = Hr5[214];
    const double H11101 = Hr5[215];
    const double Hm1m1011 = Hr5[225];
    const double H0m1011 = Hr5[226];
    const double Hm10011 = Hr5[228];
    const double H00011 = Hr5[229];
    const double H10011 = Hr5[230];
    const double H01011 = Hr5[232];
    const double H11011 = Hr5[233];
    const double Hm10111 = Hr5[237];
    const double H00111 = Hr5[238];
    const double H10111 = Hr5[239];
    const double H01111 = Hr5[241];
    const double H11111 = Hr5[242];

    delete[] Hr1;
    delete[] Hr2;
    delete[] Hr3;
    delete[] Hr4;
    delete[] Hr5;

    // double flav = 5. / 48 * fl11g(nf) * nf * nf;
    // this contribution is neglected
    // If you want to use it please check the 5/48 since
    // I'm not sure about it

    return /*flav
               * (-1216. / 75 - 512. / 75 / xm4 * zeta2 * zeta2
                  + 512. / 75 / xm3 * zeta2 * zeta2 - 256. / 5 / xm3 * zeta3
                  - 128. / 15 / xm2 * zeta2 - 1792. / 75 / xm2 * zeta2 * zeta2
                  + 128. / 5 / xm2 * zeta3 - 256. / 45 / xm * zeta2
                  - 512. / 3 / xm * zeta3 + 896. / 25 / xp4 * zeta2 * zeta2
                  - 896. / 25 / xp3 * zeta2 * zeta2 + 512. / 15 / xp3 * zeta3
                  + 128. / 15 / xp2 * zeta2 + 3136. / 25 / xp2 * zeta2 * zeta2
                  - 256. / 15 / xp2 * zeta3 + 1024. / 9 / xp * zeta3
                  - 128. / 45 / xp + 7168. / 225 * zeta2 / x
                  - 50048. / 75 * zeta2 * x - 167552. / 75 * zeta2 * x2
                  - 1152. / 5 * zeta2 * x3 + 6208. / 225 * zeta2
                  + 106496. / 75 * zeta2 * zeta2 * x
                  - 113152. / 75 * zeta2 * zeta2 * x2
                  + 43008. / 25 * zeta2 * zeta2 * x3
                  - 7616. / 75 * zeta2 * zeta2 + 896. / 15 * zeta3 / x
                  - 20608. / 15 * zeta3 * x - 63392. / 45 * zeta3 * x2
                  - 6272. / 5 * zeta3 * x3 + 4288. / 45 * zeta3
                  + 1696. / 225 / x - 45184. / 75 * x + 51968. / 75 * x2
                  - 1792. / 15 * H00m10 + 512. / 15 * H00m10 / xp4
                  - 512. / 15 * H00m10 / xp3 + 1792. / 15 * H00m10 / xp2
                  + 512. / 5 * H00m10 * x + 19456. / 15 * H00m10 * x2
                  + 1024. / 15 * H0m1 / xm4 * zeta2
                  - 1024. / 15 * H0m1 / xm3 * zeta2
                  + 3584. / 15 * H0m1 / xm2 * zeta2
                  - 22784. / 15 * H0m1 * zeta2 * x
                  + 28928. / 15 * H0m1 * zeta2 * x2 + 6784. / 15 * H0m1 * zeta2
                  + 2304. / 5 * H0m1m10 + 8704. / 15 * H0m1m10 * x
                  + 8704. / 15 * H0m1m10 * x2 - 34688. / 45 * H0m10
                  + 512. / 15 * H0m10 / xp3 - 256. / 15 * H0m10 / xp2
                  + 1024. / 9 * H0m10 / xp - 7936. / 15 * H0m10 * x
                  - 141632. / 45 * H0m10 * x2 - 12544. / 25 * H0m10 * x3
                  - 1024. / 3 * H0m100 - 512. / 15 * H0m100 / xm4
                  + 512. / 15 * H0m100 / xm3 - 1792. / 15 * H0m100 / xm2
                  + 3072. / 5 * H0m100 * x - 3328. / 3 * H0m100 * x2
                  - 3328. / 15 * H0m101 - 1024. / 15 * H0m101 / xm4
                  + 1024. / 15 * H0m101 / xm3 - 3584. / 15 * H0m101 / xm2
                  + 27136. / 15 * H0m101 * x - 8192. / 5 * H0m101 * x2
                  + 1024. / 15 * Hm1 / xm3 * zeta2
                  - 512. / 15 * Hm1 / xm2 * zeta2 + 2048. / 9 * Hm1 / xm * zeta2
                  + 1568. / 75 * Hm1 * zeta2 / x2 + 9344. / 15 * Hm1 * zeta2 * x
                  - 13216. / 45 * Hm1 * zeta2 * x2
                  + 18816. / 25 * Hm1 * zeta2 * x3 + 63296. / 45 * Hm1 * zeta2
                  + 3584. / 3 * Hm1 * zeta3 * x + 3328. / 3 * Hm1 * zeta3 * x2
                  + 1024. / 3 * Hm1 * zeta3 + 1024. / 3 * Hm10m10
                  + 3584. / 3 * Hm10m10 * x + 3328. / 3 * Hm10m10 * x2
                  - 3584. / 3 * Hm1m1 * zeta2 * x
                  - 3328. / 3 * Hm1m1 * zeta2 * x2 - 1024. / 3 * Hm1m1 * zeta2
                  - 2048. / 3 * Hm1m1m10 - 7168. / 3 * Hm1m1m10 * x
                  - 6656. / 3 * Hm1m1m10 * x2 + 69248. / 45 * Hm1m10
                  + 3136. / 225 * Hm1m10 / x2 + 22272. / 5 * Hm1m10 * x
                  + 141632. / 45 * Hm1m10 * x2 + 12544. / 25 * Hm1m10 * x3
                  + 1024. / 3 * Hm1m100 + 3584. / 3 * Hm1m100 * x
                  + 3328. / 3 * Hm1m100 * x2 - 472672. / 225 * Hm10
                  + 256. / 15 * Hm10 / xp2 + 2560. / 3 * Hm10 * zeta2 * x
                  + 2816. / 3 * Hm10 * zeta2 * x2 + 512. / 3 * Hm10 * zeta2
                  - 32. / 5 * Hm10 / x2 - 3136. / 225 * Hm10 / x
                  - 569152. / 225 * Hm10 * x - 84032. / 75 * Hm10 * x2
                  - 1152. / 5 * Hm10 * x3 - 1088 * Hm100
                  - 512. / 15 * Hm100 / xm3 + 256. / 15 * Hm100 / xm2
                  - 1024. / 9 * Hm100 / xm - 3136. / 225 * Hm100 / x2
                  - 21376. / 15 * Hm100 * x - 640 * Hm100 * x2
                  - 12544. / 25 * Hm100 * x3 - 512. / 3 * Hm1000
                  - 2560. / 3 * Hm1000 * x - 2816. / 3 * Hm1000 * x2
                  - 28672. / 45 * Hm101 - 1024. / 15 * Hm101 / xm3
                  + 512. / 15 * Hm101 / xm2 - 2048. / 9 * Hm101 / xm
                  - 3136. / 225 * Hm101 / x2 + 24064. / 15 * Hm101 * x
                  + 84032. / 45 * Hm101 * x2 - 12544. / 25 * Hm101 * x3
                  + 64. / 75 * H0 - 256. / 5 * H0 / xm4 * zeta3
                  - 256. / 5 * H0 / xm3 * zeta2 + 256. / 5 * H0 / xm3 * zeta3
                  + 128. / 5 * H0 / xm2 * zeta2 - 896. / 5 * H0 / xm2 * zeta3
                  - 512. / 3 * H0 / xm * zeta2 + 512. / 15 * H0 / xp4 * zeta3
                  + 256. / 15 * H0 / xp3 * zeta2 - 512. / 15 * H0 / xp3 * zeta3
                  - 128. / 15 * H0 / xp2 * zeta2 + 1792. / 15 * H0 / xp2 * zeta3
                  - 128. / 15 * H0 / xp2 + 512. / 9 * H0 / xp * zeta2
                  - 896. / 15 * H0 * zeta2 / x - 36224. / 15 * H0 * zeta2 * x
                  + 36160. / 9 * H0 * zeta2 * x2 - 25088. / 25 * H0 * zeta2 * x3
                  + 4544. / 45 * H0 * zeta2 + 5632. / 3 * H0 * zeta3 * x
                  - 38528. / 15 * H0 * zeta3 * x2 + 10752. / 5 * H0 * zeta3 * x3
                  + 896. / 15 * H0 * zeta3 - 1696. / 225 * H0 / x
                  - 86272. / 225 * H0 * x + 28928. / 15 * H0 * x2
                  + 2144. / 75 * H00 - 256. / 5 * H00 / xm4 * zeta2
                  + 256. / 5 * H00 / xm3 * zeta2 - 896. / 5 * H00 / xm2 * zeta2
                  + 128. / 45 * H00 / xm + 256. / 15 * H00 / xp4 * zeta2
                  - 256. / 15 * H00 / xp3 * zeta2
                  + 896. / 15 * H00 / xp2 * zeta2 - 128. / 15 * H00 / xp2
                  - 5888. / 15 * H00 * zeta2 * x
                  + 43904. / 15 * H00 * zeta2 * x2
                  - 10752. / 5 * H00 * zeta2 * x3 + 1792. / 15 * H00 * zeta2
                  + 3136. / 225 * H00 / x + 154624. / 75 * H00 * x
                  + 22784. / 25 * H00 * x2 + 1152. / 5 * H00 * x3
                  + 256. / 15 * H000 / xm3 - 128. / 15 * H000 / xm2
                  + 512. / 9 * H000 / xm - 256. / 15 * H000 / xp3
                  + 128. / 15 * H000 / xp2 - 512. / 9 * H000 / xp
                  + 7936. / 15 * H000 * x + 12544. / 25 * H000 * x3
                  + 256. / 15 * H0000 / xm4 - 256. / 15 * H0000 / xm3
                  + 896. / 15 * H0000 / xm2 - 256. / 15 * H0000 / xp4
                  + 256. / 15 * H0000 / xp3 - 896. / 15 * H0000 / xp2
                  - 512. / 5 * H0000 * x + 80704. / 225 * H1
                  + 1568. / 225 * H1 * zeta2 / x2 - 11136. / 5 * H1 * zeta2 * x
                  + 70816. / 45 * H1 * zeta2 * x2 - 6272. / 25 * H1 * zeta2 * x3
                  + 34624. / 45 * H1 * zeta2 + 2048 * H1 * zeta2 * zeta2 * x
                  - 9216. / 5 * H1 * zeta2 * zeta2 * x2
                  - 2048. / 5 * H1 * zeta2 * zeta2 - 896. / 15 * H1 * zeta3 / x2
                  + 2560. / 3 * H1 * zeta3 * x - 640 * H1 * zeta3 * x2
                  + 10752. / 5 * H1 * zeta3 * x3 - 896 * H1 * zeta3
                  + 10304. / 225 * H1 / x - 600064. / 225 * H1 * x
                  + 164992. / 75 * H1 * x2 - 512 * H10m10 * x
                  + 768 * H10m10 * x2 + 1216. / 3 * H10
                  + 896. / 15 * H10 * zeta2 / x2 - 16384. / 3 * H10 * zeta2 * x
                  + 5504 * H10 * zeta2 * x2 - 10752. / 5 * H10 * zeta2 * x3
                  + 1920 * H10 * zeta2 + 2560 * H10 * zeta3 * x
                  - 2304 * H10 * zeta3 * x2 - 512 * H10 * zeta3
                  - 6656. / 3 * H10 * x + 5824. / 3 * H10 * x2
                  - 5248. / 15 * H100 - 2560 * H100 * zeta2 * x
                  + 2304 * H100 * zeta2 * x2 + 512 * H100 * zeta2
                  - 896. / 15 * H100 / x - 2304. / 5 * H100 * x
                  + 7552. / 5 * H100 * x2 - 512. / 3 * H1000
                  + 2560. / 3 * H1000 * x - 2816. / 3 * H1000 * x2
                  + 2432. / 3 * H11 - 3584. / 3 * H11 * zeta2 * x
                  + 3328. / 3 * H11 * zeta2 * x2 + 1024. / 3 * H11 * zeta2
                  + 5120 * H11 * zeta3 * x - 4608 * H11 * zeta3 * x2
                  - 1024 * H11 * zeta3 - 13312. / 3 * H11 * x
                  + 11648. / 3 * H11 * x2 - 5120 * H110 * zeta2 * x
                  + 4608 * H110 * zeta2 * x2 + 1024 * H110 * zeta2
                  + 1408 * H1100 + 896. / 15 * H1100 / x2
                  - 10240. / 3 * H1100 * x + 3456 * H1100 * x2
                  - 10752. / 5 * H1100 * x3 + 1024 * H11100 - 5120 * H11100 * x
                  + 4608 * H11100 * x2 - 1024 * H11001 + 5120 * H11001 * x
                  - 4608 * H11001 * x2 + 512 * H10100 - 2560 * H10100 * x
                  + 2304 * H10100 * x2 - 5248. / 3 * H1001
                  - 896. / 15 * H1001 / x2 + 4608 * H1001 * x
                  - 13696. / 3 * H1001 * x2 + 10752. / 5 * H1001 * x3
                  - 512 * H10001 + 2560 * H10001 * x - 2304 * H10001 * x2
                  - 6208. / 225 * H01 + 256. / 45 * H01 / xm
                  - 4352. / 15 * H01 * zeta2 * x + 4352. / 15 * H01 * zeta2 * x2
                  + 1152. / 5 * H01 * zeta2 + 2560 * H01 * zeta3 * x
                  - 2304 * H01 * zeta3 * x2 - 512 * H01 * zeta3
                  - 10304. / 225 * H01 / x - 419008. / 225 * H01 * x
                  + 167552. / 75 * H01 * x2 - 2560 * H010 * zeta2 * x
                  + 2304 * H010 * zeta2 * x2 + 512 * H010 * zeta2
                  - 3328. / 3 * H0100 * x + 3456 * H0100 * x2
                  - 10752. / 5 * H0100 * x3 + 512 * H01100 - 2560 * H01100 * x
                  + 2304 * H01100 * x2 - 512 * H01001 + 2560 * H01001 * x
                  - 2304 * H01001 * x2 - 4544. / 45 * H001
                  + 512. / 15 * H001 / xm3 - 256. / 15 * H001 / xm2
                  + 1024. / 9 * H001 / xm + 896. / 15 * H001 / x
                  + 28288. / 15 * H001 * x - 36160. / 9 * H001 * x2
                  + 12544. / 25 * H001 * x3 - 1792. / 15 * H0001
                  + 512. / 15 * H0001 / xm4 - 512. / 15 * H0001 / xm3
                  + 1792. / 15 * H0001 / xm2 + 7424. / 15 * H0001 * x
                  - 43904. / 15 * H0001 * x2 + 10752. / 5 * H0001 * x3)*/
        +nf * CF * CA
            * (-54088421. / 48600 - 1168. / 3 * zeta2 * zeta3 * x
               - 1600. / 3 * zeta2 * zeta3 * x2 - 176 * zeta2 * zeta3
               + 1744. / 225 * zeta2 / x - 6576608. / 2025 * zeta2 * x
               - 1045922. / 2025 * zeta2 * x2 + 58736. / 75 * zeta2 * x3
               + 867496. / 2025 * zeta2 + 128. / 5 * zeta2 * zeta2 / x
               + 50084. / 45 * zeta2 * zeta2 * x
               + 24128. / 45 * zeta2 * zeta2 * x2 + 480 * zeta2 * zeta2 * x3
               - 9988. / 45 * zeta2 * zeta2 - 2032. / 45 * zeta3 / x
               - 601292. / 135 * zeta3 * x + 71092. / 45 * zeta3 * x2
               + 3824. / 5 * zeta3 * x3 + 12712. / 135 * zeta3 + 32 * zeta4 / x
               - 48 * zeta4 * x - 248 * zeta4 * x2 + 192 * zeta4
               + 24604. / 3 * zeta5 * x + 712. / 3 * zeta5 * x2
               + 926. / 3 * zeta5 + 40238. / 675 / x - 87245687. / 24300 * x
               + 3237022. / 675 * x2 - 2176. / 3 * H000m10
               + 128. / 3 * H000m10 * x - 1792. / 3 * H000m10 * x2
               - 560. / 3 * H00m1 * zeta2 * x - 3232. / 3 * H00m1 * zeta2 * x2
               - 664. / 3 * H00m1 * zeta2 + 176. / 3 * H00m1m10
               - 4256. / 3 * H00m1m10 * x - 1600. / 3 * H00m1m10 * x2
               - 2528. / 3 * H00m10 + 656 * H00m10 * x + 6224. / 3 * H00m10 * x2
               + 1536. / 5 * H00m10 * x3 - 944. / 3 * H00m100
               + 1568. / 3 * H00m100 * x + 928. / 3 * H00m100 * x2
               + 752. / 3 * H00m101 - 1568. / 3 * H00m101 * x
               + 2432. / 3 * H00m101 * x2 - 32. / 5 * H0m1 * zeta2 / x2
               - 12448. / 3 * H0m1 * zeta2 * x - 7648. / 3 * H0m1 * zeta2 * x2
               - 4032. / 5 * H0m1 * zeta2 * x3 - 776. / 3 * H0m1 * zeta2
               - 2088 * H0m1 * zeta3 * x - 2272 * H0m1 * zeta3 * x2
               - 4604. / 3 * H0m1 * zeta3 - 5248. / 3 * H0m10m10 * x
               - 1216. / 3 * H0m10m10 * x2 + 2368 * H0m1m1 * zeta2 * x
               + 2592 * H0m1m1 * zeta2 * x2 + 1792 * H0m1m1 * zeta2
               + 448. / 3 * H0m1m1m10 + 6400. / 3 * H0m1m1m10 * x
               + 640 * H0m1m1m10 * x2 + 208. / 3 * H0m1m10
               - 64. / 15 * H0m1m10 / x2 - 6656. / 3 * H0m1m10 * x
               - 6688. / 3 * H0m1m10 * x2 - 384. / 5 * H0m1m10 * x3
               - 896 * H0m1m100 - 9152. / 3 * H0m1m100 * x
               - 5168. / 3 * H0m1m100 * x2 - 5152. / 3 * H0m1m101
               - 3904. / 3 * H0m1m101 * x - 2272 * H0m1m101 * x2
               - 103424. / 45 * H0m10 - 1632 * H0m10 * zeta2 * x
               - 2224 * H0m10 * zeta2 * x2 - 4768. / 3 * H0m10 * zeta2
               + 448. / 45 * H0m10 / x2 + 128. / 15 * H0m10 / x
               - 11512. / 45 * H0m10 * x + 50828. / 45 * H0m10 * x2
               + 12928. / 25 * H0m10 * x3 - 968. / 3 * H0m100
               + 64. / 15 * H0m100 / x2 + 8528. / 3 * H0m100 * x
               + 3048 * H0m100 * x2 + 576 * H0m100 * x3 + 1472. / 3 * H0m1000
               + 4736. / 3 * H0m1000 * x + 880 * H0m1000 * x2
               + 880. / 3 * H0m101 + 64. / 15 * H0m101 / x2 + 3040 * H0m101 * x
               + 4304. / 3 * H0m101 * x2 + 768 * H0m101 * x3 + 256 * H0m1010
               + 256. / 3 * H0m1010 * x + 1120. / 3 * H0m1010 * x2
               + 288 * H0m1011 + 320. / 3 * H0m1011 * x + 448 * H0m1011 * x2
               + 4016. / 3 * H0m1001 + 736 * H0m1001 * x
               + 5264. / 3 * H0m1001 * x2 - 1276. / 75 * Hm1 * zeta2 / x2
               + 3232. / 45 * Hm1 * zeta2 / x + 144164. / 45 * Hm1 * zeta2 * x
               + 36782. / 45 * Hm1 * zeta2 * x2 - 19392. / 25 * Hm1 * zeta2 * x3
               + 83854. / 45 * Hm1 * zeta2 - 2992. / 3 * Hm1 * zeta2 * zeta2 * x
               - 10592. / 15 * Hm1 * zeta2 * zeta2 * x2
               - 2872. / 15 * Hm1 * zeta2 * zeta2 - 56. / 3 * Hm1 * zeta3 / x2
               - 320. / 3 * Hm1 * zeta3 / x - 10840. / 3 * Hm1 * zeta3 * x
               - 6668. / 3 * Hm1 * zeta3 * x2 - 672 * Hm1 * zeta3 * x3
               - 7250. / 3 * Hm1 * zeta3 - 896. / 3 * Hm100m10
               - 1792. / 3 * Hm100m10 * x - 640. / 3 * Hm100m10 * x2
               + 2800 * Hm10m1 * zeta2 * x + 1712 * Hm10m1 * zeta2 * x2
               + 1400 * Hm10m1 * zeta2 + 432 * Hm10m1m10 + 864 * Hm10m1m10 * x
               + 480 * Hm10m1m10 * x2 + 512. / 3 * Hm10m10
               - 32. / 15 * Hm10m10 / x2 - 256. / 3 * Hm10m10 / x
               - 3136. / 3 * Hm10m10 * x - 1520 * Hm10m10 * x2
               - 384. / 5 * Hm10m10 * x3 - 3176. / 3 * Hm10m100
               - 6352. / 3 * Hm10m100 * x - 3664. / 3 * Hm10m100 * x2
               - 1184 * Hm10m101 - 2368 * Hm10m101 * x - 1472 * Hm10m101 * x2
               + 112. / 5 * Hm1m1 * zeta2 / x2 + 128 * Hm1m1 * zeta2 / x
               + 12856. / 3 * Hm1m1 * zeta2 * x + 7616. / 3 * Hm1m1 * zeta2 * x2
               + 4032. / 5 * Hm1m1 * zeta2 * x3 + 8756. / 3 * Hm1m1 * zeta2
               + 11456. / 3 * Hm1m1 * zeta3 * x + 7136. / 3 * Hm1m1 * zeta3 * x2
               + 5728. / 3 * Hm1m1 * zeta3 + 1360. / 3 * Hm1m10m10
               + 2720. / 3 * Hm1m10m10 * x + 1568. / 3 * Hm1m10m10 * x2
               - 13424. / 3 * Hm1m1m1 * zeta2 * x
               - 8240. / 3 * Hm1m1m1 * zeta2 * x2 - 6712. / 3 * Hm1m1m1 * zeta2
               - 1424. / 3 * Hm1m1m1m10 - 2848. / 3 * Hm1m1m1m10 * x
               - 1696. / 3 * Hm1m1m1m10 * x2 - 808. / 3 * Hm1m1m10
               + 32. / 15 * Hm1m1m10 / x2 + 256. / 3 * Hm1m1m10 / x
               + 3248. / 3 * Hm1m1m10 * x + 1696 * Hm1m1m10 * x2
               + 384. / 5 * Hm1m1m10 * x3 + 5344. / 3 * Hm1m1m100
               + 10688. / 3 * Hm1m1m100 * x + 6272. / 3 * Hm1m1m100 * x2
               + 2000 * Hm1m1m101 + 4000 * Hm1m1m101 * x + 2464 * Hm1m1m101 * x2
               + 68468. / 45 * Hm1m10 + 11584. / 3 * Hm1m10 * zeta2 * x
               + 7168. / 3 * Hm1m10 * zeta2 * x2 + 5792. / 3 * Hm1m10 * zeta2
               - 2552. / 225 * Hm1m10 / x2 + 5504. / 45 * Hm1m10 / x
               + 6616. / 15 * Hm1m10 * x - 63356. / 45 * Hm1m10 * x2
               - 12928. / 25 * Hm1m10 * x3 - 1464 * Hm1m100 - 16 * Hm1m100 / x2
               - 512. / 3 * Hm1m100 / x - 10136. / 3 * Hm1m100 * x
               - 8728. / 3 * Hm1m100 * x2 - 576 * Hm1m100 * x3 - 968 * Hm1m1000
               - 1936 * Hm1m1000 * x - 1104 * Hm1m1000 * x2
               - 9160. / 3 * Hm1m101 - 64. / 3 * Hm1m101 / x2
               - 256. / 3 * Hm1m101 / x - 3744 * Hm1m101 * x
               - 5072. / 3 * Hm1m101 * x2 - 768 * Hm1m101 * x3
               - 752. / 3 * Hm1m1010 - 1504. / 3 * Hm1m1010 * x
               - 1120. / 3 * Hm1m1010 * x2 - 288 * Hm1m1011 - 576 * Hm1m1011 * x
               - 448 * Hm1m1011 * x2 - 4408. / 3 * Hm1m1001
               - 8816. / 3 * Hm1m1001 * x - 5552. / 3 * Hm1m1001 * x2
               - 1010516. / 225 * Hm10 - 304. / 15 * Hm10 * zeta2 / x2
               - 256. / 3 * Hm10 * zeta2 / x - 3320 * Hm10 * zeta2 * x
               - 5456. / 3 * Hm10 * zeta2 * x2 - 3648. / 5 * Hm10 * zeta2 * x3
               - 2672 * Hm10 * zeta2 - 12176. / 3 * Hm10 * zeta3 * x
               - 8528. / 3 * Hm10 * zeta3 * x2 - 4936. / 3 * Hm10 * zeta3
               + 13924. / 675 * Hm10 / x2 - 2008. / 225 * Hm10 / x
               - 4132988. / 675 * Hm10 * x - 24232. / 25 * Hm10 * x2
               + 58736. / 75 * Hm10 * x3 - 23540. / 9 * Hm100
               - 1328. / 3 * Hm100 * zeta2 * x - 176. / 3 * Hm100 * zeta2 * x2
               - 1816. / 3 * Hm100 * zeta2 + 4792. / 225 * Hm100 / x2
               - 560. / 9 * Hm100 / x - 3628 * Hm100 * x
               - 3448. / 9 * Hm100 * x2 + 21888. / 25 * Hm100 * x3
               + 796 * Hm1000 + 48. / 5 * Hm1000 / x2 + 256. / 3 * Hm1000 / x
               + 4448. / 3 * Hm1000 * x + 1384 * Hm1000 * x2
               + 1728. / 5 * Hm1000 * x3 + 552 * Hm10000 + 1104 * Hm10000 * x
               + 720 * Hm10000 * x2 - 3308. / 3 * Hm101
               - 64 * Hm101 * zeta2 * x2 + 2552. / 225 * Hm101 / x2
               - 32. / 3 * Hm101 / x - 26848. / 9 * Hm101 * x
               - 4564. / 3 * Hm101 * x2 + 12928. / 25 * Hm101 * x3
               + 1208. / 3 * Hm1010 + 32. / 15 * Hm1010 / x2
               + 1408. / 3 * Hm1010 * x + 736. / 3 * Hm1010 * x2
               + 384. / 5 * Hm1010 * x3 + 1712. / 3 * Hm10100
               + 5728. / 3 * Hm10100 * x + 4576. / 3 * Hm10100 * x2
               + 1352. / 3 * Hm1011 + 32. / 15 * Hm1011 / x2
               + 1648. / 3 * Hm1011 * x + 832. / 3 * Hm1011 * x2
               + 384. / 5 * Hm1011 * x3 + 80 * Hm10110 + 160 * Hm10110 * x
               + 160 * Hm10110 * x2 + 256. / 3 * Hm10111
               + 512. / 3 * Hm10111 * x + 512. / 3 * Hm10111 * x2
               + 272. / 3 * Hm10101 + 544. / 3 * Hm10101 * x
               + 544. / 3 * Hm10101 * x2 + 6752. / 3 * Hm1001 + 16 * Hm1001 / x2
               + 128. / 3 * Hm1001 / x + 7480. / 3 * Hm1001 * x
               + 3208. / 3 * Hm1001 * x2 + 576 * Hm1001 * x3 + 208 * Hm10010
               + 416 * Hm10010 * x + 352 * Hm10010 * x2 + 240 * Hm10011
               + 480 * Hm10011 * x + 416 * Hm10011 * x2 + 928. / 3 * Hm10001
               - 448. / 3 * Hm10001 * x - 832. / 3 * Hm10001 * x2
               - 462887. / 1215 * H0 + 112. / 15 * H0 * zeta2 / x
               - 21616. / 15 * H0 * zeta2 * x + 26818. / 45 * H0 * zeta2 * x2
               + 30256. / 25 * H0 * zeta2 * x3 - 1334. / 45 * H0 * zeta2
               + 928. / 15 * H0 * zeta2 * zeta2 * x
               - 3424. / 15 * H0 * zeta2 * zeta2 * x2
               - 1664. / 15 * H0 * zeta2 * zeta2 + 62888. / 9 * H0 * zeta3 * x
               + 21392. / 9 * H0 * zeta3 * x2 + 7872. / 5 * H0 * zeta3 * x3
               + 1316. / 9 * H0 * zeta3 + 240 * H0 * zeta4 * x
               + 48 * H0 * zeta4 * x2 + 96 * H0 * zeta4 - 9148. / 675 * H0 / x
               - 15858079. / 6075 * H0 * x + 5542108. / 2025 * H0 * x2
               - 334411. / 2025 * H00 + 256 * H00 * zeta2 * x
               + 11296. / 9 * H00 * zeta2 * x2 + 1344. / 5 * H00 * zeta2 * x3
               + 284. / 3 * H00 * zeta2 + 2032. / 3 * H00 * zeta3 * x
               + 224. / 3 * H00 * zeta3 * x2 - 320. / 3 * H00 * zeta3
               - 2872. / 225 * H00 / x + 6805949. / 2025 * H00 * x
               + 11442. / 25 * H00 * x2 - 58736. / 75 * H00 * x3
               + 18097. / 135 * H000 - 1952. / 3 * H000 * zeta2 * x
               + 160. / 3 * H000 * zeta2 * x2 - 656. / 3 * H000 * zeta2
               - 64. / 5 * H000 / x + 188434. / 135 * H000 * x
               - 75448. / 45 * H000 * x2 - 21888. / 25 * H000 * x3
               - 1636. / 9 * H0000 - 1864. / 9 * H0000 * x - 1712 * H0000 * x2
               - 2304. / 5 * H0000 * x3 + 120 * H00000 + 1168. / 3 * H00000 * x
               - 5471951. / 6075 * H1 - 1276. / 225 * H1 * zeta2 / x2
               - 2176. / 135 * H1 * zeta2 / x - 36244. / 45 * H1 * zeta2 * x
               + 3016. / 135 * H1 * zeta2 * x2 + 10864. / 25 * H1 * zeta2 * x3
               + 36224. / 45 * H1 * zeta2 + 4528. / 15 * H1 * zeta2 * zeta2 * x
               - 1696. / 15 * H1 * zeta2 * zeta2 * x2
               - 2264. / 15 * H1 * zeta2 * zeta2 - 376. / 15 * H1 * zeta3 / x2
               + 448. / 9 * H1 * zeta3 / x + 14416. / 9 * H1 * zeta3 * x
               - 6908. / 9 * H1 * zeta3 * x2 + 4512. / 5 * H1 * zeta3 * x3
               - 19994. / 9 * H1 * zeta3 + 48 * H1 * zeta4 * x
               - 48 * H1 * zeta4 * x2 - 24 * H1 * zeta4 + 66974. / 6075 * H1 / x
               - 18941999. / 6075 * H1 * x + 8390392. / 2025 * H1 * x2
               - 880. / 3 * H100m10 + 1760. / 3 * H100m10 * x
               - 1376. / 3 * H100m10 * x2 - 672 * H10m1 * zeta2 * x
               + 352 * H10m1 * zeta2 * x2 + 336 * H10m1 * zeta2
               - 544. / 3 * H10m1m10 + 1088. / 3 * H10m1m10 * x
               - 704. / 3 * H10m1m10 * x2 - 3040. / 3 * H10m10
               - 32. / 5 * H10m10 / x2 + 608 * H10m10 * x
               + 64. / 3 * H10m10 * x2 + 1152. / 5 * H10m10 * x3 - 272 * H10m100
               + 544 * H10m100 * x - 416 * H10m100 * x2 - 1280. / 3 * H10m101
               + 2560. / 3 * H10m101 * x - 1408. / 3 * H10m101 * x2
               - 145171. / 405 * H10 + 64. / 5 * H10 * zeta2 / x2
               + 128. / 9 * H10 * zeta2 / x - 464. / 3 * H10 * zeta2 * x
               - 10664. / 9 * H10 * zeta2 * x2 - 2304. / 5 * H10 * zeta2 * x3
               + 1948 * H10 * zeta2 + 2656 * H10 * zeta3 * x
               - 1568 * H10 * zeta3 * x2 - 1328 * H10 * zeta3
               + 6104. / 405 * H10 / x - 143044. / 405 * H10 * x
               + 181286. / 405 * H10 * x2 - 74332. / 135 * H100
               + 3584. / 3 * H100 * zeta2 * x - 2048. / 3 * H100 * zeta2 * x2
               - 1792. / 3 * H100 * zeta2 + 2576. / 135 * H100 / x
               + 122552. / 135 * H100 * x - 102736. / 135 * H100 * x2
               - 3188. / 9 * H1000 + 16. / 5 * H1000 / x2 - 416. / 9 * H1000 / x
               + 1576. / 9 * H1000 * x + 1480. / 9 * H1000 * x2
               - 576. / 5 * H1000 * x3 + 312 * H10000 - 624 * H10000 * x
               + 240 * H10000 * x2 - 69101. / 405 * H11
               - 16. / 15 * H11 * zeta2 / x2 - 320. / 9 * H11 * zeta2 / x
               - 1184. / 3 * H11 * zeta2 * x + 1928. / 9 * H11 * zeta2 * x2
               + 192. / 5 * H11 * zeta2 * x3 + 356 * H11 * zeta2
               + 9088. / 3 * H11 * zeta3 * x - 5344. / 3 * H11 * zeta3 * x2
               - 4544. / 3 * H11 * zeta3 + 2324. / 405 * H11 / x
               - 562724. / 405 * H11 * x + 63254. / 45 * H11 * x2
               - 1408. / 3 * H110m10 + 2816. / 3 * H110m10 * x
               - 1664. / 3 * H110m10 * x2 - 34 * H110
               + 464. / 3 * H110 * zeta2 * x - 272. / 3 * H110 * zeta2 * x2
               - 232. / 3 * H110 * zeta2 - 352. / 27 * H110 / x
               + 3820. / 9 * H110 * x - 21326. / 27 * H110 * x2
               + 176 * H110 * x3 + 4432. / 3 * H1100 + 256. / 15 * H1100 / x2
               + 640. / 9 * H1100 / x + 808 * H1100 * x
               - 17776. / 9 * H1100 * x2 - 3072. / 5 * H1100 * x3
               + 856. / 3 * H11000 - 1712. / 3 * H11000 * x
               + 368. / 3 * H11000 * x2 + 2032. / 27 * H111
               + 800. / 3 * H111 * zeta2 * x - 224. / 3 * H111 * zeta2 * x2
               - 400. / 3 * H111 * zeta2 - 784. / 27 * H111 / x
               + 2536. / 27 * H111 * x - 3556. / 9 * H111 * x2 - 220 * H1110
               + 704. / 9 * H1110 / x + 2896. / 3 * H1110 * x
               - 9632. / 9 * H1110 * x2 + 832. / 3 * H11100
               - 1664. / 3 * H11100 * x + 320. / 3 * H11100 * x2
               - 1048. / 9 * H1111 + 736. / 9 * H1111 / x
               + 5480. / 9 * H1111 * x - 2272. / 3 * H1111 * x2 - 136 * H11110
               + 272 * H11110 * x - 272 * H11110 * x2 - 104 * H11101
               + 208 * H11101 * x - 208 * H11101 * x2 - 664. / 3 * H1101
               + 704. / 9 * H1101 / x + 936 * H1101 * x - 9560. / 9 * H1101 * x2
               - 184 * H11010 + 368 * H11010 * x - 368 * H11010 * x2
               - 32 * H11011 + 64 * H11011 * x - 64 * H11011 * x2 - 384 * H11001
               + 768 * H11001 * x - 448 * H11001 * x2 - 398. / 9 * H101
               + 608. / 3 * H101 * zeta2 * x - 32. / 3 * H101 * zeta2 * x2
               - 304. / 3 * H101 * zeta2 - 1216. / 27 * H101 / x
               + 5264. / 9 * H101 * x - 19610. / 27 * H101 * x2
               - 176 * H101 * x3 - 908. / 3 * H1010 + 704. / 9 * H1010 / x
               + 3632. / 3 * H1010 * x - 11624. / 9 * H1010 * x2
               + 128. / 3 * H10100 - 256. / 3 * H10100 * x
               - 512. / 3 * H10100 * x2 - 628. / 3 * H1011
               + 256. / 3 * H1011 / x + 880 * H1011 * x - 3104. / 3 * H1011 * x2
               - 128 * H10110 + 256 * H10110 * x - 256 * H10110 * x2
               - 344. / 3 * H10101 + 688. / 3 * H10101 * x
               - 688. / 3 * H10101 * x2 - 7108. / 3 * H1001
               - 256. / 15 * H1001 / x2 + 256. / 9 * H1001 / x
               + 2944. / 3 * H1001 * x + 3920. / 9 * H1001 * x2
               + 3072. / 5 * H1001 * x3 - 80. / 3 * H10010
               + 160. / 3 * H10010 * x - 352. / 3 * H10010 * x2
               + 296. / 3 * H10011 - 592. / 3 * H10011 * x
               + 400. / 3 * H10011 * x2 + 904. / 3 * H10001
               - 1808. / 3 * H10001 * x + 1040. / 3 * H10001 * x2
               - 867496. / 2025 * H01 - 32. / 15 * H01 * zeta2 / x2
               - 104. / 3 * H01 * zeta2 * x + 768 * H01 * zeta2 * x2
               + 192. / 5 * H01 * zeta2 * x3 + 556. / 3 * H01 * zeta2
               + 696 * H01 * zeta3 * x - 2672 * H01 * zeta3 * x2
               - 3988. / 3 * H01 * zeta3 - 3752. / 225 * H01 / x
               - 5822356. / 2025 * H01 * x + 1045922. / 2025 * H01 * x2
               - 1504. / 3 * H010m10 + 128. / 3 * H010m10 * x
               - 1600. / 3 * H010m10 * x2 - 1102. / 45 * H010
               + 8320. / 3 * H010 * zeta2 * x + 3632. / 3 * H010 * zeta2 * x2
               + 1888. / 3 * H010 * zeta2 - 32. / 15 * H010 / x
               - 11792. / 45 * H010 * x - 215798. / 135 * H010 * x2
               + 176 * H010 * x3 - 352 * H0100 - 1560 * H0100 * x
               - 7168. / 3 * H0100 * x2 - 3072. / 5 * H0100 * x3
               - 344. / 3 * H01000 - 1040 * H01000 * x - 272. / 3 * H01000 * x2
               + 3494. / 135 * H011 + 976. / 3 * H011 * zeta2 * x
               - 32. / 3 * H011 * zeta2 * x2 - 40 * H011 * zeta2
               - 32. / 15 * H011 / x - 66376. / 135 * H011 * x
               - 160708. / 135 * H011 * x2 - 508. / 3 * H0110
               + 3080. / 3 * H0110 * x - 4952. / 3 * H0110 * x2 + 808 * H01100
               + 6416. / 3 * H01100 * x + 3968. / 3 * H01100 * x2
               - 1312. / 9 * H0111 + 8744. / 9 * H0111 * x
               - 14800. / 9 * H0111 * x2 - 80. / 3 * H01110
               + 2368. / 3 * H01110 * x - 976. / 3 * H01110 * x2 + 40 * H01111
               + 656 * H01111 * x - 160 * H01111 * x2 - 104. / 3 * H01101
               + 2224. / 3 * H01101 * x - 928. / 3 * H01101 * x2
               - 452. / 3 * H0101 + 1144 * H0101 * x - 5648. / 3 * H0101 * x2
               - 80. / 3 * H01010 + 2368. / 3 * H01010 * x
               - 976. / 3 * H01010 * x2 + 112. / 3 * H01011
               + 2080. / 3 * H01011 * x - 176 * H01011 * x2 - 880 * H01001
               - 5632. / 3 * H01001 * x - 1680 * H01001 * x2 + 1334. / 45 * H001
               - 160 * H001 * zeta2 * x - 64 * H001 * zeta2
               + 16. / 15 * H001 / x + 53336. / 45 * H001 * x
               - 26818. / 45 * H001 * x2 - 17328. / 25 * H001 * x3
               - 244. / 3 * H0010 + 816 * H0010 * x - 19216. / 9 * H0010 * x2
               - 384. / 5 * H0010 * x3 - 304. / 3 * H00100
               + 896. / 3 * H00100 * x - 1184. / 3 * H00100 * x2
               - 356. / 3 * H0011 + 2224. / 3 * H0011 * x
               - 21008. / 9 * H0011 * x2 - 384. / 5 * H0011 * x3 + 24 * H00110
               + 2384. / 3 * H00110 * x - 1120. / 3 * H00110 * x2
               + 248. / 3 * H00111 + 816 * H00111 * x - 736. / 3 * H00111 * x2
               + 280. / 3 * H00101 + 2608. / 3 * H00101 * x
               - 800. / 3 * H00101 * x2 - 284. / 3 * H0001 + 400 * H0001 * x
               - 11296. / 9 * H0001 * x2 + 192. / 5 * H0001 * x3
               + 472. / 3 * H00010 + 2288. / 3 * H00010 * x - 160 * H00010 * x2
               + 168 * H00011 + 2608. / 3 * H00011 * x - 544. / 3 * H00011 * x2
               + 656. / 3 * H00001 + 2080. / 3 * H00001 * x
               - 160. / 3 * H00001 * x2)
        + nf * CF * CF
              * (+149023. / 600 + 256. / 3 * zeta2 * zeta3 * x
                 - 832 * zeta2 * zeta3 * x2 - 160 * zeta2 * zeta3
                 - 592. / 45 * zeta2 / x + 166993. / 45 * zeta2 * x
                 + 14444. / 15 * zeta2 * x2 - 8528. / 25 * zeta2 * x3
                 - 4598. / 45 * zeta2 + 1844. / 15 * zeta2 * zeta2 * x
                 - 20176. / 15 * zeta2 * zeta2 * x2
                 - 12672. / 25 * zeta2 * zeta2 * x3 + 1736. / 15 * zeta2 * zeta2
                 - 368. / 15 * zeta3 / x + 267682. / 45 * zeta3 * x
                 + 3956. / 15 * zeta3 * x2 - 352 * zeta3 * x3
                 + 11107. / 15 * zeta3 + 144 * zeta4 * x - 108 * zeta4
                 - 23032. / 3 * zeta5 * x + 368 * zeta5 * x2 - 124 * zeta5
                 - 964. / 75 / x + 399401. / 300 * x - 115586. / 75 * x2
                 + 1600. / 3 * H000m10 + 4096. / 3 * H000m10 * x2
                 - 32 * H00m1 * zeta2 * x - 2624. / 3 * H00m1 * zeta2 * x2
                 - 1264. / 3 * H00m1 * zeta2 - 1184. / 3 * H00m1m10
                 + 576 * H00m1m10 * x - 2176. / 3 * H00m1m10 * x2
                 + 3488. / 3 * H00m10 - 384 * H00m10 * x
                 + 560. / 3 * H00m10 * x2 - 768. / 5 * H00m10 * x3
                 + 1888. / 3 * H00m100 - 128 * H00m100 * x
                 + 4672. / 3 * H00m100 * x2 + 224 * H00m101 + 320 * H00m101 * x
                 + 512 * H00m101 * x2 + 1728 * H0m1 * zeta2 * x
                 + 136 * H0m1 * zeta2 * x2 + 1536. / 5 * H0m1 * zeta2 * x3
                 - 1160 * H0m1 * zeta2 + 320 * H0m1 * zeta3 * x
                 - 352. / 3 * H0m1 * zeta3 * x2 + 96 * H0m1 * zeta3
                 - 256 * H0m10m10 + 1280. / 3 * H0m10m10 * x
                 - 1280. / 3 * H0m10m10 * x2 - 416 * H0m1m1 * zeta2 * x
                 - 128. / 3 * H0m1m1 * zeta2 * x2 - 208 * H0m1m1 * zeta2
                 + 288 * H0m1m1m10 - 448 * H0m1m1m10 * x
                 + 1280. / 3 * H0m1m1m10 * x2 - 1072 * H0m1m10
                 + 1024 * H0m1m10 * x + 80 * H0m1m10 * x2 - 224 * H0m1m100
                 + 1984. / 3 * H0m1m100 * x - 640 * H0m1m100 * x2
                 + 352 * H0m1m101 + 192 * H0m1m101 * x + 256 * H0m1m101 * x2
                 + 39088. / 15 * H0m10 + 640. / 3 * H0m10 * zeta2 * x
                 - 192 * H0m10 * zeta2 * x2 + 160 * H0m10 * zeta2
                 - 64. / 15 * H0m10 / x2 - 64. / 15 * H0m10 / x
                 + 29776. / 45 * H0m10 * x - 2872. / 5 * H0m10 * x2
                 - 704. / 5 * H0m10 * x3 + 4592. / 3 * H0m100
                 - 3040. / 3 * H0m100 * x - 280. / 3 * H0m100 * x2
                 - 768. / 5 * H0m100 * x3 + 592. / 3 * H0m1000
                 - 608. / 3 * H0m1000 * x + 768 * H0m1000 * x2 + 624 * H0m101
                 - 1216 * H0m101 * x - 96 * H0m101 * x2
                 - 1536. / 5 * H0m101 * x3 + 32 * H0m1010 + 64 * H0m1010 * x
                 + 128 * H0m1010 * x2 + 32 * H0m1011 + 64 * H0m1011 * x
                 + 128 * H0m1011 * x2 - 96 * H0m1001 + 64 * H0m1001 * x
                 + 192 * H0m1001 * x2 + 88. / 15 * Hm1 * zeta2 / x2
                 + 128. / 15 * Hm1 * zeta2 / x - 73544. / 15 * Hm1 * zeta2 * x
                 - 8364. / 5 * Hm1 * zeta2 * x2 + 1056. / 5 * Hm1 * zeta2 * x3
                 - 49604. / 15 * Hm1 * zeta2
                 + 3712. / 3 * Hm1 * zeta2 * zeta2 * x
                 + 12608. / 15 * Hm1 * zeta2 * zeta2 * x2
                 + 4672. / 15 * Hm1 * zeta2 * zeta2 + 32. / 5 * Hm1 * zeta3 / x2
                 - 352 * Hm1 * zeta3 * x + 104 * Hm1 * zeta3 * x2
                 + 1152. / 5 * Hm1 * zeta3 * x3 - 184 * Hm1 * zeta3
                 - 640. / 3 * Hm100m10 - 1280. / 3 * Hm100m10 * x
                 - 1664. / 3 * Hm100m10 * x2 - 64. / 3 * Hm10m1 * zeta2 * x
                 + 704. / 3 * Hm10m1 * zeta2 * x2 - 32. / 3 * Hm10m1 * zeta2
                 + 704. / 3 * Hm10m1m10 + 1408. / 3 * Hm10m1m10 * x
                 + 1408. / 3 * Hm10m1m10 * x2 - 848 * Hm10m10
                 - 2656. / 3 * Hm10m10 * x - 112. / 3 * Hm10m10 * x2
                 - 288 * Hm10m100 - 576 * Hm10m100 * x - 704 * Hm10m100 * x2
                 + 128 * Hm10m101 + 256 * Hm10m101 * x
                 - 128. / 15 * Hm1m1 * zeta2 / x2 + 320. / 3 * Hm1m1 * zeta2 * x
                 - 264 * Hm1m1 * zeta2 * x2 - 1536. / 5 * Hm1m1 * zeta2 * x3
                 + 8 * Hm1m1 * zeta2 - 800. / 3 * Hm1m1 * zeta3 * x
                 + 352. / 3 * Hm1m1 * zeta3 * x2 - 400. / 3 * Hm1m1 * zeta3
                 + 640. / 3 * Hm1m10m10 + 1280. / 3 * Hm1m10m10 * x
                 + 1280. / 3 * Hm1m10m10 * x2 + 1664. / 3 * Hm1m1m1 * zeta2 * x
                 + 128. / 3 * Hm1m1m1 * zeta2 * x2 + 832. / 3 * Hm1m1m1 * zeta2
                 - 640. / 3 * Hm1m1m1m10 - 1280. / 3 * Hm1m1m1m10 * x
                 - 1280. / 3 * Hm1m1m1m10 * x2 + 1040 * Hm1m1m10
                 + 1216 * Hm1m1m10 * x + 176 * Hm1m1m10 * x2 + 192 * Hm1m1m100
                 + 384 * Hm1m1m100 * x + 640 * Hm1m1m100 * x2 - 384 * Hm1m1m101
                 - 768 * Hm1m1m101 * x - 256 * Hm1m1m101 * x2
                 - 7288. / 3 * Hm1m10 - 192 * Hm1m10 * zeta2 * x
                 + 192 * Hm1m10 * zeta2 * x2 - 96 * Hm1m10 * zeta2
                 + 176. / 45 * Hm1m10 / x2 - 15920. / 9 * Hm1m10 * x
                 + 728 * Hm1m10 * x2 + 704. / 5 * Hm1m10 * x3 - 984 * Hm1m100
                 + 64. / 15 * Hm1m100 / x2 - 1072 * Hm1m100 * x
                 + 280. / 3 * Hm1m100 * x2 + 768. / 5 * Hm1m100 * x3
                 - 320 * Hm1m1000 - 640 * Hm1m1000 * x - 768 * Hm1m1000 * x2
                 + 512 * Hm1m101 + 128. / 15 * Hm1m101 / x2
                 + 1504. / 3 * Hm1m101 * x + 352 * Hm1m101 * x2
                 + 1536. / 5 * Hm1m101 * x3 - 64 * Hm1m1010 - 128 * Hm1m1010 * x
                 - 128 * Hm1m1010 * x2 - 64 * Hm1m1011 - 128 * Hm1m1011 * x
                 - 128 * Hm1m1011 * x2 + 32 * Hm1m1001 + 64 * Hm1m1001 * x
                 - 192 * Hm1m1001 * x2 + 224536. / 45 * Hm10
                 + 32. / 5 * Hm10 * zeta2 / x2 - 496. / 3 * Hm10 * zeta2 * x
                 + 120 * Hm10 * zeta2 * x2 + 1152. / 5 * Hm10 * zeta2 * x3
                 + 248. / 3 * Hm10 * zeta2 + 6112. / 3 * Hm10 * zeta3 * x
                 + 4000. / 3 * Hm10 * zeta3 * x2 + 1904. / 3 * Hm10 * zeta3
                 - 1852. / 225 * Hm10 / x2 + 16. / 45 * Hm10 / x
                 + 281416. / 45 * Hm10 * x + 17032. / 15 * Hm10 * x2
                 - 8528. / 25 * Hm10 * x3 + 57512. / 15 * Hm100
                 - 1280 * Hm100 * zeta2 * x - 1152 * Hm100 * zeta2 * x2
                 - 256 * Hm100 * zeta2 - 368. / 45 * Hm100 / x2
                 - 64. / 15 * Hm100 / x + 216856. / 45 * Hm100 * x
                 + 14296. / 15 * Hm100 * x2 - 1472. / 5 * Hm100 * x3
                 + 640 * Hm1000 - 32. / 15 * Hm1000 / x2
                 + 3088. / 3 * Hm1000 * x + 608. / 3 * Hm1000 * x2
                 - 384. / 5 * Hm1000 * x3 - 112. / 3 * Hm10000
                 - 224. / 3 * Hm10000 * x + 160. / 3 * Hm10000 * x2
                 + 31384. / 15 * Hm101 + 128 * Hm101 * zeta2 * x
                 + 128 * Hm101 * zeta2 * x2 + 64 * Hm101 * zeta2
                 - 176. / 45 * Hm101 / x2 - 128. / 15 * Hm101 / x
                 + 180832. / 45 * Hm101 * x + 10184. / 5 * Hm101 * x2
                 - 704. / 5 * Hm101 * x3 + 64 * Hm1010 + 160 * Hm1010 * x
                 + 96 * Hm1010 * x2 - 1376. / 3 * Hm10100
                 - 5056. / 3 * Hm10100 * x - 3904. / 3 * Hm10100 * x2
                 + 64 * Hm1011 + 160 * Hm1011 * x + 96 * Hm1011 * x2
                 - 96 * Hm1001 - 64. / 15 * Hm1001 / x2 + 448. / 3 * Hm1001 * x
                 + 64 * Hm1001 * x2 - 768. / 5 * Hm1001 * x3 + 32 * Hm10010
                 + 64 * Hm10010 * x + 64 * Hm10010 * x2 + 32 * Hm10011
                 + 64 * Hm10011 * x + 64 * Hm10011 * x2 + 384 * Hm10001
                 + 1536 * Hm10001 * x + 1280 * Hm10001 * x2 + 13594. / 225 * H0
                 + 16. / 3 * H0 * zeta2 / x + 14602. / 9 * H0 * zeta2 * x
                 + 13516. / 5 * H0 * zeta2 * x2 - 1408. / 5 * H0 * zeta2 * x3
                 + 1511. / 15 * H0 * zeta2 + 6904. / 15 * H0 * zeta2 * zeta2 * x
                 - 448. / 5 * H0 * zeta2 * zeta2 * x2
                 + 116. / 15 * H0 * zeta2 * zeta2 - 12284. / 3 * H0 * zeta3 * x
                 + 2456. / 3 * H0 * zeta3 * x2 - 4416. / 5 * H0 * zeta3 * x3
                 + 300 * H0 * zeta3 - 48 * H0 * zeta4 * x - 24 * H0 * zeta4
                 + 644. / 75 * H0 / x + 37658. / 75 * H0 * x
                 - 61902. / 25 * H0 * x2 - 872. / 9 * H00
                 - 540 * H00 * zeta2 * x + 8336. / 3 * H00 * zeta2 * x2
                 + 192 * H00 * zeta2 * x3 + 284. / 3 * H00 * zeta2
                 - 760 * H00 * zeta3 * x + 3712. / 3 * H00 * zeta3 * x2
                 + 284 * H00 * zeta3 + 176. / 45 * H00 / x
                 - 36848. / 9 * H00 * x - 5008. / 3 * H00 * x2
                 + 8528. / 25 * H00 * x3 - 178. / 3 * H000
                 - 664. / 3 * H000 * zeta2 * x + 2624. / 3 * H000 * zeta2 * x2
                 + 524. / 3 * H000 * zeta2 + 64. / 15 * H000 / x
                 - 74476. / 45 * H000 * x - 1128 * H000 * x2
                 + 1472. / 5 * H000 * x3 + 44. / 3 * H0000
                 + 244. / 3 * H0000 * x - 3680. / 3 * H0000 * x2
                 + 768. / 5 * H0000 * x3 - 60 * H00000 - 184. / 3 * H00000 * x
                 - 480 * H00000 * x2 + 32906. / 45 * H1
                 + 88. / 45 * H1 * zeta2 / x2 - 5324. / 9 * H1 * zeta2 * x
                 + 1376 * H1 * zeta2 * x2 - 352. / 5 * H1 * zeta2 * x3
                 - 2230. / 3 * H1 * zeta2 - 1312. / 15 * H1 * zeta2 * zeta2 * x
                 - 2336. / 15 * H1 * zeta2 * zeta2 * x2
                 + 656. / 15 * H1 * zeta2 * zeta2 + 272. / 15 * H1 * zeta3 / x2
                 - 3088. / 3 * H1 * zeta3 * x + 2384. / 3 * H1 * zeta3 * x2
                 - 3264. / 5 * H1 * zeta3 * x3 + 4912. / 3 * H1 * zeta3
                 - 608. / 45 * H1 / x + 10906. / 5 * H1 * x
                 - 14918. / 5 * H1 * x2 + 1408. / 3 * H100m10
                 - 2816. / 3 * H100m10 * x + 2432. / 3 * H100m10 * x2
                 + 896 * H10m1 * zeta2 * x - 640 * H10m1 * zeta2 * x2
                 - 448 * H10m1 * zeta2 - 128 * H10m1m10 + 256 * H10m1m10 * x
                 - 256 * H10m1m10 * x2 + 2464. / 3 * H10m10
                 + 64. / 15 * H10m10 / x2 - 2560. / 3 * H10m10 * x
                 + 1216. / 3 * H10m10 * x2 - 768. / 5 * H10m10 * x3
                 + 1472. / 3 * H10m100 - 2944. / 3 * H10m100 * x
                 + 2560. / 3 * H10m100 * x2 + 384 * H10m101 - 768 * H10m101 * x
                 + 512 * H10m101 * x2 - 239. / 3 * H10
                 - 176. / 15 * H10 * zeta2 / x2 - 4336. / 3 * H10 * zeta2 * x
                 + 7928. / 3 * H10 * zeta2 * x2 + 2112. / 5 * H10 * zeta2 * x3
                 - 1276 * H10 * zeta2 - 5792. / 3 * H10 * zeta3 * x
                 + 4256. / 3 * H10 * zeta3 * x2 + 2896. / 3 * H10 * zeta3
                 + 3650. / 3 * H10 * x - 3560. / 3 * H10 * x2
                 - 1066. / 15 * H100 - 928 * H100 * zeta2 * x
                 + 864 * H100 * zeta2 * x2 + 464 * H100 * zeta2
                 + 208. / 15 * H100 / x + 9316. / 15 * H100 * x
                 - 7808. / 15 * H100 * x2 + 92. / 3 * H1000
                 - 32. / 15 * H1000 / x2 + 1952. / 3 * H1000 * x
                 - 1024 * H1000 * x2 + 384. / 5 * H1000 * x3 - 832. / 3 * H10000
                 + 1664. / 3 * H10000 * x - 1280. / 3 * H10000 * x2
                 - 55. / 3 * H11 - 832 * H11 * zeta2 * x
                 + 3520. / 3 * H11 * zeta2 * x2 - 124. / 3 * H11 * zeta2
                 - 5440. / 3 * H11 * zeta3 * x + 3520. / 3 * H11 * zeta3 * x2
                 + 2720. / 3 * H11 * zeta3 + 4198. / 3 * H11 * x
                 - 4348. / 3 * H11 * x2 + 1024. / 3 * H110m10
                 - 2048. / 3 * H110m10 * x + 1280. / 3 * H110m10 * x2
                 - 482 * H110 - 256. / 3 * H110 * zeta2 * x
                 + 1024. / 3 * H110 * zeta2 * x2 + 128. / 3 * H110 * zeta2
                 + 4564. / 3 * H110 * x - 3268. / 3 * H110 * x2 - 1604 * H1100
                 - 208. / 15 * H1100 / x2 + 1792. / 3 * H1100 * x
                 + 760. / 3 * H1100 * x2 + 2496. / 5 * H1100 * x3
                 - 656. / 3 * H11000 + 1312. / 3 * H11000 * x
                 - 928. / 3 * H11000 * x2 - 522 * H111 - 544 * H111 * zeta2 * x
                 + 544 * H111 * zeta2 * x2 + 272 * H111 * zeta2
                 + 4948. / 3 * H111 * x - 1212 * H111 * x2 - 1516. / 3 * H1110
                 + 1488 * H1110 * x - 3800. / 3 * H1110 * x2 - 656 * H11100
                 + 1312 * H11100 * x - 928 * H11100 * x2 - 1588. / 3 * H1111
                 + 4688. / 3 * H1111 * x - 4096. / 3 * H1111 * x2
                 - 1072. / 3 * H11110 + 2144. / 3 * H11110 * x
                 - 2144. / 3 * H11110 * x2 - 400 * H11111 + 800 * H11111 * x
                 - 800 * H11111 * x2 - 1136. / 3 * H11101
                 + 2272. / 3 * H11101 * x - 2272. / 3 * H11101 * x2
                 - 1436. / 3 * H1101 + 1440 * H1101 * x - 3784. / 3 * H1101 * x2
                 - 1024. / 3 * H11010 + 2048. / 3 * H11010 * x
                 - 2048. / 3 * H11010 * x2 - 1328. / 3 * H11011
                 + 2656. / 3 * H11011 * x - 2656. / 3 * H11011 * x2
                 + 64. / 3 * H11001 - 128. / 3 * H11001 * x
                 - 1024. / 3 * H11001 * x2 - 1414. / 3 * H101
                 - 1888. / 3 * H101 * zeta2 * x + 1888. / 3 * H101 * zeta2 * x2
                 + 944. / 3 * H101 * zeta2 + 1476 * H101 * x - 1012 * H101 * x2
                 - 1472. / 3 * H1010 + 4256. / 3 * H1010 * x - 1184 * H1010 * x2
                 - 1312. / 3 * H10100 + 2624. / 3 * H10100 * x
                 - 2048. / 3 * H10100 * x2 - 1652. / 3 * H1011
                 + 4928. / 3 * H1011 * x - 1408 * H1011 * x2
                 - 1184. / 3 * H10110 + 2368. / 3 * H10110 * x
                 - 2368. / 3 * H10110 * x2 - 1408. / 3 * H10111
                 + 2816. / 3 * H10111 * x - 2816. / 3 * H10111 * x2
                 - 432 * H10101 + 864 * H10101 * x - 864 * H10101 * x2
                 + 3788. / 3 * H1001 + 208. / 15 * H1001 / x2
                 + 4384. / 3 * H1001 * x - 7376. / 3 * H1001 * x2
                 - 2496. / 5 * H1001 * x3 - 1168. / 3 * H10010
                 + 2336. / 3 * H10010 * x - 2336. / 3 * H10010 * x2
                 - 496 * H10011 + 992 * H10011 * x - 992 * H10011 * x2
                 - 336 * H10001 + 672 * H10001 * x - 736 * H10001 * x2
                 + 4598. / 45 * H01 - 1576 * H01 * zeta2 * x
                 + 3904. / 3 * H01 * zeta2 * x2 - 700. / 3 * H01 * zeta2
                 - 160. / 3 * H01 * zeta3 * x + 6976. / 3 * H01 * zeta3 * x2
                 + 3152. / 3 * H01 * zeta3 + 608. / 45 * H01 / x
                 + 38141. / 15 * H01 * x - 14444. / 15 * H01 * x2
                 + 384 * H010m10 - 128 * H010m10 * x + 1280. / 3 * H010m10 * x2
                 - 464. / 3 * H010 - 7856. / 3 * H010 * zeta2 * x
                 - 800 * H010 * zeta2 * x2 - 376 * H010 * zeta2
                 + 1592 * H010 * x - 3268. / 3 * H010 * x2 + 140. / 3 * H0100
                 + 2184 * H0100 * x + 760. / 3 * H0100 * x2
                 + 2496. / 5 * H0100 * x3 - 344. / 3 * H01000
                 + 2480. / 3 * H01000 * x - 928. / 3 * H01000 * x2
                 - 628. / 3 * H011 - 784 * H011 * zeta2 * x
                 + 544 * H011 * zeta2 * x2 + 376. / 3 * H011 * zeta2
                 + 5272. / 3 * H011 * x - 1212 * H011 * x2 - 964. / 3 * H0110
                 + 1032 * H0110 * x - 3800. / 3 * H0110 * x2
                 - 2392. / 3 * H01100 - 3152. / 3 * H01100 * x
                 - 6208. / 3 * H01100 * x2 - 364 * H0111 + 3488. / 3 * H0111 * x
                 - 4096. / 3 * H0111 * x2 - 808. / 3 * H01110
                 + 1552. / 3 * H01110 * x - 2176. / 3 * H01110 * x2
                 - 904. / 3 * H01111 + 1808. / 3 * H01111 * x
                 - 800 * H01111 * x2 - 808. / 3 * H01101 + 560 * H01101 * x
                 - 2272. / 3 * H01101 * x2 - 908. / 3 * H0101 + 1064 * H0101 * x
                 - 3784. / 3 * H0101 * x2 - 704. / 3 * H01010 + 448 * H01010 * x
                 - 2048. / 3 * H01010 * x2 - 968. / 3 * H01011
                 + 1936. / 3 * H01011 * x - 2624. / 3 * H01011 * x2
                 + 440 * H01001 + 7024. / 3 * H01001 * x + 800 * H01001 * x2
                 - 1511. / 15 * H001 - 2368. / 3 * H001 * zeta2 * x
                 + 1504. / 3 * H001 * zeta2 * x2 + 160. / 3 * H001 * zeta2
                 - 48. / 5 * H001 / x - 43234. / 45 * H001 * x
                 - 13516. / 5 * H001 * x2 + 704. / 5 * H001 * x3
                 - 532. / 3 * H0010 + 1880. / 3 * H0010 * x - 1280 * H0010 * x2
                 - 224. / 3 * H00100 + 1472. / 3 * H00100 * x
                 - 1568. / 3 * H00100 * x2 - 248 * H0011 + 2536. / 3 * H0011 * x
                 - 1504 * H0011 * x2 - 704. / 3 * H00110
                 + 1408. / 3 * H00110 * x - 800 * H00110 * x2
                 - 872. / 3 * H00111 + 1744. / 3 * H00111 * x
                 - 2816. / 3 * H00111 * x2 - 752. / 3 * H00101
                 + 1504. / 3 * H00101 * x - 864 * H00101 * x2 - 284. / 3 * H0001
                 + 156 * H0001 * x - 8336. / 3 * H0001 * x2
                 - 1728. / 5 * H0001 * x3 - 608. / 3 * H00010
                 + 1024. / 3 * H00010 * x - 2528. / 3 * H00010 * x2
                 - 808. / 3 * H00011 + 1424. / 3 * H00011 * x
                 - 3136. / 3 * H00011 * x2 - 524. / 3 * H00001
                 + 664. / 3 * H00001 * x - 2624. / 3 * H00001 * x2)
        + nf * CA * CA
              * (+271171. / 486 + 1448 * zeta2 * zeta3 * x
                 - 688. / 3 * zeta2 * zeta3 * x2 + 644. / 3 * zeta2 * zeta3
                 + 80588. / 405 * zeta2 / x - 3029684. / 405 * zeta2 * x
                 + 2479774. / 405 * zeta2 * x2 + 228. / 5 * zeta2 * x3
                 - 1628. / 3 * zeta2 + 4544. / 45 * zeta2 * zeta2 / x
                 + 37108. / 45 * zeta2 * zeta2 * x
                 - 2272. / 3 * zeta2 * zeta2 * x2
                 - 5184. / 25 * zeta2 * zeta2 * x3 + 958. / 15 * zeta2 * zeta2
                 + 1904. / 27 * zeta3 / x - 30292. / 9 * zeta3 * x
                 + 82799. / 15 * zeta3 * x2 + 60 * zeta3 * x3
                 - 5302. / 15 * zeta3 - 32 * zeta4 / x - 96 * zeta4 * x
                 + 248 * zeta4 * x2 - 84 * zeta4 - 1624 * zeta5 * x
                 - 1744. / 3 * zeta5 * x2 - 2612. / 3 * zeta5
                 - 5017249. / 3645 / x + 2812798. / 243 * x
                 - 79819747. / 7290 * x2 - 328. / 3 * H000m10
                 + 400. / 3 * H000m10 * x + 256. / 3 * H000m10 * x2
                 - 280. / 3 * H00m1 * zeta2 * x - 192 * H00m1 * zeta2 * x2
                 + 92. / 3 * H00m1 * zeta2 - 104 * H00m1m10
                 + 1520. / 3 * H00m1m10 * x - 640. / 3 * H00m1m10 * x2
                 + 384 * H00m10 + 8456. / 9 * H00m10 * x2
                 - 384. / 5 * H00m10 * x3 - 272. / 3 * H00m100
                 + 208. / 3 * H00m100 * x + 448. / 3 * H00m100 * x2
                 - 248. / 3 * H00m101 + 1040. / 3 * H00m101 * x
                 + 256. / 3 * H00m101 * x2 + 64 * H0m1 * zeta2 / x
                 + 2552. / 3 * H0m1 * zeta2 * x - 3676. / 3 * H0m1 * zeta2 * x2
                 + 768. / 5 * H0m1 * zeta2 * x3 - 160. / 3 * H0m1 * zeta2
                 + 1016. / 3 * H0m1 * zeta3 * x - 112. / 3 * H0m1 * zeta3 * x2
                 + 676. / 3 * H0m1 * zeta3 - 184. / 3 * H0m10m10
                 + 1712. / 3 * H0m10m10 * x - 32. / 3 * H0m10m10 * x2
                 - 640 * H0m1m1 * zeta2 * x - 256. / 3 * H0m1m1 * zeta2 * x2
                 - 288 * H0m1m1 * zeta2 + 16 * H0m1m1m10 - 800 * H0m1m1m10 * x
                 - 64 * H0m1m1m10 * x2 - 352. / 3 * H0m1m10
                 + 128. / 3 * H0m1m10 / x + 2176. / 3 * H0m1m10 * x
                 - 1208. / 3 * H0m1m10 * x2 + 172. / 3 * H0m1m100
                 + 1864. / 3 * H0m1m100 * x - 368. / 3 * H0m1m100 * x2
                 + 296 * H0m1m101 + 240 * H0m1m101 * x
                 + 160. / 3 * H0m1m101 * x2 - 812. / 5 * H0m10
                 - 272 * H0m10 * zeta2 * x - 496. / 3 * H0m10 * zeta2 * x2
                 + 824. / 3 * H0m10 * zeta2 + 1984. / 45 * H0m10 / x
                 + 2992. / 5 * H0m10 * x + 19118. / 135 * H0m10 * x2
                 + 24 * H0m10 * x3 + 824. / 3 * H0m100 - 256. / 3 * H0m100 / x
                 - 644. / 3 * H0m100 * x + 14092. / 9 * H0m100 * x2
                 - 384. / 5 * H0m100 * x3 - 440. / 3 * H0m1000
                 + 432 * H0m1000 * x + 512. / 3 * H0m1000 * x2
                 - 16. / 3 * H0m101 - 128. / 3 * H0m101 / x - 488 * H0m101 * x
                 + 1024 * H0m101 * x2 - 768. / 5 * H0m101 * x3 - 24 * H0m1010
                 + 304 * H0m1010 * x + 416. / 3 * H0m1010 * x2
                 - 32. / 3 * H0m1011 + 1088. / 3 * H0m1011 * x
                 + 512. / 3 * H0m1011 * x2 - 644. / 3 * H0m1001
                 + 1640. / 3 * H0m1001 * x + 784. / 3 * H0m1001 * x2
                 - Hm1 * zeta2 / x2 + 42976. / 135 * Hm1 * zeta2 / x
                 - 15512. / 15 * Hm1 * zeta2 * x
                 - 105449. / 135 * Hm1 * zeta2 * x2 - 36 * Hm1 * zeta2 * x3
                 - 8836. / 45 * Hm1 * zeta2
                 + 3664. / 15 * Hm1 * zeta2 * zeta2 * x
                 + 2968. / 15 * Hm1 * zeta2 * zeta2 * x2
                 + 136. / 3 * Hm1 * zeta2 * zeta2 + 16. / 5 * Hm1 * zeta3 / x2
                 + 48 * Hm1 * zeta3 / x + 772. / 3 * Hm1 * zeta3 * x
                 - 532 * Hm1 * zeta3 * x2 + 576. / 5 * Hm1 * zeta3 * x3
                 + 2120. / 3 * Hm1 * zeta3 + 48 * Hm100m10 + 96 * Hm100m10 * x
                 - 64 * Hm100m10 * x2 - 2144. / 3 * Hm10m1 * zeta2 * x
                 - 704. / 3 * Hm10m1 * zeta2 * x2 - 1072. / 3 * Hm10m1 * zeta2
                 - 512. / 3 * Hm10m1m10 - 1024. / 3 * Hm10m1m10 * x
                 - 448. / 3 * Hm10m1m10 * x2 - 112. / 3 * Hm10m10
                 + 448. / 9 * Hm10m10 / x + 472. / 3 * Hm10m10 * x
                 + 1504. / 9 * Hm10m10 * x2 + 808. / 3 * Hm10m100
                 + 1616. / 3 * Hm10m100 * x + 368. / 3 * Hm10m100 * x2
                 + 272 * Hm10m101 + 544 * Hm10m101 * x + 160 * Hm10m101 * x2
                 - 64. / 15 * Hm1m1 * zeta2 / x2 - 896. / 9 * Hm1m1 * zeta2 / x
                 - 2192. / 3 * Hm1m1 * zeta2 * x
                 + 1876. / 9 * Hm1m1 * zeta2 * x2
                 - 768. / 5 * Hm1m1 * zeta2 * x3 - 936 * Hm1m1 * zeta2
                 - 2912. / 3 * Hm1m1 * zeta3 * x
                 - 1040. / 3 * Hm1m1 * zeta3 * x2 - 1456. / 3 * Hm1m1 * zeta3
                 - 512. / 3 * Hm1m10m10 - 1024. / 3 * Hm1m10m10 * x
                 - 448. / 3 * Hm1m10m10 * x2 + 1328 * Hm1m1m1 * zeta2 * x
                 + 592 * Hm1m1m1 * zeta2 * x2 + 664 * Hm1m1m1 * zeta2
                 + 208 * Hm1m1m1m10 + 416 * Hm1m1m1m10 * x
                 + 224 * Hm1m1m1m10 * x2 + 80 * Hm1m1m10
                 - 512. / 9 * Hm1m1m10 / x - 304 * Hm1m1m10 * x
                 - 3416. / 9 * Hm1m1m10 * x2 - 1600. / 3 * Hm1m1m100
                 - 3200. / 3 * Hm1m1m100 * x - 1184. / 3 * Hm1m1m100 * x2
                 - 560 * Hm1m1m101 - 1120 * Hm1m1m101 * x - 480 * Hm1m1m101 * x2
                 - 3464. / 9 * Hm1m10 - 672 * Hm1m10 * zeta2 * x
                 - 32 * Hm1m10 * zeta2 * x2 - 336 * Hm1m10 * zeta2
                 - 2. / 3 * Hm1m10 / x2 + 3584. / 27 * Hm1m10 / x
                 - 256. / 3 * Hm1m10 * x + 5618. / 27 * Hm1m10 * x2
                 - 24 * Hm1m10 * x3 + 528 * Hm1m100 + 32. / 15 * Hm1m100 / x2
                 + 256. / 3 * Hm1m100 / x + 580. / 3 * Hm1m100 * x
                 - 1412. / 3 * Hm1m100 * x2 + 384. / 5 * Hm1m100 * x3
                 + 80 * Hm1m1000 + 160 * Hm1m1000 * x - 224 * Hm1m1000 * x2
                 + 976 * Hm1m101 + 64. / 15 * Hm1m101 / x2
                 + 640. / 9 * Hm1m101 / x + 1736. / 3 * Hm1m101 * x
                 - 3584. / 9 * Hm1m101 * x2 + 768. / 5 * Hm1m101 * x3
                 - 80 * Hm1m1010 - 160 * Hm1m1010 * x - 224 * Hm1m1010 * x2
                 - 256. / 3 * Hm1m1011 - 512. / 3 * Hm1m1011 * x
                 - 704. / 3 * Hm1m1011 * x2 + 536. / 3 * Hm1m1001
                 + 1072. / 3 * Hm1m1001 * x - 368. / 3 * Hm1m1001 * x2
                 + 487361. / 405 * Hm10 + 16. / 5 * Hm10 * zeta2 / x2
                 - 224. / 9 * Hm10 * zeta2 / x - 424 * Hm10 * zeta2 * x
                 - 12164. / 9 * Hm10 * zeta2 * x2 + 576. / 5 * Hm10 * zeta2 * x3
                 + 748 * Hm10 * zeta2 + 1904. / 3 * Hm10 * zeta3 * x
                 + 608. / 3 * Hm10 * zeta3 * x2 + 664. / 3 * Hm10 * zeta3
                 + 19. / 15 * Hm10 / x2 + 15646. / 45 * Hm10 / x
                 + 740146. / 405 * Hm10 * x + 487654. / 405 * Hm10 * x2
                 + 228. / 5 * Hm10 * x3 - 5516. / 135 * Hm100
                 - 320 * Hm100 * zeta2 * x - 544 * Hm100 * zeta2 * x2
                 - 64 * Hm100 * zeta2 + 2. / 3 * Hm100 / x2
                 - 67088. / 135 * Hm100 / x + 148964. / 135 * Hm100 * x
                 + 44954. / 45 * Hm100 * x2 + 24 * Hm100 * x3
                 - 1964. / 9 * Hm1000 - 16. / 15 * Hm1000 / x2
                 + 224. / 3 * Hm1000 / x + 9248. / 9 * Hm1000 * x
                 + 13784. / 9 * Hm1000 * x2 - 192. / 5 * Hm1000 * x3
                 + 32. / 3 * Hm10000 + 64. / 3 * Hm10000 * x
                 + 544. / 3 * Hm10000 * x2 + 176. / 45 * Hm101
                 - 48 * Hm101 * zeta2 * x - 16 * Hm101 * zeta2 * x2
                 - 24 * Hm101 * zeta2 + 2. / 3 * Hm101 / x2
                 - 34016. / 135 * Hm101 / x + 14872. / 15 * Hm101 * x
                 + 119494. / 135 * Hm101 * x2 + 24 * Hm101 * x3 - 16 * Hm1010
                 + 128. / 3 * Hm1010 / x + 1384. / 3 * Hm1010 * x
                 + 1744. / 3 * Hm1010 * x2 - 80 * Hm10100 - 352 * Hm10100 * x
                 - 256 * Hm10100 * x2 + 8 * Hm1011 + 512. / 9 * Hm1011 / x
                 + 528 * Hm1011 * x + 5744. / 9 * Hm1011 * x2
                 + 112. / 3 * Hm10110 + 224. / 3 * Hm10110 * x
                 + 224. / 3 * Hm10110 * x2 + 32. / 3 * Hm10111
                 + 64. / 3 * Hm10111 * x + 64. / 3 * Hm10111 * x2
                 + 80. / 3 * Hm10101 + 160. / 3 * Hm10101 * x
                 + 160. / 3 * Hm10101 * x2 - 1796. / 3 * Hm1001
                 - 32. / 15 * Hm1001 / x2 + 128. / 3 * Hm1001 / x
                 + 628 * Hm1001 * x + 4376. / 3 * Hm1001 * x2
                 - 384. / 5 * Hm1001 * x3 + 112 * Hm10010 + 224 * Hm10010 * x
                 + 256 * Hm10010 * x2 + 400. / 3 * Hm10011
                 + 800. / 3 * Hm10011 * x + 896. / 3 * Hm10011 * x2
                 + 440. / 3 * Hm10001 + 1456. / 3 * Hm10001 * x
                 + 1840. / 3 * Hm10001 * x2 - 3220177. / 1215 * H0
                 + 2128. / 45 * H0 * zeta2 / x - 242624. / 45 * H0 * zeta2 * x
                 + 881432. / 135 * H0 * zeta2 * x2 + 48 * H0 * zeta2 * x3
                 + 33848. / 45 * H0 * zeta2 + 2224. / 3 * H0 * zeta2 * zeta2 * x
                 - 1552. / 15 * H0 * zeta2 * zeta2 * x2
                 - 784. / 15 * H0 * zeta2 * zeta2 - 128. / 9 * H0 * zeta3 / x
                 - 53144. / 9 * H0 * zeta3 * x + 19496. / 9 * H0 * zeta3 * x2
                 - 384 * H0 * zeta3 * x3 - 416. / 9 * H0 * zeta3
                 - 192 * H0 * zeta4 * x - 48 * H0 * zeta4 * x2 - 72 * H0 * zeta4
                 - 195577. / 1215 * H0 / x + 13320758. / 1215 * H0 * x
                 - 13259914. / 1215 * H0 * x2 + 423067. / 405 * H00
                 - 12932. / 3 * H00 * zeta2 * x + 14704. / 9 * H00 * zeta2 * x2
                 + 192. / 5 * H00 * zeta2 * x3 - 32 * H00 * zeta2
                 - 4432. / 3 * H00 * zeta3 * x - 1136. / 3 * H00 * zeta3
                 - 14. / 5 * H00 / x + 1068164. / 135 * H00 * x
                 - 689798. / 135 * H00 * x2 - 228. / 5 * H00 * x3
                 - 32618. / 27 * H000 - 1312 * H000 * zeta2 * x
                 + 192 * H000 * zeta2 * x2 + 56 * H000 * zeta2
                 + 32. / 15 * H000 / x + 743056. / 135 * H000 * x
                 - 152792. / 27 * H000 * x2 - 24 * H000 * x3 + 1436. / 9 * H0000
                 + 11224. / 3 * H0000 * x - 360 * H0000 * x2
                 + 384. / 5 * H0000 * x3 - 240 * H00000 + 1056 * H00000 * x
                 - 311806. / 243 * H1 - 1. / 3 * H1 * zeta2 / x2
                 - 3704. / 27 * H1 * zeta2 / x - 94852. / 27 * H1 * zeta2 * x
                 + 34045. / 9 * H1 * zeta2 * x2 + 12 * H1 * zeta2 * x3
                 - 2110. / 27 * H1 * zeta2 - 32. / 5 * H1 * zeta2 * zeta2 * x
                 - 136. / 5 * H1 * zeta2 * zeta2 * x2
                 + 16. / 5 * H1 * zeta2 * zeta2 + 112. / 15 * H1 * zeta3 / x2
                 - 944. / 9 * H1 * zeta3 / x - 6152. / 3 * H1 * zeta3 * x
                 + 16340. / 9 * H1 * zeta3 * x2 - 1344. / 5 * H1 * zeta3 * x3
                 + 2444. / 3 * H1 * zeta3 - 48 * H1 * zeta4 * x
                 + 48 * H1 * zeta4 * x2 + 24 * H1 * zeta4
                 + 306254. / 405 * H1 / x + 2653214. / 243 * H1 * x
                 - 12973432. / 1215 * H1 * x2 + 352. / 3 * H100m10
                 - 704. / 3 * H100m10 * x + 608. / 3 * H100m10 * x2
                 + 208 * H10m1 * zeta2 * x - 112 * H10m1 * zeta2 * x2
                 - 104 * H10m1 * zeta2 - 16. / 3 * H10m1m10
                 + 32. / 3 * H10m1m10 * x - 224. / 3 * H10m1m10 * x2
                 + 336 * H10m10 + 32. / 15 * H10m10 / x2 + 128. / 9 * H10m10 / x
                 - 752. / 3 * H10m10 * x + 424. / 9 * H10m10 * x2
                 - 384. / 5 * H10m10 * x3 + 96 * H10m100 - 192 * H10m100 * x
                 + 160 * H10m100 * x2 + 304. / 3 * H10m101
                 - 608. / 3 * H10m101 * x + 224. / 3 * H10m101 * x2
                 + 602. / 27 * H10 - 64. / 15 * H10 * zeta2 / x2
                 - 352. / 3 * H10 * zeta2 / x - 11176. / 9 * H10 * zeta2 * x
                 + 16180. / 9 * H10 * zeta2 * x2 + 768. / 5 * H10 * zeta2 * x3
                 - 4468. / 9 * H10 * zeta2 - 3728. / 3 * H10 * zeta3 * x
                 + 2480. / 3 * H10 * zeta3 * x2 + 1864. / 3 * H10 * zeta3
                 - 25084. / 81 * H10 / x + 16262. / 3 * H10 * x
                 - 435044. / 81 * H10 * x2 + 3758. / 27 * H100
                 - 2240. / 3 * H100 * zeta2 * x + 1520. / 3 * H100 * zeta2 * x2
                 + 1120. / 3 * H100 * zeta2 + 136. / 3 * H100 / x
                 + 91832. / 27 * H100 * x - 100910. / 27 * H100 * x2
                 + 64. / 3 * H1000 - 16. / 15 * H1000 / x2
                 + 1088. / 9 * H1000 / x + 3584. / 3 * H1000 * x
                 - 13256. / 9 * H1000 * x2 + 192. / 5 * H1000 * x3
                 - 688. / 3 * H10000 + 1376. / 3 * H10000 * x
                 - 896. / 3 * H10000 * x2 - 22454. / 81 * H11
                 - 896. / 9 * H11 * zeta2 / x - 11672. / 9 * H11 * zeta2 * x
                 + 13048. / 9 * H11 * zeta2 * x2 + 472. / 9 * H11 * zeta2
                 - 1312 * H11 * zeta3 * x + 848 * H11 * zeta3 * x2
                 + 656 * H11 * zeta3 - 332 * H11 / x + 507886. / 81 * H11 * x
                 - 486118. / 81 * H11 * x2 + 144 * H110m10 - 288 * H110m10 * x
                 + 160 * H110m10 * x2 - 2314. / 27 * H110
                 - 1456. / 3 * H110 * zeta2 * x + 1168. / 3 * H110 * zeta2 * x2
                 + 728. / 3 * H110 * zeta2 + 968. / 27 * H110 / x
                 + 100208. / 27 * H110 * x - 107626. / 27 * H110 * x2
                 - 5884. / 9 * H1100 - 16. / 3 * H1100 / x2
                 + 896. / 9 * H1100 / x + 14156. / 9 * H1100 * x
                 - 12388. / 9 * H1100 * x2 + 192 * H1100 * x3
                 - 896. / 3 * H11000 + 1792. / 3 * H11000 * x
                 - 1216. / 3 * H11000 * x2 - 4300. / 27 * H111
                 - 1184. / 3 * H111 * zeta2 * x + 896. / 3 * H111 * zeta2 * x2
                 + 592. / 3 * H111 * zeta2 + 2168. / 27 * H111 / x
                 + 94736. / 27 * H111 * x - 98716. / 27 * H111 * x2
                 - 1028. / 9 * H1110 + 640. / 9 * H1110 / x
                 + 10456. / 9 * H1110 * x - 11648. / 9 * H1110 * x2
                 - 784. / 3 * H11100 + 1568. / 3 * H11100 * x
                 - 1184. / 3 * H11100 * x2 - 860. / 9 * H1111
                 + 416. / 9 * H1111 / x + 8488. / 9 * H1111 * x
                 - 1024 * H1111 * x2 - 248. / 3 * H11110 + 496. / 3 * H11110 * x
                 - 496. / 3 * H11110 * x2 - 80 * H11111 + 160 * H11111 * x
                 - 160 * H11111 * x2 - 280. / 3 * H11101 + 560. / 3 * H11101 * x
                 - 560. / 3 * H11101 * x2 - 832. / 9 * H1101
                 + 640. / 9 * H1101 / x + 10304. / 9 * H1101 * x
                 - 1260 * H1101 * x2 - 344. / 3 * H11010 + 688. / 3 * H11010 * x
                 - 688. / 3 * H11010 * x2 - 304. / 3 * H11011
                 + 608. / 3 * H11011 * x - 608. / 3 * H11011 * x2
                 - 256. / 3 * H11001 + 512. / 3 * H11001 * x
                 - 704. / 3 * H11001 * x2 - 3086. / 27 * H101
                 - 1072. / 3 * H101 * zeta2 * x + 784. / 3 * H101 * zeta2 * x2
                 + 536. / 3 * H101 * zeta2 + 1912. / 27 * H101 / x
                 + 96004. / 27 * H101 * x - 99326. / 27 * H101 * x2
                 - 412. / 9 * H1010 + 256. / 3 * H1010 / x
                 + 9836. / 9 * H1010 * x - 11432. / 9 * H1010 * x2
                 - 488. / 3 * H10100 + 976. / 3 * H10100 * x
                 - 736. / 3 * H10100 * x2 - 520. / 9 * H1011 + 64 * H1011 / x
                 + 9440. / 9 * H1011 * x - 10460. / 9 * H1011 * x2
                 - 352. / 3 * H10110 + 704. / 3 * H10110 * x
                 - 704. / 3 * H10110 * x2 - 320. / 3 * H10111
                 + 640. / 3 * H10111 * x - 640. / 3 * H10111 * x2
                 - 280. / 3 * H10101 + 560. / 3 * H10101 * x
                 - 560. / 3 * H10101 * x2 + 5812. / 9 * H1001
                 + 16. / 3 * H1001 / x2 + 896. / 9 * H1001 / x
                 + 9340. / 9 * H1001 * x - 5072. / 3 * H1001 * x2
                 - 192 * H1001 * x3 - 448. / 3 * H10010 + 896. / 3 * H10010 * x
                 - 800. / 3 * H10010 * x2 - 584. / 3 * H10011
                 + 1168. / 3 * H10011 * x - 1072. / 3 * H10011 * x2
                 - 872. / 3 * H10001 + 1744. / 3 * H10001 * x
                 - 1312. / 3 * H10001 * x2 + 1628. / 3 * H01
                 - 224. / 3 * H01 * zeta2 / x - 8528. / 3 * H01 * zeta2 * x
                 + 3280. / 3 * H01 * zeta2 * x2 - 136. / 3 * H01 * zeta2
                 - 3440. / 3 * H01 * zeta3 * x + 992 * H01 * zeta3 * x2
                 + 1240. / 3 * H01 * zeta3 + 60226. / 405 * H01 / x
                 + 83774. / 9 * H01 * x - 2479774. / 405 * H01 * x2
                 + 544. / 3 * H010m10 + 64. / 3 * H010m10 * x
                 + 608. / 3 * H010m10 * x2 - 964. / 3 * H010
                 - 5360. / 3 * H010 * zeta2 * x - 176. / 3 * H010 * zeta2 * x2
                 - 712. / 3 * H010 * zeta2 - 832. / 9 * H010 / x
                 + 43292. / 9 * H010 * x - 4906 * H010 * x2 + 72 * H0100
                 + 128. / 3 * H0100 / x + 9628. / 3 * H0100 * x
                 - 9212. / 9 * H0100 * x2 + 192 * H0100 * x3 + 464. / 3 * H01000
                 + 3616. / 3 * H01000 * x - 352. / 3 * H01000 * x2
                 - 4436. / 9 * H011 - 3152. / 3 * H011 * zeta2 * x
                 + 736. / 3 * H011 * zeta2 * x2 - 8 * H011 * zeta2
                 - 3280. / 27 * H011 / x + 46388. / 9 * H011 * x
                 - 140176. / 27 * H011 * x2 - 24 * H0110 + 160. / 3 * H0110 / x
                 + 2672 * H0110 * x - 15080. / 9 * H0110 * x2
                 - 652. / 3 * H01100 + 440 * H01100 * x
                 - 1856. / 3 * H01100 * x2 + 512. / 9 * H0111 / x
                 + 7264. / 3 * H0111 * x - 10808. / 9 * H0111 * x2 - 8 * H01110
                 + 1936. / 3 * H01110 * x - 688. / 3 * H01110 * x2
                 - 80. / 3 * H01111 + 1408. / 3 * H01111 * x - 192 * H01111 * x2
                 + 1952. / 3 * H01101 * x - 640. / 3 * H01101 * x2
                 - 40. / 3 * H0101 + 160. / 3 * H0101 / x + 2480 * H0101 * x
                 - 3884. / 3 * H0101 * x2 + 124. / 3 * H01010
                 + 2072. / 3 * H01010 * x - 176 * H01010 * x2 - 56. / 3 * H01011
                 + 1840. / 3 * H01011 * x - 688. / 3 * H01011 * x2
                 + 892. / 3 * H01001 + 1512 * H01001 * x
                 + 464. / 3 * H01001 * x2 - 33848. / 45 * H001
                 - 1368 * H001 * zeta2 * x + 320. / 3 * H001 * zeta2 * x2
                 - 244. / 3 * H001 * zeta2 - 16. / 5 * H001 / x
                 + 269552. / 45 * H001 * x - 881432. / 135 * H001 * x2
                 - 24 * H001 * x3 - 4 * H0010 + 3316 * H0010 * x
                 - 1296 * H0010 * x2 + 460. / 3 * H00100
                 + 3544. / 3 * H00100 * x - 32 * H00100 * x2 - 100. / 3 * H0011
                 + 11140. / 3 * H0011 * x - 13540. / 9 * H0011 * x2
                 + 296. / 3 * H00110 + 3568. / 3 * H00110 * x
                 - 512. / 3 * H00110 * x2 + 32 * H00111 + 3200. / 3 * H00111 * x
                 - 192 * H00111 * x2 + 88. / 3 * H00101 + 3344. / 3 * H00101 * x
                 - 640. / 3 * H00101 * x2 + 32 * H0001 + 12932. / 3 * H0001 * x
                 - 14704. / 9 * H0001 * x2 - 576. / 5 * H0001 * x3 + 36 * H00010
                 + 1256 * H00010 * x - 160 * H00010 * x2 + 12 * H00011
                 + 1432 * H00011 * x - 288 * H00011 * x2 - 56 * H00001
                 + 4336. / 3 * H00001 * x - 192 * H00001 * x2)
        + nf * nf * CF
              * (+54984209. / 24300 - 64 * zeta2 * zeta3 * x
                 + 32 * zeta2 * zeta3 - 2656. / 135 * zeta2 / x
                 + 90604. / 405 * zeta2 * x - 82984. / 405 * zeta2 * x2
                 - 3424. / 75 * zeta2 * x3 - 282608. / 405 * zeta2
                 - 6752. / 45 * zeta2 * zeta2 * x
                 - 1264. / 45 * zeta2 * zeta2 * x2 - 908. / 45 * zeta2 * zeta2
                 + 128. / 27 * zeta3 / x + 18364. / 27 * zeta3 * x
                 - 1112. / 3 * zeta3 * x2 + 64 * zeta3 * x3
                 - 10436. / 27 * zeta3 + 48 * zeta5 * x - 24 * zeta5
                 + 1100288. / 18225 / x - 39073417. / 12150 * x
                 + 16071037. / 18225 * x2 - 160 * H00m10 + 64 * H00m10 * x
                 - 128 * H00m10 * x2 - 224. / 3 * H0m1 * zeta2 * x
                 + 224. / 3 * H0m1 * zeta2 * x2 + 272. / 3 * H0m1 * zeta2
                 + 224. / 3 * H0m1m10 - 64. / 3 * H0m1m10 * x
                 + 64 * H0m1m10 * x2 - 2480. / 9 * H0m10 - 32. / 45 * H0m10 / x2
                 + 2176. / 9 * H0m10 * x - 176. / 9 * H0m10 * x2
                 - 64. / 5 * H0m10 * x3 - 512. / 3 * H0m100
                 + 512. / 3 * H0m100 * x - 416. / 3 * H0m100 * x2
                 - 160. / 3 * H0m101 + 64 * H0m101 * x - 128. / 3 * H0m101 * x2
                 + 8. / 15 * Hm1 * zeta2 / x2 - 32. / 3 * Hm1 * zeta2 / x
                 + 800. / 9 * Hm1 * zeta2 * x - 440. / 9 * Hm1 * zeta2 * x2
                 + 96. / 5 * Hm1 * zeta2 * x3 + 1384. / 9 * Hm1 * zeta2
                 + 224. / 3 * Hm1 * zeta3 * x + 112. / 3 * Hm1 * zeta3 * x2
                 + 112. / 3 * Hm1 * zeta3 + 128. / 3 * Hm10m10
                 + 256. / 3 * Hm10m10 * x + 128. / 3 * Hm10m10 * x2
                 - 256. / 3 * Hm1m1 * zeta2 * x - 128. / 3 * Hm1m1 * zeta2 * x2
                 - 128. / 3 * Hm1m1 * zeta2 - 128. / 3 * Hm1m1m10
                 - 256. / 3 * Hm1m1m10 * x - 128. / 3 * Hm1m1m10 * x2
                 + 1328. / 9 * Hm1m10 + 16. / 45 * Hm1m10 / x2
                 - 64. / 9 * Hm1m10 / x + 448. / 3 * Hm1m10 * x
                 + 112. / 9 * Hm1m10 * x2 + 64. / 5 * Hm1m10 * x3
                 + 224. / 3 * Hm1m100 + 448. / 3 * Hm1m100 * x
                 + 224. / 3 * Hm1m100 * x2 + 64. / 3 * Hm1m101
                 + 128. / 3 * Hm1m101 * x + 64. / 3 * Hm1m101 * x2
                 - 23312. / 45 * Hm10 + 128. / 3 * Hm10 * zeta2 * x
                 + 64. / 3 * Hm10 * zeta2 * x2 + 64. / 3 * Hm10 * zeta2
                 - 1016. / 675 * Hm10 / x2 - 2608. / 135 * Hm10 / x
                 - 38176. / 135 * Hm10 * x + 20912. / 135 * Hm10 * x2
                 - 3424. / 75 * Hm10 * x3 - 2464. / 9 * Hm100
                 - 16. / 15 * Hm100 / x2 + 64. / 3 * Hm100 / x
                 - 992. / 9 * Hm100 * x + 1184. / 9 * Hm100 * x2
                 - 192. / 5 * Hm100 * x3 - 32 * Hm1000 - 64 * Hm1000 * x
                 - 32 * Hm1000 * x2 - 80 * Hm101 - 16. / 45 * Hm101 / x2
                 + 64. / 9 * Hm101 / x - 128. / 9 * Hm101 * x
                 + 496. / 9 * Hm101 * x2 - 64. / 5 * Hm101 * x3
                 - 32. / 3 * Hm1001 - 64. / 3 * Hm1001 * x
                 - 32. / 3 * Hm1001 * x2 + 10607848. / 6075 * H0
                 + 5524. / 9 * H0 * zeta2 * x - 304. / 3 * H0 * zeta2 * x2
                 - 288. / 5 * H0 * zeta2 * x3 - 4408. / 9 * H0 * zeta2
                 - 112. / 5 * H0 * zeta2 * zeta2 * x
                 + 56. / 5 * H0 * zeta2 * zeta2 + 3376. / 9 * H0 * zeta3 * x
                 - 64. / 3 * H0 * zeta3 * x2 - 1268. / 9 * H0 * zeta3
                 + 776. / 675 * H0 / x - 8642522. / 6075 * H0 * x
                 + 1049188. / 2025 * H0 * x2 + 385274. / 405 * H00
                 + 1520. / 3 * H00 * zeta2 * x + 160. / 3 * H00 * zeta2 * x2
                 - 436. / 3 * H00 * zeta2 + 304. / 3 * H00 * zeta3 * x
                 - 152. / 3 * H00 * zeta3 + 16. / 15 * H00 / x
                 - 277906. / 405 * H00 * x - 28088. / 405 * H00 * x2
                 + 3424. / 75 * H00 * x3 + 17782. / 27 * H000
                 + 176 * H000 * zeta2 * x - 88 * H000 * zeta2
                 - 27200. / 27 * H000 * x + 64. / 9 * H000 * x2
                 + 192. / 5 * H000 * x3 + 1636. / 9 * H0000
                 - 7232. / 9 * H0000 * x - 928. / 9 * H0000 * x2 + 120 * H00000
                 - 240 * H00000 * x + 1348858. / 1215 * H1
                 + 8. / 45 * H1 * zeta2 / x2 + 32. / 3 * H1 * zeta2 / x
                 + 216 * H1 * zeta2 * x - 1544. / 9 * H1 * zeta2 * x2
                 - 192. / 5 * H1 * zeta2 * x3 - 608. / 9 * H1 * zeta2
                 + 224. / 9 * H1 * zeta3 * x - 272. / 9 * H1 * zeta3 * x2
                 - 112. / 9 * H1 * zeta3 - 13904. / 405 * H1 / x
                 - 1830218. / 1215 * H1 * x + 507692. / 1215 * H1 * x2
                 - 64. / 3 * H10m10 + 128. / 3 * H10m10 * x
                 - 64. / 3 * H10m10 * x2 + 33758. / 81 * H10
                 + 64. / 3 * H10 * zeta2 * x - 160. / 3 * H10 * zeta2 * x2
                 - 32. / 3 * H10 * zeta2 + 512. / 27 * H10 / x
                 - 58504. / 81 * H10 * x + 26728. / 81 * H10 * x2
                 + 3080. / 27 * H100 - 64. / 9 * H100 / x
                 - 6472. / 27 * H100 * x + 5344. / 27 * H100 * x2
                 + 80. / 9 * H1000 - 160. / 9 * H1000 * x
                 + 448. / 9 * H1000 * x2 + 33766. / 81 * H11
                 - 64. / 3 * H11 * zeta2 * x - 32. / 3 * H11 * zeta2 * x2
                 + 32. / 3 * H11 * zeta2 + 512. / 27 * H11 / x
                 - 55556. / 81 * H11 * x + 23384. / 81 * H11 * x2
                 + 1048. / 9 * H110 - 64. / 9 * H110 / x - 1976. / 9 * H110 * x
                 + 544. / 3 * H110 * x2 - 32 * H110 * x3 + 16. / 3 * H1100
                 - 32. / 3 * H1100 * x + 64. / 3 * H1100 * x2
                 + 3224. / 27 * H111 - 64. / 9 * H111 / x
                 - 6184. / 27 * H111 * x + 4144. / 27 * H111 * x2
                 + 64. / 3 * H1110 - 128. / 3 * H1110 * x + 32 * H1110 * x2
                 + 136. / 9 * H1111 - 272. / 9 * H1111 * x
                 + 272. / 9 * H1111 * x2 + 32. / 3 * H1101 - 64. / 3 * H1101 * x
                 + 32 * H1101 * x2 + 424. / 3 * H101 - 64. / 9 * H101 / x
                 - 872. / 3 * H101 * x + 1600. / 9 * H101 * x2 + 32 * H101 * x3
                 + 32 * H1010 - 64 * H1010 * x + 160. / 3 * H1010 * x2
                 + 80. / 3 * H1011 - 160. / 3 * H1011 * x
                 + 160. / 3 * H1011 * x2 + 64. / 3 * H1001
                 - 128. / 3 * H1001 * x + 64 * H1001 * x2 + 282608. / 405 * H01
                 + 32 * H01 * zeta2 * x + 64. / 3 * H01 * zeta2 * x2
                 + 8. / 3 * H01 * zeta2 + 16. / 45 * H01 / x
                 - 205132. / 405 * H01 * x + 82984. / 405 * H01 * x2
                 + 812. / 3 * H010 - 400. / 3 * H010 * x + 1568. / 9 * H010 * x2
                 - 32 * H010 * x3 + 24 * H0100 + 7388. / 27 * H011
                 - 3856. / 27 * H011 * x + 3952. / 27 * H011 * x2
                 + 104. / 3 * H0110 - 64. / 3 * H0110 * x + 32. / 3 * H0110 * x2
                 + 304. / 9 * H0111 - 176. / 9 * H0111 * x
                 + 80. / 9 * H0111 * x2 + 104. / 3 * H0101 - 64. / 3 * H0101 * x
                 + 32. / 3 * H0101 * x2 + 4408. / 9 * H001
                 + 32 * H001 * zeta2 * x - 16 * H001 * zeta2 - 372 * H001 * x
                 + 304. / 3 * H001 * x2 + 224. / 5 * H001 * x3 + 88 * H0010
                 - 160 * H0010 * x - 32. / 3 * H0010 * x2 + 16 * H00100
                 - 32 * H00100 * x + 88 * H0011 - 160 * H0011 * x
                 - 32. / 3 * H0011 * x2 + 16 * H00110 - 32 * H00110 * x
                 + 16 * H00111 - 32 * H00111 * x + 16 * H00101 - 32 * H00101 * x
                 + 436. / 3 * H0001 - 1328. / 3 * H0001 * x
                 - 160. / 3 * H0001 * x2 + 48 * H00010 - 96 * H00010 * x
                 + 48 * H00011 - 96 * H00011 * x + 88 * H00001
                 - 176 * H00001 * x)
        + nf * nf * CA
              * (-78049. / 1215 + 160. / 27 * zeta2 / x
                 + 48460. / 81 * zeta2 * x - 360 * zeta2 * x2
                 - 48. / 5 * zeta2 * x3 - 100. / 9 * zeta2
                 - 688. / 45 * zeta2 * zeta2 * x
                 + 104. / 15 * zeta2 * zeta2 * x2 + 56. / 15 * zeta2 * zeta2
                 + 64. / 27 * zeta3 / x + 2768. / 9 * zeta3 * x
                 - 2992. / 27 * zeta3 * x2 - 268. / 9 * zeta3 - 3832. / 3645 / x
                 - 1441136. / 1215 * x + 4708987. / 3645 * x2 + 32. / 3 * H00m10
                 - 32. / 3 * H00m10 * x2 - 32. / 3 * H0m1 * zeta2 * x
                 + 32. / 3 * H0m1 * zeta2 * x2 - 16. / 3 * H0m1 * zeta2
                 - 32. / 3 * H0m1m10 - 64. / 3 * H0m1m10 * x + 112. / 3 * H0m10
                 + 136. / 9 * H0m10 * x - 272. / 9 * H0m10 * x2
                 + 16. / 3 * H0m100 + 32. / 3 * H0m100 * x
                 - 80. / 3 * H0m100 * x2 - 32. / 3 * H0m101 * x2
                 - 16. / 3 * Hm1 * zeta2 / x - 136. / 9 * Hm1 * zeta2 * x
                 + 256. / 9 * Hm1 * zeta2 * x2 - 248. / 9 * Hm1 * zeta2
                 - 16 * Hm1 * zeta3 * x + 8. / 3 * Hm1 * zeta3 * x2
                 - 8 * Hm1 * zeta3 - 16 * Hm10m10 - 32 * Hm10m10 * x
                 - 32. / 3 * Hm10m10 * x2 + 80. / 3 * Hm1m1 * zeta2 * x
                 + 16. / 3 * Hm1m1 * zeta2 * x2 + 40. / 3 * Hm1m1 * zeta2
                 + 16 * Hm1m1m10 + 32 * Hm1m1m10 * x + 32. / 3 * Hm1m1m10 * x2
                 - 368. / 9 * Hm1m10 - 32. / 9 * Hm1m10 / x
                 - 496. / 9 * Hm1m10 * x - 32. / 9 * Hm1m10 * x2
                 - 64. / 3 * Hm1m100 - 128. / 3 * Hm1m100 * x
                 - 16. / 3 * Hm1m100 * x2 - 16. / 3 * Hm1m101
                 - 32. / 3 * Hm1m101 * x + 5888. / 81 * Hm10
                 + 16 * Hm10 * zeta2 * x + 80. / 3 * Hm10 * zeta2 * x2
                 + 8 * Hm10 * zeta2 - 4. / 15 * Hm10 / x2 + 160. / 27 * Hm10 / x
                 + 2656. / 81 * Hm10 * x - 5036. / 81 * Hm10 * x2
                 - 48. / 5 * Hm10 * x3 + 1192. / 27 * Hm100
                 + 32. / 3 * Hm100 / x + 128. / 27 * Hm100 * x
                 - 1768. / 27 * Hm100 * x2 - 112. / 9 * Hm1000
                 - 224. / 9 * Hm1000 * x - 368. / 9 * Hm1000 * x2
                 + 64. / 9 * Hm101 + 32. / 9 * Hm101 / x - 112. / 9 * Hm101 * x
                 - 272. / 9 * Hm101 * x2 - 16. / 3 * Hm1010
                 - 32. / 3 * Hm1010 * x - 32. / 3 * Hm1010 * x2
                 - 16. / 3 * Hm1011 - 32. / 3 * Hm1011 * x
                 - 32. / 3 * Hm1011 * x2 - 32. / 3 * Hm1001
                 - 64. / 3 * Hm1001 * x - 80. / 3 * Hm1001 * x2
                 - 58412. / 1215 * H0 + 932. / 3 * H0 * zeta2 * x
                 - 952. / 9 * H0 * zeta2 * x2 + 28. / 9 * H0 * zeta2
                 + 704. / 9 * H0 * zeta3 * x - 16 * H0 * zeta3 * x2
                 - 64. / 9 * H0 * zeta3 + 4. / 15 * H0 / x
                 - 1363772. / 1215 * H0 * x + 285488. / 405 * H0 * x2
                 - 1088. / 81 * H00 + 80 * H00 * zeta2 * x
                 - 32. / 3 * H00 * zeta2 * x2 - 8. / 3 * H00 * zeta2
                 - 7244. / 9 * H00 * x + 20348. / 81 * H00 * x2
                 + 48. / 5 * H00 * x3 - 532. / 27 * H000
                 - 10280. / 27 * H000 * x + 1880. / 27 * H000 * x2
                 - 8. / 9 * H0000 - 272. / 3 * H0000 * x + 862. / 243 * H1
                 + 16. / 3 * H1 * zeta2 / x + 1432. / 27 * H1 * zeta2 * x
                 - 904. / 27 * H1 * zeta2 * x2 - 560. / 27 * H1 * zeta2
                 - 16 * H1 * zeta3 * x + 56. / 3 * H1 * zeta3 * x2
                 + 8 * H1 * zeta3 + 80. / 27 * H1 / x - 181160. / 243 * H1 * x
                 + 190376. / 243 * H1 * x2 + 32. / 3 * H10m10
                 - 64. / 3 * H10m10 * x + 32. / 3 * H10m10 * x2
                 + 320. / 27 * H10 + 64. / 9 * H10 * zeta2 * x
                 + 80. / 9 * H10 * zeta2 * x2 - 32. / 9 * H10 * zeta2
                 - 160. / 27 * H10 / x - 7168. / 27 * H10 * x
                 + 7640. / 27 * H10 * x2 + 688. / 27 * H100 - 32. / 9 * H100 / x
                 - 2720. / 27 * H100 * x + 2432. / 27 * H100 * x2
                 + 56. / 3 * H1000 - 112. / 3 * H1000 * x + 64. / 3 * H1000 * x2
                 + 1112. / 81 * H11 + 320. / 9 * H11 * zeta2 * x
                 - 176. / 9 * H11 * zeta2 * x2 - 160. / 9 * H11 * zeta2
                 - 160. / 27 * H11 / x - 21832. / 81 * H11 * x
                 + 23392. / 81 * H11 * x2 + 616. / 27 * H110
                 - 32. / 9 * H110 / x - 3152. / 27 * H110 * x
                 + 3320. / 27 * H110 * x2 + 184. / 9 * H1100
                 - 368. / 9 * H1100 * x + 320. / 9 * H1100 * x2
                 + 232. / 27 * H111 - 32. / 9 * H111 / x - 1808. / 27 * H111 * x
                 + 1976. / 27 * H111 * x2 + 56. / 9 * H1110
                 - 112. / 9 * H1110 * x + 160. / 9 * H1110 * x2
                 + 56. / 9 * H1111 - 112. / 9 * H1111 * x
                 + 112. / 9 * H1111 * x2 + 88. / 9 * H1101
                 - 176. / 9 * H1101 * x + 128. / 9 * H1101 * x2 + 8. / 27 * H101
                 - 32. / 9 * H101 / x - 688. / 27 * H101 * x
                 + 856. / 27 * H101 * x2 + 40. / 9 * H1010 - 80. / 9 * H1010 * x
                 + 128. / 9 * H1010 * x2 - 8. / 9 * H1011 + 16. / 9 * H1011 * x
                 - 16. / 9 * H1011 * x2 + 8. / 9 * H1001 - 16. / 9 * H1001 * x
                 - 80. / 9 * H1001 * x2 + 100. / 9 * H01 + 32 * H01 * zeta2 * x
                 - 15268. / 27 * H01 * x + 360 * H01 * x2 + 8. / 3 * H010
                 - 1556. / 9 * H010 * x + 280. / 3 * H010 * x2 + 16. / 3 * H0100
                 - 128. / 3 * H0100 * x + 16 * H0100 * x2 + 8. / 3 * H011
                 - 1576. / 9 * H011 * x + 280. / 3 * H011 * x2 + 16. / 3 * H0110
                 - 128. / 3 * H0110 * x + 64. / 3 * H0110 * x2 - 32 * H0111 * x
                 + 32. / 3 * H0111 * x2 - 16. / 3 * H0101 - 64. / 3 * H0101 * x
                 - 28. / 9 * H001 - 2660. / 9 * H001 * x + 952. / 9 * H001 * x2
                 - 160. / 3 * H0010 * x + 32. / 3 * H0010 * x2
                 - 160. / 3 * H0011 * x + 32. / 3 * H0011 * x2 + 8. / 3 * H0001
                 - 80 * H0001 * x + 32. / 3 * H0001 * x2);
}

//==========================================================================================//
//  Massless quark coefficient functions for F2 at
//  O(alpha_s^3) for mu=Q.
//  The term fl_ps_11 is put to zero for the reason
//  explained in page 15 of arXiv:1205.5727
//
//  Eq. (B.10) of Ref. [hep-ph/0504242v1]
//------------------------------------------------------------------------------------------//

double C2_ps3_massless(double x, int nf) {

    // the commented varibles are needed for the flav term

    double x2 = x * x;
    double x3 = x2 * x;

    // double xm = 1 - x;

    // double xp = 1 + x;

    // Allocate pointers for the harmonic polylogs
    double wx = x;
    int nw = 5;
    int n1 = -1;
    int n2 = 1;
    int sz = n2 - n1 + 1;
    double *Hr1 = new double[sz];
    double *Hr2 = new double[sz * sz];
    double *Hr3 = new double[sz * sz * sz];
    double *Hr4 = new double[sz * sz * sz * sz];
    double *Hr5 = new double[sz * sz * sz * sz * sz];

    // Call polylogs
    apf_hplog_(&wx, &nw, Hr1, Hr2, Hr3, Hr4, Hr5, &n1, &n2);

    // weight 1
    const double Hm1 = Hr1[0];
    const double H0 = Hr1[1];
    const double H1 = Hr1[2];

    // weight 2
    const double Hm1m1 = Hr2[0];
    const double H0m1 = Hr2[1];
    const double Hm10 = Hr2[3];
    const double H00 = Hr2[4];
    const double H10 = Hr2[5];
    const double H01 = Hr2[7];
    const double H11 = Hr2[8];

    // weight 3
    const double H0m1m1 = Hr3[1];
    const double H00m1 = Hr3[4];
    const double Hm1m10 = Hr3[9];
    const double H0m10 = Hr3[10];
    const double Hm100 = Hr3[12];
    const double H000 = Hr3[13];
    const double H100 = Hr3[14];
    const double H010 = Hr3[16];
    const double H110 = Hr3[17];
    const double Hm101 = Hr3[21];
    const double H001 = Hr3[22];
    const double H101 = Hr3[23];
    const double H011 = Hr3[25];
    const double H111 = Hr3[26];

    // weight 4
    const double Hm1m1m10 = Hr4[27];
    const double H0m1m10 = Hr4[28];
    const double Hm10m10 = Hr4[30];
    const double H00m10 = Hr4[31];
    const double H10m10 = Hr4[32];
    const double Hm1m100 = Hr4[36];
    const double H0m100 = Hr4[37];
    const double Hm1000 = Hr4[39];
    const double H0000 = Hr4[40];
    const double H1000 = Hr4[41];
    const double H0100 = Hr4[43];
    const double H1100 = Hr4[44];
    const double Hm1010 = Hr4[48];
    const double H0010 = Hr4[49];
    const double H1010 = Hr4[50];
    const double H0110 = Hr4[52];
    const double H1110 = Hr4[53];
    const double Hm1m101 = Hr4[63];
    const double H0m101 = Hr4[64];
    const double Hm1001 = Hr4[66];
    const double H0001 = Hr4[67];
    const double H1001 = Hr4[68];
    const double H0101 = Hr4[70];
    const double H1101 = Hr4[71];
    const double Hm1011 = Hr4[75];
    const double H0011 = Hr4[76];
    const double H1011 = Hr4[77];
    const double H0111 = Hr4[79];
    const double H1111 = Hr4[80];

    //  weight 5
    const double H0m1m1m10 = Hr5[82];
    const double H00m1m10 = Hr5[85];
    const double H0m10m10 = Hr5[91];
    const double H000m10 = Hr5[94];
    const double H010m10 = Hr5[97];
    const double H0m1m100 = Hr5[109];
    const double H00m100 = Hr5[112];
    const double H0m1000 = Hr5[118];
    const double H00000 = Hr5[121];
    const double H01000 = Hr5[124];
    // const double Hm10100 = Hr5[129];
    const double H00100 = Hr5[130];
    const double H01100 = Hr5[133];
    const double H0m1010 = Hr5[145];
    const double H00010 = Hr5[148];
    const double H01010 = Hr5[151];
    const double H00110 = Hr5[157];
    const double H01110 = Hr5[160];
    const double H0m1m101 = Hr5[190];
    const double H00m101 = Hr5[193];
    const double H0m1001 = Hr5[199];
    // const double Hm10001 = Hr5[201];
    const double H00001 = Hr5[202];
    const double H01001 = Hr5[205];
    const double H00101 = Hr5[211];
    const double H01101 = Hr5[214];
    const double H0m1011 = Hr5[226];
    const double H00011 = Hr5[229];
    const double H01011 = Hr5[232];
    const double H00111 = Hr5[238];
    const double H01111 = Hr5[241];

    delete[] Hr1;
    delete[] Hr2;
    delete[] Hr3;
    delete[] Hr4;
    delete[] Hr5;

    // double flav = 5. / 18 * fl11ps(nf) * nf;
    // this contribution is neglected
    // If you want to use it please check the 5/48 since
    // I'm not sure about it

    return /*flav
               * (+192 / xm * zeta3 + 512 / xp * zeta2 + 512 * zeta2 * zeta3 * x
                  + 512. / 3 * zeta2 / x + 11648. / 15 * zeta2 * x
                  - 768 * zeta2 * x2 - 4992. / 5 * zeta2 * x3
                  - 7552. / 15 * zeta2 + 2624. / 5 * zeta2 * zeta2 * x
                  + 3072. / 5 * zeta2 * zeta2 * x2
                  - 4608. / 25 * zeta2 * zeta2 * x3 + 1024. / 15 * zeta3 / x
                  - 51712. / 15 * zeta3 * x + 1344. / 5 * zeta3 * x2
                  - 1920 * zeta3 * x3 - 3136. / 15 * zeta3 + 5120 * zeta5 * x
                  + 128. / 5 / x - 1728. / 5 * x + 1152. / 5 * x2 - 192. / 5
                  - 896. / 3 * H0m1 * zeta2 * x + 3072. / 5 * H0m1 * zeta2 * x3
                  + 256 * H0m1m10 * x + 128 * H0m10 / xm - 128. / 3 * H0m10 * x
                  - 1152 * H0m10 * x2 - 768 * H0m10 * x3 - 384 * H0m10
                  + 256. / 3 * H0m100 * x - 1536. / 5 * H0m100 * x3
                  + 1280. / 3 * H0m101 * x - 3072. / 5 * H0m101 * x3
                  - 128 * Hm1 * zeta2 / x2 - 1024. / 15 * Hm1 * zeta2 / x
                  + 8064. / 5 * Hm1 * zeta2 * x + 192. / 5 * Hm1 * zeta2 * x2
                  + 1152 * Hm1 * zeta2 * x3 + 11072. / 15 * Hm1 * zeta2
                  - 4096. / 5 * Hm1 * zeta2 * zeta2 * x
                  - 2048. / 5 * Hm1 * zeta2 * zeta2
                  - 256. / 5 * Hm1 * zeta3 / x2 - 512 * Hm1 * zeta3 * x
                  + 2304. / 5 * Hm1 * zeta3 * x3
                  + 1024. / 15 * Hm1m1 * zeta2 / x2
                  + 2048. / 3 * Hm1m1 * zeta2 * x
                  - 3072. / 5 * Hm1m1 * zeta2 * x3 - 256. / 3 * Hm1m10 / x2
                  + 1280. / 3 * Hm1m10 * x + 1152 * Hm1m10 * x2
                  + 768 * Hm1m10 * x3 + 128 * Hm1m10 - 512. / 15 * Hm1m100 / x2
                  - 1024. / 3 * Hm1m100 * x + 1536. / 5 * Hm1m100 * x3
                  - 1024. / 15 * Hm1m101 / x2 - 2048. / 3 * Hm1m101 * x
                  + 3072. / 5 * Hm1m101 * x3 - 512. / 15 * Hm10 * zeta2 / x2
                  - 1408. / 3 * Hm10 * zeta2 * x + 1536. / 5 * Hm10 * zeta2 * x3
                  - 1024 * Hm10 * zeta3 * x - 512 * Hm10 * zeta3
                  + 1664. / 15 * Hm10 / x2 + 256. / 3 * Hm10 / x
                  - 3200. / 3 * Hm10 * x - 768 * Hm10 * x2
                  - 4992. / 5 * Hm10 * x3 - 1280. / 3 * Hm10
                  + 1024 * Hm100 * zeta2 * x + 512 * Hm100 * zeta2
                  + 256. / 3 * Hm100 / x2 + 512. / 15 * Hm100 / x
                  - 13696. / 15 * Hm100 * x - 1536. / 5 * Hm100 * x2
                  - 768 * Hm100 * x3 - 6016. / 15 * Hm100 + 128 * Hm1000 * x
                  + 256. / 3 * Hm101 / x2 + 1024. / 15 * Hm101 / x
                  - 20992. / 15 * Hm101 * x + 2688. / 5 * Hm101 * x2
                  - 768 * Hm101 * x3 - 10112. / 15 * Hm101 + 1024 * Hm10100 * x
                  + 512 * Hm10100 + 512. / 15 * Hm1001 / x2
                  + 1024. / 3 * Hm1001 * x - 1536. / 5 * Hm1001 * x3
                  - 1024 * Hm10001 * x - 512 * Hm10001 + 64 * H0 / xm * zeta2
                  - 64 * H0 / xp * zeta2 - 320 * H0 / xp
                  - 2944. / 3 * H0 * zeta2 * x + 2688. / 5 * H0 * zeta2 * x2
                  - 1536 * H0 * zeta2 * x3 - 512. / 15 * H0 * zeta2
                  + 2048. / 3 * H0 * zeta3 * x + 768 * H0 * zeta3 * x2
                  - 3072. / 5 * H0 * zeta3 * x3 - 128. / 5 * H0 / x
                  - 384. / 5 * H0 * x + 4992. / 5 * H0 * x2 + 4672. / 15 * H0
                  + 128 * H00 / xm - 256 * H00 / xp - 512 * H00 * zeta2 * x
                  - 768 * H00 * zeta2 * x2 - 256. / 3 * H00 / x
                  - 5888. / 15 * H00 * x + 768 * H00 * x2 + 4992. / 5 * H00 * x3
                  + 512. / 3 * H00 - 64 * H000 / xm + 64 * H000 / xp
                  + 128. / 3 * H000 * x + 768 * H000 * x3
                  - 128. / 3 * H1 * zeta2 / x2 - 640. / 3 * H1 * zeta2 * x
                  + 576 * H1 * zeta2 * x2 - 384 * H1 * zeta2 * x3
                  + 64 * H1 * zeta2 - 256. / 15 * H1 * zeta3 / x2
                  + 2432. / 3 * H1 * zeta3 * x + 768 * H1 * zeta3 * x2
                  - 768. / 5 * H1 * zeta3 * x3 - 1024 * H1 * zeta3
                  + 256. / 3 * H1 / x - 768. / 5 * H1 * x + 768 * H1 * x2
                  - 10496. / 15 * H1 + 256 * Hm10m10 * x
                  + 512. / 15 * H10 * zeta2 / x2 - 1408. / 3 * H10 * zeta2 * x
                  - 768 * H10 * zeta2 * x2 + 1536. / 5 * H10 * zeta2 * x3
                  + 1024 * H10 * zeta2 - 512. / 15 * H100 / x
                  - 128. / 5 * H100 * x - 1536. / 5 * H100 * x2
                  + 3584. / 15 * H100 - 128 * H1000 * x + 512. / 15 * H1100 / x2
                  - 1792. / 3 * H1100 * x - 768 * H1100 * x2
                  + 1536. / 5 * H1100 * x3 + 1024 * H1100
                  - 512. / 15 * H1001 / x2 + 1792. / 3 * H1001 * x
                  + 768 * H1001 * x2 - 1536. / 5 * H1001 * x3 - 1024 * H1001
                  - 512 * H01 / xp - 128 * H01 * zeta2 * x
                  - 1024 * H01 * zeta3 * x - 512 * H01 * zeta3
                  - 256. / 3 * H01 / x - 9216. / 5 * H01 * x + 768 * H01 * x2
                  + 7552. / 15 * H01 + 1024 * H010 * zeta2 * x
                  + 512 * H010 * zeta2 - 1792. / 3 * H0100 * x
                  - 768 * H0100 * x2 + 1536. / 5 * H0100 * x3
                  + 1024 * H01100 * x + 512 * H01100 - 1024 * H01001 * x
                  - 512 * H01001 + 2816. / 3 * H001 * x - 2688. / 5 * H001 * x2
                  + 768 * H001 * x3 + 512. / 15 * H001 + 512 * H0001 * x
                  + 768 * H0001 * x2)*/
        +nf * CF * CA
            * (+488366. / 243 + 904. / 3 * zeta2 * zeta3 * x
               + 296 * zeta2 * zeta3 + 75776. / 405 * zeta2 / x
               - 240892. / 405 * zeta2 * x - 123728. / 405 * zeta2 * x2
               + 16. / 5 * zeta2 * x3 - 65828. / 81 * zeta2
               + 4544. / 45 * zeta2 * zeta2 / x - 908. / 15 * zeta2 * zeta2 * x
               + 1232. / 15 * zeta2 * zeta2 * x2
               - 768. / 25 * zeta2 * zeta2 * x3 + 1072. / 15 * zeta2 * zeta2
               + 3952. / 45 * zeta3 / x + 88132. / 135 * zeta3 * x
               + 57544. / 135 * zeta3 * x2 + 32 * zeta3 * x3
               - 73232. / 135 * zeta3 - 32 * zeta4 / x + 24 * zeta4 * x
               + 32 * zeta4 * x2 - 24 * zeta4 - 132 * zeta5 * x - 620 * zeta5
               - 4852532. / 3645 / x - 130832. / 243 * x - 510478. / 3645 * x2
               - 496. / 3 * H000m10 + 112. / 3 * H000m10 * x
               - 88. / 3 * H00m1 * zeta2 * x + 88. / 3 * H00m1 * zeta2
               - 176. / 3 * H00m1m10 + 368. / 3 * H00m1m10 * x
               + 904. / 3 * H00m10 + 8. / 3 * H00m10 * x
               + 1664. / 9 * H00m10 * x2 - 320. / 3 * H00m100
               + 112. / 3 * H00m100 * x - 176. / 3 * H00m101
               + 272. / 3 * H00m101 * x + 64 * H0m1 * zeta2 / x
               + 860. / 3 * H0m1 * zeta2 * x - 352. / 3 * H0m1 * zeta2 * x2
               - 260 * H0m1 * zeta2 + 128. / 3 * H0m1 * zeta3 * x
               - 128. / 3 * H0m1 * zeta3 - 304. / 3 * H0m10m10
               + 304. / 3 * H0m10m10 * x - 96 * H0m1m1 * zeta2 * x
               + 96 * H0m1m1 * zeta2 + 352. / 3 * H0m1m1m10
               - 352. / 3 * H0m1m1m10 * x - 184 * H0m1m10
               + 128. / 3 * H0m1m10 / x + 88 * H0m1m10 * x
               - 64. / 3 * H0m1m10 * x2 - 344. / 3 * H0m1m100
               + 344. / 3 * H0m1m100 * x - 112. / 3 * H0m1m101
               + 112. / 3 * H0m1m101 * x - 568 * H0m10
               - 208. / 3 * H0m10 * zeta2 * x + 208. / 3 * H0m10 * zeta2
               + 416. / 9 * H0m10 / x + 1256. / 9 * H0m10 * x
               - 5360. / 27 * H0m10 * x2 + 64. / 5 * H0m10 * x3
               + 1136. / 3 * H0m100 - 256. / 3 * H0m100 / x
               - 496. / 3 * H0m100 * x + 1792. / 9 * H0m100 * x2
               - 320. / 3 * H0m1000 + 320. / 3 * H0m1000 * x + 168 * H0m101
               - 128. / 3 * H0m101 / x - 728. / 3 * H0m101 * x
               + 320. / 3 * H0m101 * x2 - 176. / 3 * H0m1010
               + 176. / 3 * H0m1010 * x - 64 * H0m1011 + 64 * H0m1011 * x
               - 344. / 3 * H0m1001 + 344. / 3 * H0m1001 * x
               + 32. / 15 * Hm1 * zeta2 / x2 + 8624. / 27 * Hm1 * zeta2 / x
               - 1348. / 9 * Hm1 * zeta2 * x + 3656. / 27 * Hm1 * zeta2 * x2
               - 96. / 5 * Hm1 * zeta2 * x3 + 116. / 9 * Hm1 * zeta2
               + 16. / 3 * Hm1 * zeta3 / x + 40 * Hm1 * zeta3 * x
               + 16. / 3 * Hm1 * zeta3 * x2 + 40 * Hm1 * zeta3
               - 328. / 3 * Hm10m10 + 256. / 9 * Hm10m10 / x
               - 328. / 3 * Hm10m10 * x + 256. / 9 * Hm10m10 * x2
               - 608. / 9 * Hm1m1 * zeta2 / x - 256. / 3 * Hm1m1 * zeta2 * x
               - 608. / 9 * Hm1m1 * zeta2 * x2 - 256. / 3 * Hm1m1 * zeta2
               + 368. / 3 * Hm1m1m10 - 320. / 9 * Hm1m1m10 / x
               + 368. / 3 * Hm1m1m10 * x - 320. / 9 * Hm1m1m10 * x2
               - 1288. / 9 * Hm1m10 + 64. / 45 * Hm1m10 / x2
               + 3680. / 27 * Hm1m10 / x - 1304. / 9 * Hm1m10 * x
               + 3248. / 27 * Hm1m10 * x2 - 64. / 5 * Hm1m10 * x3 - 76 * Hm1m100
               + 64. / 3 * Hm1m100 / x - 76 * Hm1m100 * x
               + 64. / 3 * Hm1m100 * x2 + 440. / 3 * Hm1m101
               + 448. / 9 * Hm1m101 / x + 440. / 3 * Hm1m101 * x
               + 448. / 9 * Hm1m101 * x2 + 22076. / 135 * Hm10
               - 896. / 9 * Hm10 * zeta2 / x + 56. / 3 * Hm10 * zeta2 * x
               - 896. / 9 * Hm10 * zeta2 * x2 + 56. / 3 * Hm10 * zeta2
               - 16. / 45 * Hm10 / x2 + 44608. / 135 * Hm10 / x
               - 32044. / 135 * Hm10 * x - 9032. / 135 * Hm10 * x2
               + 16. / 5 * Hm10 * x3 - 328 * Hm100 - 64. / 45 * Hm100 / x2
               - 13264. / 27 * Hm100 / x - 776. / 9 * Hm100 * x
               - 6352. / 27 * Hm100 * x2 + 64. / 5 * Hm100 * x3 + 48 * Hm1000
               + 416. / 3 * Hm1000 / x + 48 * Hm1000 * x
               + 416. / 3 * Hm1000 * x2 - 760. / 9 * Hm101
               - 64. / 45 * Hm101 / x2 - 6784. / 27 * Hm101 / x
               + 232. / 3 * Hm101 * x - 2032. / 27 * Hm101 * x2
               + 64. / 5 * Hm101 * x3 + 8 * Hm1010 + 64 * Hm1010 / x
               + 8 * Hm1010 * x + 64 * Hm1010 * x2 + 64. / 3 * Hm1011
               + 704. / 9 * Hm1011 / x + 64. / 3 * Hm1011 * x
               + 704. / 9 * Hm1011 * x2 - 76 * Hm1001 + 320. / 3 * Hm1001 / x
               - 76 * Hm1001 * x + 320. / 3 * Hm1001 * x2 - 1630676. / 1215 * H0
               + 1888. / 45 * H0 * zeta2 / x + 32348. / 135 * H0 * zeta2 * x
               + 91736. / 135 * H0 * zeta2 * x2 + 128. / 5 * H0 * zeta2 * x3
               + 101852. / 135 * H0 * zeta2
               + 1832. / 15 * H0 * zeta2 * zeta2 * x - 56 * H0 * zeta2 * zeta2
               - 128. / 9 * H0 * zeta3 / x - 6920. / 9 * H0 * zeta3 * x
               + 3488. / 9 * H0 * zeta3 * x2 - 192. / 5 * H0 * zeta3 * x3
               - 320. / 9 * H0 * zeta3 - 48 * H0 * zeta4 * x - 48 * H0 * zeta4
               - 198736. / 1215 * H0 / x + 847354. / 1215 * H0 * x
               - 1199092. / 1215 * H0 * x2 + 801616. / 405 * H00
               - 2552. / 9 * H00 * zeta2 * x - 128. / 9 * H00 * zeta2 * x2
               + 192. / 5 * H00 * zeta2 * x3 - 920. / 9 * H00 * zeta2
               - 352 * H00 * zeta3 * x - 1024. / 3 * H00 * zeta3
               + 64. / 45 * H00 / x + 468976. / 405 * H00 * x
               + 22408. / 405 * H00 * x2 - 16. / 5 * H00 * x3
               - 23416. / 27 * H000 - 688. / 3 * H000 * zeta2 * x
               + 320. / 3 * H000 * zeta2 + 15644. / 27 * H000 * x
               - 2792. / 3 * H000 * x2 - 64. / 5 * H000 * x3 + 3416. / 9 * H0000
               + 5360. / 9 * H0000 * x - 240 * H00000 + 240 * H00000 * x
               - 3310. / 3 * H1 + 32. / 45 * H1 * zeta2 / x2
               - 3320. / 27 * H1 * zeta2 / x - 332. / 3 * H1 * zeta2 * x
               + 6128. / 27 * H1 * zeta2 * x2 + 32. / 5 * H1 * zeta2 * x3
               - 4. / 9 * H1 * zeta2 - 64. / 15 * H1 * zeta3 / x2
               - 208. / 3 * H1 * zeta3 / x + 356. / 3 * H1 * zeta3 * x
               + 400. / 3 * H1 * zeta3 * x2 - 192. / 5 * H1 * zeta3 * x3
               - 140 * H1 * zeta3 + 1004296. / 1215 * H1 / x
               + 16250. / 27 * H1 * x - 394996. / 1215 * H1 * x2
               - 16. / 3 * H10m10 + 128. / 9 * H10m10 / x + 16. / 3 * H10m10 * x
               - 128. / 9 * H10m10 * x2 + 6784. / 27 * H10
               + 64. / 15 * H10 * zeta2 / x2 - 256. / 3 * H10 * zeta2 / x
               + 976. / 3 * H10 * zeta2 * x + 64. / 3 * H10 * zeta2 * x2
               + 192. / 5 * H10 * zeta2 * x3 - 304 * H10 * zeta2
               - 27628. / 81 * H10 / x - 2092. / 27 * H10 * x
               + 13552. / 81 * H10 * x2 + 2444. / 45 * H100
               + 1768. / 45 * H100 / x + 5956. / 45 * H100 * x
               - 10168. / 45 * H100 * x2 + 608. / 3 * H1000
               + 896. / 9 * H1000 / x - 608. / 3 * H1000 * x
               - 896. / 9 * H1000 * x2 + 668. / 27 * H11
               - 704. / 9 * H11 * zeta2 / x + 320. / 3 * H11 * zeta2 * x
               + 704. / 9 * H11 * zeta2 * x2 - 320. / 3 * H11 * zeta2
               - 29348. / 81 * H11 / x + 5008. / 27 * H11 * x
               + 12320. / 81 * H11 * x2 + 488. / 9 * H110 + 536. / 27 * H110 / x
               + 232. / 9 * H110 * x - 2696. / 27 * H110 * x2 - 220. / 3 * H1100
               + 64. / 15 * H1100 / x2 + 704. / 9 * H1100 / x
               + 284. / 3 * H1100 * x - 1280. / 9 * H1100 * x2
               + 192. / 5 * H1100 * x3 - 40. / 3 * H111 + 776. / 27 * H111 / x
               + 328. / 3 * H111 * x - 3368. / 27 * H111 * x2 + 136. / 3 * H1110
               + 544. / 9 * H1110 / x - 136. / 3 * H1110 * x
               - 544. / 9 * H1110 * x2 + 112. / 3 * H1111 + 448. / 9 * H1111 / x
               - 112. / 3 * H1111 * x - 448. / 9 * H1111 * x2 + 136. / 3 * H1101
               + 544. / 9 * H1101 / x - 136. / 3 * H1101 * x
               - 544. / 9 * H1101 * x2 - 640. / 9 * H101 + 1480. / 27 * H101 / x
               + 1648. / 9 * H101 * x - 4504. / 27 * H101 * x2 + 52 * H1010
               + 64 * H1010 / x - 52 * H1010 * x - 64 * H1010 * x2 + 40 * H1011
               + 160. / 3 * H1011 / x - 40 * H1011 * x - 160. / 3 * H1011 * x2
               + 740. / 3 * H1001 - 64. / 15 * H1001 / x2 + 704. / 9 * H1001 / x
               - 268 * H1001 * x - 128. / 9 * H1001 * x2 - 192. / 5 * H1001 * x3
               + 65828. / 81 * H01 - 224. / 3 * H01 * zeta2 / x
               - 172. / 3 * H01 * zeta2 * x + 64 * H01 * zeta2 * x2
               - 236. / 3 * H01 * zeta2 - 168 * H01 * zeta3 * x
               - 168 * H01 * zeta3 + 58048. / 405 * H01 / x
               + 28952. / 81 * H01 * x + 123728. / 405 * H01 * x2
               + 32. / 3 * H010m10 + 32. / 3 * H010m10 * x - 2788. / 9 * H010
               - 784. / 3 * H010 * zeta2 * x - 784. / 3 * H010 * zeta2
               - 832. / 9 * H010 / x - 232. / 9 * H010 * x
               - 3880. / 9 * H010 * x2 + 68. / 3 * H0100 + 128. / 3 * H0100 / x
               + 724. / 3 * H0100 * x - 1568. / 9 * H0100 * x2
               + 192. / 5 * H0100 * x3 + 704. / 3 * H01000
               + 704. / 3 * H01000 * x - 11800. / 27 * H011
               - 448. / 3 * H011 * zeta2 * x - 448. / 3 * H011 * zeta2
               - 3280. / 27 * H011 / x - 2260. / 27 * H011 * x
               - 10552. / 27 * H011 * x2 + 40. / 3 * H0110
               + 160. / 3 * H0110 / x + 56. / 3 * H0110 * x
               - 992. / 9 * H0110 * x2 + 56 * H01100 + 56 * H01100 * x
               + 8. / 9 * H0111 + 512. / 9 * H0111 / x + 176. / 9 * H0111 * x
               - 800. / 9 * H0111 * x2 + 272. / 3 * H01110
               + 272. / 3 * H01110 * x + 224. / 3 * H01111
               + 224. / 3 * H01111 * x + 272. / 3 * H01101
               + 272. / 3 * H01101 * x - 40. / 3 * H0101 + 160. / 3 * H0101 / x
               + 40. / 3 * H0101 * x - 224. / 3 * H0101 * x2 + 296. / 3 * H01010
               + 296. / 3 * H01010 * x + 80 * H01011 + 80 * H01011 * x
               + 216 * H01001 + 216 * H01001 * x - 101852. / 135 * H001
               - 712. / 3 * H001 * zeta2 * x - 296. / 3 * H001 * zeta2
               + 64. / 15 * H001 / x - 13508. / 135 * H001 * x
               - 91736. / 135 * H001 * x2 - 64. / 5 * H001 * x3 + 28 * H0010
               + 476. / 3 * H0010 * x - 256. / 3 * H0010 * x2 + 104 * H00100
               + 664. / 3 * H00100 * x - 340. / 9 * H0011
               + 1076. / 9 * H0011 * x - 832. / 9 * H0011 * x2
               + 368. / 3 * H00110 + 208 * H00110 * x + 304. / 3 * H00111
               + 592. / 3 * H00111 * x + 208. / 3 * H00101 + 176 * H00101 * x
               + 920. / 9 * H0001 + 2576. / 9 * H0001 * x
               + 128. / 9 * H0001 * x2 - 192. / 5 * H0001 * x3 + 8 * H00010
               + 232 * H00010 * x + 88. / 3 * H00011 + 760. / 3 * H00011 * x
               - 320. / 3 * H00001 + 800. / 3 * H00001 * x)
        + nf * CF * CF
              * (+186998. / 2025 + 16 * zeta2 * zeta3 * x + 144 * zeta2 * zeta3
                 - 208. / 9 * zeta2 / x - 41204. / 135 * zeta2 * x
                 + 58312. / 405 * zeta2 * x2 + 544. / 75 * zeta2 * x3
                 - 123476. / 135 * zeta2 + 128. / 5 * zeta2 * zeta2 / x
                 + 1384. / 5 * zeta2 * zeta2 * x
                 - 1888. / 9 * zeta2 * zeta2 * x2
                 + 1536. / 25 * zeta2 * zeta2 * x3 + 64 * zeta2 * zeta2
                 - 4384. / 45 * zeta3 / x - 36796. / 15 * zeta3 * x
                 - 26768. / 135 * zeta3 * x2 - 16 * zeta3 * x3
                 - 17752. / 45 * zeta3 + 32 * zeta4 / x - 24 * zeta4 * x
                 - 32 * zeta4 * x2 + 24 * zeta4 + 176 * zeta5 * x + 176 * zeta5
                 + 24146. / 675 / x - 1089248. / 2025 * x + 276604. / 675 * x2
                 - 384 * H000m10 + 448 * H00m1 * zeta2 + 384 * H00m1m10
                 - 256 * H00m1m10 * x - 352 * H00m10 - 160. / 3 * H00m10 * x
                 + 320. / 3 * H00m10 * x2 - 576 * H00m100 + 64 * H00m100 * x
                 - 256 * H00m101 - 128 * H00m101 * x - 176 * H0m1 * zeta2 * x
                 - 224. / 3 * H0m1 * zeta2 * x2 + 624 * H0m1 * zeta2
                 - 224 * H0m1 * zeta3 * x + 224 * H0m1 * zeta3 + 256 * H0m10m10
                 - 256 * H0m10m10 * x + 256 * H0m1m1 * zeta2 * x
                 - 256 * H0m1m1 * zeta2 - 256 * H0m1m1m10 + 256 * H0m1m1m10 * x
                 + 416 * H0m1m10 + 160. / 3 * H0m1m10 * x
                 - 320. / 3 * H0m1m10 * x2 + 448 * H0m1m100 - 448 * H0m1m100 * x
                 + 128 * H0m1m101 - 128 * H0m1m101 * x - 912 * H0m10
                 - 192 * H0m10 * zeta2 * x + 192 * H0m10 * zeta2
                 - 128. / 45 * H0m10 / x2 - 6736. / 9 * H0m10 * x
                 - 992. / 9 * H0m10 * x2 - 32. / 5 * H0m10 * x3 - 736 * H0m100
                 - 32 * H0m100 * x + 320. / 3 * H0m100 * x2 - 256 * H0m1000
                 + 256 * H0m1000 * x - 416 * H0m101 + 608. / 3 * H0m101 * x
                 + 64. / 3 * H0m101 * x2 - 64 * H0m1001 + 64 * H0m1001 * x
                 - 16. / 15 * Hm1 * zeta2 / x2 + 656. / 9 * Hm1 * zeta2 / x
                 + 4744. / 3 * Hm1 * zeta2 * x + 1232. / 9 * Hm1 * zeta2 * x2
                 + 48. / 5 * Hm1 * zeta2 * x3 + 1528 * Hm1 * zeta2
                 - 64 * Hm1 * zeta3 / x + 256 * Hm1 * zeta3 * x
                 - 64 * Hm1 * zeta3 * x2 + 256 * Hm1 * zeta3 + 320 * Hm10m10
                 - 64 * Hm10m10 / x + 320 * Hm10m10 * x - 64 * Hm10m10 * x2
                 + 96 * Hm1m1 * zeta2 / x - 224 * Hm1m1 * zeta2 * x
                 + 96 * Hm1m1 * zeta2 * x2 - 224 * Hm1m1 * zeta2
                 - 320 * Hm1m1m10 + 64 * Hm1m1m10 / x - 320 * Hm1m1m10 * x
                 + 64 * Hm1m1m10 * x2 + 1104 * Hm1m10 - 32. / 45 * Hm1m10 / x2
                 + 992. / 9 * Hm1m10 / x + 9872. / 9 * Hm1m10 * x
                 + 992. / 9 * Hm1m10 * x2 + 32. / 5 * Hm1m10 * x3
                 + 576 * Hm1m100 - 320. / 3 * Hm1m100 / x + 576 * Hm1m100 * x
                 - 320. / 3 * Hm1m100 * x2 + 64 * Hm1m101 - 64 * Hm1m101 / x
                 + 64 * Hm1m101 * x - 64 * Hm1m101 * x2 - 47216. / 45 * Hm10
                 - 32. / 3 * Hm10 * zeta2 / x + 352 * Hm10 * zeta2 * x
                 - 32. / 3 * Hm10 * zeta2 * x2 + 352 * Hm10 * zeta2
                 - 704. / 675 * Hm10 / x2 - 688. / 45 * Hm10 / x
                 - 134368. / 135 * Hm10 * x + 704. / 15 * Hm10 * x2
                 + 544. / 75 * Hm10 * x3 - 5392. / 3 * Hm100
                 - 32. / 15 * Hm100 / x2 - 64 * Hm100 / x
                 - 5552. / 3 * Hm100 * x - 96 * Hm100 * x2
                 + 96. / 5 * Hm100 * x3 - 448 * Hm1000 + 64. / 3 * Hm1000 / x
                 - 448 * Hm1000 * x + 64. / 3 * Hm1000 * x2 - 976 * Hm101
                 + 32. / 45 * Hm101 / x2 - 160. / 9 * Hm101 / x
                 - 9296. / 9 * Hm101 * x - 736. / 9 * Hm101 * x2
                 - 32. / 5 * Hm101 * x3 - 64 * Hm1010 - 64. / 3 * Hm1010 / x
                 - 64 * Hm1010 * x - 64. / 3 * Hm1010 * x2 - 64 * Hm1011
                 - 64. / 3 * Hm1011 / x - 64 * Hm1011 * x
                 - 64. / 3 * Hm1011 * x2 - 192 * Hm1001 - 64. / 3 * Hm1001 / x
                 - 192 * Hm1001 * x - 64. / 3 * Hm1001 * x2
                 + 514634. / 2025 * H0 + 128. / 15 * H0 * zeta2 / x
                 - 74732. / 45 * H0 * zeta2 * x + 13216. / 45 * H0 * zeta2 * x2
                 - 64. / 5 * H0 * zeta2 * x3 - 48028. / 45 * H0 * zeta2
                 + 328. / 3 * H0 * zeta2 * zeta2 * x
                 + 232. / 3 * H0 * zeta2 * zeta2 + 848. / 3 * H0 * zeta3 * x
                 + 1088. / 9 * H0 * zeta3 * x2 + 384. / 5 * H0 * zeta3 * x3
                 - 80. / 3 * H0 * zeta3 + 48 * H0 * zeta4 * x + 48 * H0 * zeta4
                 + 1184. / 675 * H0 / x - 3262136. / 2025 * H0 * x
                 + 303352. / 2025 * H0 * x2 - 1066. / 27 * H00
                 - 2320. / 3 * H00 * zeta2 * x + 5824. / 9 * H00 * zeta2 * x2
                 - 384. / 5 * H00 * zeta2 * x3 - 512. / 3 * H00 * zeta2
                 - 256. / 3 * H00 * zeta3 * x - 640. / 3 * H00 * zeta3
                 + 32. / 15 * H00 / x + 190. / 3 * H00 * x
                 + 872. / 81 * H00 * x2 - 544. / 75 * H00 * x3
                 + 5180. / 9 * H000 - 1328. / 3 * H000 * zeta2 * x
                 - 1328. / 3 * H000 * zeta2 + 3212. / 3 * H000 * x
                 - 1952. / 9 * H000 * x2 - 96. / 5 * H000 * x3
                 - 440. / 3 * H0000 + 880. / 3 * H0000 * x
                 - 3808. / 9 * H0000 * x2 + 240 * H00000 + 240 * H00000 * x
                 + 593936. / 405 * H1 - 16. / 45 * H1 * zeta2 / x2
                 - 320. / 27 * H1 * zeta2 / x + 116. / 3 * H1 * zeta2 * x
                 - 112. / 27 * H1 * zeta2 * x2 - 16. / 5 * H1 * zeta2 * x3
                 - 172. / 9 * H1 * zeta2 + 128. / 15 * H1 * zeta3 / x2
                 + 128. / 9 * H1 * zeta3 / x - 184 * H1 * zeta3 * x
                 - 1280. / 9 * H1 * zeta3 * x2 + 384. / 5 * H1 * zeta3 * x3
                 + 680. / 3 * H1 * zeta3 - 51638. / 405 * H1 / x
                 - 529946. / 405 * H1 * x - 12352. / 405 * H1 * x2
                 + 18868. / 27 * H10 - 128. / 15 * H10 * zeta2 / x2
                 - 160. / 9 * H10 * zeta2 / x - 1408. / 3 * H10 * zeta2 * x
                 + 1312. / 9 * H10 * zeta2 * x2 - 384. / 5 * H10 * zeta2 * x3
                 + 1280. / 3 * H10 * zeta2 + 5704. / 81 * H10 / x
                 - 24436. / 27 * H10 * x + 11000. / 81 * H10 * x2
                 + 16792. / 45 * H100 + 4112. / 135 * H100 / x
                 - 22792. / 45 * H100 * x + 13888. / 135 * H100 * x2
                 - 632. / 3 * H1000 - 224. / 9 * H1000 / x
                 + 632. / 3 * H1000 * x + 224. / 9 * H1000 * x2
                 + 8684. / 9 * H11 - 512. / 9 * H11 * zeta2 / x
                 - 280. / 3 * H11 * zeta2 * x + 512. / 9 * H11 * zeta2 * x2
                 + 280. / 3 * H11 * zeta2 + 5516. / 81 * H11 / x
                 - 9760. / 9 * H11 * x + 4168. / 81 * H11 * x2
                 + 3892. / 9 * H110 - 304. / 27 * H110 / x
                 - 3748. / 9 * H110 * x - 128. / 27 * H110 * x2
                 + 832. / 3 * H1100 - 128. / 15 * H1100 / x2
                 + 832. / 9 * H1100 / x - 320 * H1100 * x
                 + 320. / 9 * H1100 * x2 - 384. / 5 * H1100 * x3
                 + 3904. / 9 * H111 + 64. / 9 * H111 / x - 3904. / 9 * H111 * x
                 - 64. / 9 * H111 * x2 + 200. / 3 * H1110 + 800. / 9 * H1110 / x
                 - 200. / 3 * H1110 * x - 800. / 9 * H1110 * x2
                 + 176. / 3 * H1111 + 704. / 9 * H1111 / x
                 - 176. / 3 * H1111 * x - 704. / 9 * H1111 * x2
                 + 200. / 3 * H1101 + 800. / 9 * H1101 / x
                 - 200. / 3 * H1101 * x - 800. / 9 * H1101 * x2
                 + 5140. / 9 * H101 - 1168. / 27 * H101 / x
                 - 5284. / 9 * H101 * x + 1600. / 27 * H101 * x2
                 + 224. / 3 * H1010 + 896. / 9 * H1010 / x
                 - 224. / 3 * H1010 * x - 896. / 9 * H1010 * x2 + 72 * H1011
                 + 96 * H1011 / x - 72 * H1011 * x - 96 * H1011 * x2
                 - 800. / 3 * H1001 + 128. / 15 * H1001 / x2
                 + 448. / 9 * H1001 / x + 928. / 3 * H1001 * x
                 - 1600. / 9 * H1001 * x2 + 384. / 5 * H1001 * x3
                 + 123476. / 135 * H01 - 72 * H01 * zeta2 * x
                 + 640. / 3 * H01 * zeta2 * x2 + 40. / 3 * H01 * zeta2
                 + 496. / 3 * H01 * zeta3 * x + 496. / 3 * H01 * zeta3
                 + 352. / 45 * H01 / x - 93164. / 135 * H01 * x
                 - 58312. / 405 * H01 * x2 + 6256. / 9 * H010
                 + 640. / 3 * H010 * zeta2 * x + 640. / 3 * H010 * zeta2
                 - 3136. / 9 * H010 * x - 3712. / 27 * H010 * x2
                 + 448. / 3 * H0100 - 304 * H0100 * x - 64 * H0100 * x2
                 - 384. / 5 * H0100 * x3 - 496. / 3 * H01000
                 - 496. / 3 * H01000 * x + 2540. / 3 * H011
                 - 16. / 3 * H011 * zeta2 * x - 16. / 3 * H011 * zeta2
                 - 736. / 3 * H011 * x - 4288. / 27 * H011 * x2
                 + 488. / 3 * H0110 + 40. / 3 * H0110 * x - 224 * H0110 * x2
                 + 704. / 3 * H01100 + 704. / 3 * H01100 * x + 192 * H0111
                 + 160. / 3 * H0111 * x - 1856. / 9 * H0111 * x2
                 + 400. / 3 * H01110 + 400. / 3 * H01110 * x + 352. / 3 * H01111
                 + 352. / 3 * H01111 * x + 400. / 3 * H01101
                 + 400. / 3 * H01101 * x + 584. / 3 * H0101
                 + 136. / 3 * H0101 * x - 800. / 3 * H0101 * x2
                 + 448. / 3 * H01010 + 448. / 3 * H01010 * x + 144 * H01011
                 + 144 * H01011 * x - 256. / 3 * H01001 - 256. / 3 * H01001 * x
                 + 48028. / 45 * H001 - 160 * H001 * zeta2 * x
                 - 96 * H001 * zeta2 - 128. / 15 * H001 / x
                 + 13684. / 15 * H001 * x - 13216. / 45 * H001 * x2
                 + 32. / 5 * H001 * x3 + 584. / 3 * H0010 + 216 * H0010 * x
                 - 3328. / 9 * H0010 * x2 + 544. / 3 * H00100
                 + 544. / 3 * H00100 * x + 880. / 3 * H0011 + 320 * H0011 * x
                 - 3488. / 9 * H0011 * x2 + 736. / 3 * H00110
                 + 736. / 3 * H00110 * x + 704. / 3 * H00111
                 + 704. / 3 * H00111 * x + 288 * H00101 + 288 * H00101 * x
                 + 512. / 3 * H0001 + 720 * H0001 * x - 5824. / 9 * H0001 * x2
                 + 384. / 5 * H0001 * x3 + 1072. / 3 * H00010
                 + 1072. / 3 * H00010 * x + 400 * H00011 + 400 * H00011 * x
                 + 1328. / 3 * H00001 + 1328. / 3 * H00001 * x)
        + nf * nf * CF
              * (-311984. / 1215 - 32. / 9 * zeta2 / x + 9176. / 81 * zeta2 * x
                 - 176. / 9 * zeta2 * x2 + 32. / 5 * zeta2 * x3
                 + 1832. / 81 * zeta2 + 56. / 15 * zeta2 * zeta2 * x
                 + 56. / 15 * zeta2 * zeta2 + 128. / 27 * zeta3 / x
                 + 136. / 27 * zeta3 * x + 160. / 27 * zeta3 * x2
                 + 328. / 27 * zeta3 + 107968. / 3645 / x + 224864. / 1215 * x
                 + 153392. / 3645 * x2 + 128. / 3 * H0m10
                 + 128. / 9 * H0m10 * x2 + 512. / 9 * Hm10
                 - 32. / 45 * Hm10 / x2 - 32. / 9 * Hm10 / x
                 + 640. / 9 * Hm10 * x + 160. / 9 * Hm10 * x2
                 + 32. / 5 * Hm10 * x3 + 128. / 3 * Hm100 + 128. / 9 * Hm100 / x
                 + 128. / 3 * Hm100 * x + 128. / 9 * Hm100 * x2
                 - 219832. / 1215 * H0 + 992. / 27 * H0 * zeta2 * x
                 + 32. / 9 * H0 * zeta2 * x2 + 848. / 27 * H0 * zeta2
                 - 16. / 9 * H0 * zeta3 * x - 16. / 9 * H0 * zeta3
                 + 32. / 45 * H0 / x - 51352. / 1215 * H0 * x
                 + 46688. / 405 * H0 * x2 - 12608. / 81 * H00
                 + 32. / 9 * H00 * zeta2 * x + 32. / 9 * H00 * zeta2
                 - 18080. / 81 * H00 * x + 1744. / 27 * H00 * x2
                 - 32. / 5 * H00 * x3 - 1784. / 27 * H000
                 - 4376. / 27 * H000 * x + 256. / 9 * H000 * x2
                 - 368. / 9 * H0000 - 368. / 9 * H0000 * x - 4280. / 81 * H1
                 - 32. / 9 * H1 * zeta2 / x + 8. / 3 * H1 * zeta2 * x
                 + 32. / 9 * H1 * zeta2 * x2 - 8. / 3 * H1 * zeta2
                 + 112. / 81 * H1 / x + 2120. / 81 * H1 * x
                 + 2048. / 81 * H1 * x2 - 8 * H10 - 160. / 27 * H10 / x
                 + 8. / 3 * H10 * x + 304. / 27 * H10 * x2 + 200. / 9 * H11
                 - 944. / 81 * H11 / x - 200. / 9 * H11 * x
                 + 944. / 81 * H11 * x2 + 8. / 3 * H110 + 32. / 9 * H110 / x
                 - 8. / 3 * H110 * x - 32. / 9 * H110 * x2 + 32. / 9 * H111
                 + 128. / 27 * H111 / x - 32. / 9 * H111 * x
                 - 128. / 27 * H111 * x2 + 8. / 3 * H101 + 32. / 9 * H101 / x
                 - 8. / 3 * H101 * x - 32. / 9 * H101 * x2 - 1832. / 81 * H01
                 - 16. / 3 * H01 * zeta2 * x - 16. / 3 * H01 * zeta2
                 - 3416. / 81 * H01 * x + 176. / 9 * H01 * x2 - 128. / 9 * H010
                 - 176. / 9 * H010 * x - 32. / 9 * H010 * x2 - 152. / 27 * H011
                 - 296. / 27 * H011 * x - 128. / 9 * H011 * x2 + 16. / 3 * H0110
                 + 16. / 3 * H0110 * x + 64. / 9 * H0111 + 64. / 9 * H0111 * x
                 + 16. / 3 * H0101 + 16. / 3 * H0101 * x - 848. / 27 * H001
                 - 992. / 27 * H001 * x - 32. / 9 * H001 * x2 + 112. / 9 * H0011
                 + 112. / 9 * H0011 * x - 32. / 9 * H0001
                 - 32. / 9 * H0001 * x);
}

//==========================================================================================//
//  Massless gluon coefficient functions for FL at
//  O(alpha_s^3) for mu=Q.
//  The term fl_ps_11 is put to zero for the reason
//  explained in page 15 of arXiv:1205.5727
//
//  Eq. (B.17) of Ref. [hep-ph/0504242v1]
//------------------------------------------------------------------------------------------//

double CL_g3_massless(double x, int nf) {

    // the commented varibles are needed for the flav term

    double x2 = x * x;
    double x3 = x2 * x;

    // double xm = 1 - x;
    // double xm2 = xm * xm;
    // double xm3 = xm2 * xm;

    // double xp = 1 + x;
    // double xp2 = xp * xp;
    // double xp3 = xp2 * xp;

    // Allocate pointers for the harmonic polylogs
    double wx = x;
    int nw = 5;
    int n1 = -1;
    int n2 = 1;
    int sz = n2 - n1 + 1;
    double *Hr1 = new double[sz];
    double *Hr2 = new double[sz * sz];
    double *Hr3 = new double[sz * sz * sz];
    double *Hr4 = new double[sz * sz * sz * sz];
    double *Hr5 = new double[sz * sz * sz * sz * sz];

    // Call polylogs
    apf_hplog_(&wx, &nw, Hr1, Hr2, Hr3, Hr4, Hr5, &n1, &n2);

    // weight 1
    const double Hm1 = Hr1[0];
    const double H0 = Hr1[1];
    const double H1 = Hr1[2];

    // weight 2
    const double Hm1m1 = Hr2[0];
    const double H0m1 = Hr2[1];
    const double Hm10 = Hr2[3];
    const double H00 = Hr2[4];
    const double H10 = Hr2[5];
    const double H01 = Hr2[7];
    const double H11 = Hr2[8];

    // weight 3
    const double H0m1m1 = Hr3[1];
    const double Hm1m10 = Hr3[9];
    const double H0m10 = Hr3[10];
    const double Hm100 = Hr3[12];
    const double H000 = Hr3[13];
    const double H100 = Hr3[14];
    const double H010 = Hr3[16];
    const double H110 = Hr3[17];
    const double Hm101 = Hr3[21];
    const double H001 = Hr3[22];
    const double H101 = Hr3[23];
    const double H011 = Hr3[25];
    const double H111 = Hr3[26];

    // weight 4
    const double Hm1m1m10 = Hr4[27];
    const double H0m1m10 = Hr4[28];
    const double Hm10m10 = Hr4[30];
    const double H00m10 = Hr4[31];
    const double H10m10 = Hr4[32];
    const double Hm1m100 = Hr4[36];
    const double H0m100 = Hr4[37];
    const double Hm1000 = Hr4[39];
    const double H0000 = Hr4[40];
    const double H1000 = Hr4[41];
    const double H0100 = Hr4[43];
    const double H1100 = Hr4[44];
    const double Hm1010 = Hr4[48];
    const double H0010 = Hr4[49];
    const double H1010 = Hr4[50];
    const double H0110 = Hr4[52];
    const double H1110 = Hr4[53];
    const double Hm1m101 = Hr4[63];
    const double H0m101 = Hr4[64];
    const double Hm1001 = Hr4[66];
    const double H0001 = Hr4[67];
    const double H1001 = Hr4[68];
    const double H0101 = Hr4[70];
    const double H1101 = Hr4[71];
    const double Hm1011 = Hr4[75];
    const double H0011 = Hr4[76];
    const double H1011 = Hr4[77];
    const double H0111 = Hr4[79];
    const double H1111 = Hr4[80];

    //  weight 5
    const double H0m1m1m10 = Hr5[82];
    const double H00m1m10 = Hr5[85];
    const double H0m10m10 = Hr5[91];
    const double H0m1m100 = Hr5[109];
    const double H00m100 = Hr5[112];
    const double H0m1000 = Hr5[118];
    const double H01000 = Hr5[124];
    const double Hm10100 = Hr5[129];
    const double H00100 = Hr5[130];
    // const double H10100 = Hr5[131];
    const double H01100 = Hr5[133];
    // const double H11100 = Hr5[134];
    const double H0m1m101 = Hr5[190];
    const double H00m101 = Hr5[193];
    const double H0m1001 = Hr5[199];
    const double Hm10001 = Hr5[201];
    // const double H10001 = Hr5[203];
    const double H01001 = Hr5[205];
    // const double H11001 = Hr5[206];

    delete[] Hr1;
    delete[] Hr2;
    delete[] Hr3;
    delete[] Hr4;
    delete[] Hr5;

    // double flav = 5. / 48 * fl11g(nf) * nf * nf;
    // this contribution is neglected
    // If you want to use it please check the 5/48 since
    // I'm not sure about it
    return /*flav
               * (-17152. / 225 - 512. / 75 / xm3 * zeta2 * zeta2
                  + 512. / 75 / xm2 * zeta2 * zeta2 - 256. / 5 / xm2 * zeta3
                  - 128. / 15 / xm * zeta2 + 128. / 5 / xm * zeta3
                  + 896. / 25 / xp3 * zeta2 * zeta2
                  - 896. / 25 / xp2 * zeta2 * zeta2 + 512. / 15 / xp2 * zeta3
                  + 128. / 15 / xp * zeta2 - 256. / 15 / xp * zeta3
                  + 28672. / 225 * zeta2 / x - 6144. / 25 * zeta2 * x
                  - 128768. / 75 * zeta2 * x2 - 768. / 5 * zeta2 * x3
                  + 36352. / 225 * zeta2 + 36992. / 75 * zeta2 * zeta2 * x
                  - 60928. / 75 * zeta2 * zeta2 * x2
                  + 28672. / 25 * zeta2 * zeta2 * x3 + 3584. / 15 * zeta3 / x
                  - 8576. / 45 * zeta3 * x - 57472. / 45 * zeta3 * x2
                  - 12544. / 15 * zeta3 * x3 + 384 * zeta3 + 6784. / 225 / x
                  - 20352. / 25 * x + 21504. / 25 * x2
                  + 512. / 15 * H00m10 / xp3 - 512. / 15 * H00m10 / xp2
                  - 2048. / 15 * H00m10 * x + 16384. / 15 * H00m10 * x2
                  + 1024. / 15 * H0m1 / xm3 * zeta2
                  - 1024. / 15 * H0m1 / xm2 * zeta2
                  - 4096. / 5 * H0m1 * zeta2 * x
                  + 15872. / 15 * H0m1 * zeta2 * x2 + 2048. / 15 * H0m1m10 * x
                  - 1024. / 15 * H0m1m10 * x2 - 256. / 15 * H0m10
                  + 512. / 15 * H0m10 / xp2 - 256. / 15 * H0m10 / xp
                  - 6656. / 45 * H0m10 * x - 18944. / 9 * H0m10 * x2
                  - 25088. / 75 * H0m10 * x3 - 512. / 15 * H0m100 / xm3
                  + 512. / 15 * H0m100 / xm2 + 5632. / 15 * H0m100 * x
                  - 512 * H0m100 * x2 - 1024. / 15 * H0m101 / xm3
                  + 1024. / 15 * H0m101 / xm2 + 13312. / 15 * H0m101 * x
                  - 16384. / 15 * H0m101 * x2 + 1024. / 15 * Hm1 / xm2 * zeta2
                  - 512. / 15 * Hm1 / xm * zeta2 + 6272. / 75 * Hm1 * zeta2 / x2
                  + 4096. / 45 * Hm1 * zeta2 / x + 1408. / 15 * Hm1 * zeta2 * x
                  - 256. / 9 * Hm1 * zeta2 * x2 + 12544. / 25 * Hm1 * zeta2 * x3
                  + 9472. / 15 * Hm1 * zeta2 + 512 * Hm1 * zeta3 * x
                  + 512 * Hm1 * zeta3 * x2 + 512 * Hm10m10 * x
                  + 512 * Hm10m10 * x2 - 512 * Hm1m1 * zeta2 * x
                  - 512 * Hm1m1 * zeta2 * x2 - 1024 * Hm1m1m10 * x
                  - 1024 * Hm1m1m10 * x2 + 6656. / 15 * Hm1m10
                  + 12544. / 225 * Hm1m10 / x2 - 8192. / 45 * Hm1m10 / x
                  + 110336. / 45 * Hm1m10 * x + 18944. / 9 * Hm1m10 * x2
                  + 25088. / 75 * Hm1m10 * x3 + 512 * Hm1m100 * x
                  + 512 * Hm1m100 * x2 - 71936. / 75 * Hm10
                  + 256. / 15 * Hm10 / xp + 512 * Hm10 * zeta2 * x
                  + 512 * Hm10 * zeta2 * x2 - 128. / 5 * Hm10 / x2
                  - 12544. / 225 * Hm10 / x - 316928. / 225 * Hm10 * x
                  - 47488. / 75 * Hm10 * x2 - 768. / 5 * Hm10 * x3
                  - 1280. / 3 * Hm100 - 512. / 15 * Hm100 / xm2
                  + 256. / 15 * Hm100 / xm - 12544. / 225 * Hm100 / x2
                  - 29696. / 45 * Hm100 * x - 512 * Hm100 * x2
                  - 25088. / 75 * Hm100 * x3 - 512 * Hm1000 * x
                  - 512 * Hm1000 * x2 - 2048. / 5 * Hm101
                  - 1024. / 15 * Hm101 / xm2 + 512. / 15 * Hm101 / xm
                  - 12544. / 225 * Hm101 / x2 - 8192. / 45 * Hm101 / x
                  + 50944. / 45 * Hm101 * x + 9728. / 9 * Hm101 * x2
                  - 25088. / 75 * Hm101 * x3 - 1664. / 75 * H0
                  - 256. / 5 * H0 / xm3 * zeta3 - 256. / 5 * H0 / xm2 * zeta2
                  + 256. / 5 * H0 / xm2 * zeta3 + 128. / 5 * H0 / xm * zeta2
                  + 512. / 15 * H0 / xp3 * zeta3 + 256. / 15 * H0 / xp2 * zeta2
                  - 512. / 15 * H0 / xp2 * zeta3 - 128. / 15 * H0 / xp * zeta2
                  - 128. / 15 * H0 / xp - 3584. / 15 * H0 * zeta2 / x
                  - 50176. / 45 * H0 * zeta2 * x
                  + 113152. / 45 * H0 * zeta2 * x2
                  - 50176. / 75 * H0 * zeta2 * x3 - 1792. / 5 * H0 * zeta2
                  + 768 * H0 * zeta3 * x - 27392. / 15 * H0 * zeta3 * x2
                  + 7168. / 5 * H0 * zeta3 * x3 - 6784. / 225 * H0 / x
                  - 147968. / 225 * H0 * x + 25472. / 15 * H0 * x2
                  + 22528. / 225 * H00 - 256. / 5 * H00 / xm3 * zeta2
                  + 256. / 5 * H00 / xm2 * zeta2 + 256. / 15 * H00 / xp3 * zeta2
                  - 256. / 15 * H00 / xp2 * zeta2 - 128. / 15 * H00 / xp
                  - 1792. / 5 * H00 * zeta2 * x + 25856. / 15 * H00 * zeta2 * x2
                  - 7168. / 5 * H00 * zeta2 * x3 + 12544. / 225 * H00 / x
                  + 82816. / 75 * H00 * x + 45568. / 75 * H00 * x2
                  + 768. / 5 * H00 * x3 + 256. / 15 * H000 / xm2
                  - 128. / 15 * H000 / xm - 256. / 15 * H000 / xp2
                  + 128. / 15 * H000 / xp + 6656. / 45 * H000 * x
                  + 25088. / 75 * H000 * x3 + 256. / 15 * H0000 / xm3
                  - 256. / 15 * H0000 / xm2 - 256. / 15 * H0000 / xp3
                  + 256. / 15 * H0000 / xp2 + 2048. / 15 * H0000 * x
                  + 2944. / 25 * H1 + 6272. / 225 * H1 * zeta2 / x2
                  + 4096. / 45 * H1 * zeta2 / x - 55168. / 45 * H1 * zeta2 * x
                  + 9472. / 9 * H1 * zeta2 * x2 - 12544. / 75 * H1 * zeta2 * x3
                  + 3328. / 15 * H1 * zeta2 + 6144. / 5 * H1 * zeta2 * zeta2 * x
                  - 6144. / 5 * H1 * zeta2 * zeta2 * x2
                  - 3584. / 15 * H1 * zeta3 / x2 - 256 * H1 * zeta3 / x
                  + 1792. / 3 * H1 * zeta3 * x - 768 * H1 * zeta3 * x2
                  + 7168. / 5 * H1 * zeta3 * x3 - 768 * H1 * zeta3
                  + 41216. / 225 * H1 / x - 490496. / 225 * H1 * x
                  + 46976. / 25 * H1 * x2 - 512 * H10m10 * x + 512 * H10m10 * x2
                  + 3584. / 15 * H10 * zeta2 / x2 + 256 * H10 * zeta2 / x
                  - 9472. / 3 * H10 * zeta2 * x + 3328 * H10 * zeta2 * x2
                  - 7168. / 5 * H10 * zeta2 * x3 + 768 * H10 * zeta2
                  + 1536 * H10 * zeta3 * x - 1536 * H10 * zeta3 * x2
                  - 1408 * H10 * x + 1408 * H10 * x2 - 5632. / 15 * H100
                  - 1536 * H100 * zeta2 * x + 1536 * H100 * zeta2 * x2
                  - 3584. / 15 * H100 / x - 1536. / 5 * H100 * x
                  + 4608. / 5 * H100 * x2 + 512 * H1000 * x - 512 * H1000 * x2
                  - 512 * H11 * zeta2 * x + 512 * H11 * zeta2 * x2
                  + 3072 * H11 * zeta3 * x - 3072 * H11 * zeta3 * x2
                  - 2816 * H11 * x + 2816 * H11 * x2 - 3072 * H110 * zeta2 * x
                  + 3072 * H110 * zeta2 * x2 + 768 * H1100
                  + 3584. / 15 * H1100 / x2 + 256 * H1100 / x
                  - 6400. / 3 * H1100 * x + 2304 * H1100 * x2
                  - 7168. / 5 * H1100 * x3 - 3072 * H11100 * x
                  + 3072 * H11100 * x2 + 3072 * H11001 * x - 3072 * H11001 * x2
                  - 1536 * H10100 * x + 1536 * H10100 * x2 - 768 * H1001
                  - 3584. / 15 * H1001 / x2 - 256 * H1001 / x
                  + 7936. / 3 * H1001 * x - 2816 * H1001 * x2
                  + 7168. / 5 * H1001 * x3 + 1536 * H10001 * x
                  - 1536 * H10001 * x2 - 36352. / 225 * H01
                  - 1024. / 15 * H01 * zeta2 * x - 512. / 15 * H01 * zeta2 * x2
                  + 1536 * H01 * zeta3 * x - 1536 * H01 * zeta3 * x2
                  - 41216. / 225 * H01 / x - 261632. / 225 * H01 * x
                  + 128768. / 75 * H01 * x2 - 1536 * H010 * zeta2 * x
                  + 1536 * H010 * zeta2 * x2 - 1792. / 3 * H0100 * x
                  + 2304 * H0100 * x2 - 7168. / 5 * H0100 * x3
                  - 1536 * H01100 * x + 1536 * H01100 * x2 + 1536 * H01001 * x
                  - 1536 * H01001 * x2 + 1792. / 5 * H001
                  + 512. / 15 * H001 / xm2 - 256. / 15 * H001 / xm
                  + 3584. / 15 * H001 / x + 8704. / 9 * H001 * x
                  - 113152. / 45 * H001 * x2 + 25088. / 75 * H001 * x3
                  + 512. / 15 * H0001 / xm3 - 512. / 15 * H0001 / xm2
                  + 3328. / 15 * H0001 * x - 25856. / 15 * H0001 * x2
                  + 7168. / 5 * H0001 * x3)*/
        +nf * CF * CA
            * (-137812. / 675 + 256 * zeta2 * zeta3 * x
               + 20896. / 225 * zeta2 / x - 123824. / 675 * zeta2 * x
               - 38208. / 25 * zeta2 * x2 + 117472. / 225 * zeta2 * x3
               - 46304. / 225 * zeta2 - 3712. / 15 * zeta2 * zeta2 * x
               + 14176. / 15 * zeta2 * zeta2 * x2 + 320 * zeta2 * zeta2 * x3
               + 1984. / 15 * zeta3 / x - 121616. / 45 * zeta3 * x
               + 69248. / 45 * zeta3 * x2 + 7648. / 15 * zeta3 * x3
               - 416. / 5 * zeta3 + 3920 * zeta5 * x + 100592. / 675 / x
               - 84264. / 25 * x + 2312348. / 675 * x2 - 256 * H00m1m10 * x
               - 128 * H00m10 - 128 * H00m10 * x + 2432. / 3 * H00m10 * x2
               + 1024. / 5 * H00m10 * x3 + 64 * H00m100 * x - 128 * H00m101 * x
               - 128. / 5 * H0m1 * zeta2 / x2 - 1056 * H0m1 * zeta2 * x
               - 800. / 3 * H0m1 * zeta2 * x2 - 2688. / 5 * H0m1 * zeta2 * x3
               + 192 * H0m1 * zeta2 - 224 * H0m1 * zeta3 * x
               - 256 * H0m10m10 * x + 256 * H0m1m1 * zeta2 * x
               + 256 * H0m1m1m10 * x + 128 * H0m1m10 - 256. / 15 * H0m1m10 / x2
               - 1856 * H0m1m10 * x - 1472. / 3 * H0m1m10 * x2
               - 256. / 5 * H0m1m10 * x3 - 448 * H0m1m100 * x
               - 128 * H0m1m101 * x - 3264. / 5 * H0m10
               - 192 * H0m10 * zeta2 * x + 1792. / 45 * H0m10 / x2
               + 512. / 15 * H0m10 / x - 17728. / 45 * H0m10 * x
               + 51392. / 45 * H0m10 * x2 + 25856. / 75 * H0m10 * x3
               - 256 * H0m100 + 256. / 15 * H0m100 / x2 + 2656. / 3 * H0m100 * x
               + 512 * H0m100 * x2 + 384 * H0m100 * x3 + 256 * H0m1000 * x
               - 128 * H0m101 + 256. / 15 * H0m101 / x2 + 128 * H0m101 * x
               + 64. / 3 * H0m101 * x2 + 512 * H0m101 * x3 + 64 * H0m1001 * x
               - 5584. / 75 * Hm1 * zeta2 / x2 - 2192. / 45 * Hm1 * zeta2 / x
               + 26032. / 15 * Hm1 * zeta2 * x + 10048. / 45 * Hm1 * zeta2 * x2
               - 12928. / 25 * Hm1 * zeta2 * x3 + 5104. / 5 * Hm1 * zeta2
               - 3072. / 5 * Hm1 * zeta2 * zeta2 * x
               - 3072. / 5 * Hm1 * zeta2 * zeta2 * x2
               - 224. / 3 * Hm1 * zeta3 / x2 + 112. / 3 * Hm1 * zeta3 * x
               - 448 * Hm1 * zeta3 * x2 - 448 * Hm1 * zeta3 * x3
               + 112 * Hm1 * zeta3 + 128 * Hm10m10 - 128. / 15 * Hm10m10 / x2
               - 640. / 3 * Hm10m10 * x - 384 * Hm10m10 * x2
               - 256. / 5 * Hm10m10 * x3 + 448. / 5 * Hm1m1 * zeta2 / x2
               - 64 * Hm1m1 * zeta2 * x + 512 * Hm1m1 * zeta2 * x2
               + 2688. / 5 * Hm1m1 * zeta2 * x3 - 128 * Hm1m1 * zeta2
               - 128 * Hm1m1m10 + 128. / 15 * Hm1m1m10 / x2
               + 1024. / 3 * Hm1m1m10 * x + 512 * Hm1m1m10 * x2
               + 256. / 5 * Hm1m1m10 * x3 + 13664. / 15 * Hm1m10
               - 11168. / 225 * Hm1m10 / x2 - 224. / 45 * Hm1m10 / x
               - 6688. / 45 * Hm1m10 * x - 61184. / 45 * Hm1m10 * x2
               - 25856. / 75 * Hm1m10 * x3 + 224 * Hm1m100 - 64 * Hm1m100 / x2
               - 96 * Hm1m100 * x - 640 * Hm1m100 * x2 - 384 * Hm1m100 * x3
               + 64 * Hm1m101 - 256. / 3 * Hm1m101 / x2 + 704. / 3 * Hm1m101 * x
               - 256 * Hm1m101 * x2 - 512 * Hm1m101 * x3 - 388144. / 225 * Hm10
               - 1216. / 15 * Hm10 * zeta2 / x2 + 1216. / 3 * Hm10 * zeta2 * x
               - 96 * Hm10 * zeta2 * x2 - 2432. / 5 * Hm10 * zeta2 * x3
               + 96 * Hm10 * zeta2 - 768 * Hm10 * zeta3 * x
               - 768 * Hm10 * zeta3 * x2 + 62416. / 675 * Hm10 / x2
               + 4928. / 225 * Hm10 / x - 1862192. / 675 * Hm10 * x
               - 43664. / 75 * Hm10 * x2 + 117472. / 225 * Hm10 * x3
               - 992 * Hm100 + 768 * Hm100 * zeta2 * x
               + 768 * Hm100 * zeta2 * x2 + 20128. / 225 * Hm100 / x2
               + 128. / 3 * Hm100 / x - 12512. / 9 * Hm100 * x
               + 416. / 3 * Hm100 * x2 + 14592. / 25 * Hm100 * x3 - 128 * Hm1000
               + 192. / 5 * Hm1000 / x2 - 224 * Hm1000 * x + 96 * Hm1000 * x2
               + 1152. / 5 * Hm1000 * x3 - 1696. / 3 * Hm101
               + 11168. / 225 * Hm101 / x2 + 416. / 9 * Hm101 / x
               - 16288. / 9 * Hm101 * x - 8128. / 9 * Hm101 * x2
               + 25856. / 75 * Hm101 * x3 + 128. / 15 * Hm1010 / x2
               - 128. / 3 * Hm1010 * x + 256. / 5 * Hm1010 * x3
               + 768 * Hm10100 * x + 768 * Hm10100 * x2
               + 128. / 15 * Hm1011 / x2 - 128. / 3 * Hm1011 * x
               + 256. / 5 * Hm1011 * x3 - 32 * Hm1001 + 64 * Hm1001 / x2
               - 224 * Hm1001 * x + 128 * Hm1001 * x2 + 384 * Hm1001 * x3
               - 768 * Hm10001 * x - 768 * Hm10001 * x2 - 6572. / 135 * H0
               + 448. / 15 * H0 * zeta2 / x - 3824. / 15 * H0 * zeta2 * x
               - 33488. / 45 * H0 * zeta2 * x2 + 60512. / 75 * H0 * zeta2 * x3
               - 1232. / 15 * H0 * zeta2 - 192. / 5 * H0 * zeta2 * zeta2 * x
               + 2208 * H0 * zeta3 * x + 2624. / 3 * H0 * zeta3 * x2
               + 5248. / 5 * H0 * zeta3 * x3 - 40432. / 675 * H0 / x
               - 231856. / 75 * H0 * x + 703408. / 225 * H0 * x2
               - 1336. / 225 * H00 + 224 * H00 * zeta2 * x
               - 2240. / 3 * H00 * zeta2 * x2 + 896. / 5 * H00 * zeta2 * x3
               + 128 * H00 * zeta3 * x - 12448. / 225 * H00 / x
               - 123328. / 675 * H00 * x + 41728. / 25 * H00 * x2
               - 117472. / 225 * H00 * x3 + 464. / 5 * H000
               - 256. / 5 * H000 / x - 2352. / 5 * H000 * x
               + 192. / 5 * H000 * x2 - 14592. / 25 * H000 * x3
               - 1312. / 3 * H0000 * x - 1536. / 5 * H0000 * x3
               - 28892. / 225 * H1 - 5584. / 225 * H1 * zeta2 / x2
               + 112. / 45 * H1 * zeta2 / x + 10064. / 45 * H1 * zeta2 * x
               - 43312. / 45 * H1 * zeta2 * x2 + 21728. / 75 * H1 * zeta2 * x3
               + 7072. / 15 * H1 * zeta2 - 1504. / 15 * H1 * zeta3 / x2
               + 128 * H1 * zeta3 / x - 2800. / 3 * H1 * zeta3 * x
               + 544 * H1 * zeta3 * x2 + 3008. / 5 * H1 * zeta3 * x3
               - 240 * H1 * zeta3 + 7936. / 75 * H1 / x - 737108. / 225 * H1 * x
               + 742192. / 225 * H1 * x2 - 128. / 5 * H10m10 / x2
               - 576 * H10m10 * x + 448 * H10m10 * x2 + 768. / 5 * H10m10 * x3
               - 864. / 5 * H10 + 256. / 5 * H10 * zeta2 / x2
               - 128 * H10 * zeta2 / x + 448 * H10 * zeta2 * x
               - 608 * H10 * zeta2 * x2 - 1536. / 5 * H10 * zeta2 * x3
               + 544 * H10 * zeta2 + 608. / 15 * H10 / x - 17968. / 15 * H10 * x
               + 19952. / 15 * H10 * x2 - 3392. / 15 * H100
               - 704. / 15 * H100 / x - 2176. / 5 * H100 * x
               + 10624. / 15 * H100 * x2 - 128 * H1000 + 64. / 5 * H1000 / x2
               + 96 * H1000 * x + 96 * H1000 * x2 - 384. / 5 * H1000 * x3
               - 864. / 5 * H11 - 64. / 15 * H11 * zeta2 / x2
               + 512. / 3 * H11 * zeta2 * x - 256 * H11 * zeta2 * x2
               + 128. / 5 * H11 * zeta2 * x3 + 64 * H11 * zeta2
               + 608. / 15 * H11 / x - 22208. / 15 * H11 * x
               + 8064. / 5 * H11 * x2 - 48 * H110 - 896. / 3 * H110 * x
               + 688. / 3 * H110 * x2 + 352. / 3 * H110 * x3 + 416 * H1100
               + 1024. / 15 * H1100 / x2 - 128 * H1100 / x
               + 1312. / 3 * H1100 * x - 384 * H1100 * x2
               - 2048. / 5 * H1100 * x3 - 16 * H111 - 208 * H111 * x
               + 224 * H111 * x2 - 16 * H101 - 448. / 3 * H101 * x
               + 848. / 3 * H101 * x2 - 352. / 3 * H101 * x3 - 480 * H1001
               - 1024. / 15 * H1001 / x2 + 128 * H1001 / x
               - 1888. / 3 * H1001 * x + 640 * H1001 * x2
               + 2048. / 5 * H1001 * x3 + 46304. / 225 * H01
               - 128. / 15 * H01 * zeta2 / x2 + 960 * H01 * zeta2 * x
               - 736. / 3 * H01 * zeta2 * x2 + 128. / 5 * H01 * zeta2 * x3
               + 64 * H01 * zeta2 - 480 * H01 * zeta3 * x
               - 768 * H01 * zeta3 * x2 - 15968. / 225 * H01 / x
               - 64384. / 25 * H01 * x + 38208. / 25 * H01 * x2
               + 64. / 15 * H010 + 1088 * H010 * zeta2 * x
               + 768 * H010 * zeta2 * x2 - 128. / 15 * H010 / x
               - 4192. / 5 * H010 * x + 1232. / 15 * H010 * x2
               + 352. / 3 * H010 * x3 - 3488. / 3 * H0100 * x - 512 * H0100 * x2
               - 2048. / 5 * H0100 * x3 - 256 * H01000 * x + 64. / 15 * H011
               + 128 * H011 * zeta2 * x - 128. / 15 * H011 / x
               - 12656. / 15 * H011 * x + 1024. / 5 * H011 * x2 - 96 * H0110 * x
               + 832 * H01100 * x + 768 * H01100 * x2 - 32 * H0111 * x
               - 32 * H0101 * x - 960 * H01001 * x - 768 * H01001 * x2
               + 1232. / 15 * H001 + 128 * H001 * zeta2 * x
               + 64. / 15 * H001 / x - 6256. / 45 * H001 * x
               + 33488. / 45 * H001 * x2 - 11552. / 25 * H001 * x3
               - 640. / 3 * H0010 * x - 256. / 5 * H0010 * x3 - 64 * H00100 * x
               - 640. / 3 * H0011 * x - 256. / 5 * H0011 * x3 - 352 * H0001 * x
               + 2240. / 3 * H0001 * x2 + 128. / 5 * H0001 * x3)
        + nf * CF * CF
              * (+1504. / 75 - 128 * zeta2 * zeta3 * x - 2368. / 45 * zeta2 / x
                 + 56872. / 45 * zeta2 * x - 2944. / 15 * zeta2 * x2
                 - 17056. / 75 * zeta2 * x3 + 3856. / 15 * zeta2
                 + 1888. / 5 * zeta2 * zeta2 * x
                 - 2624. / 3 * zeta2 * zeta2 * x2
                 - 8448. / 25 * zeta2 * zeta2 * x3 - 1472. / 15 * zeta3 / x
                 + 127408. / 45 * zeta3 * x - 63008. / 45 * zeta3 * x2
                 - 704. / 3 * zeta3 * x3 + 4448. / 15 * zeta3 - 4000 * zeta5 * x
                 - 1392. / 25 / x + 134096. / 75 * x - 43808. / 25 * x2
                 + 512 * H00m1m10 * x + 256 * H00m10 - 256. / 3 * H00m10 * x
                 - 2816. / 3 * H00m10 * x2 - 512. / 5 * H00m10 * x3
                 - 128 * H00m100 * x + 256 * H00m101 * x
                 + 3904. / 3 * H0m1 * zeta2 * x + 2048. / 3 * H0m1 * zeta2 * x2
                 + 1024. / 5 * H0m1 * zeta2 * x3 - 384 * H0m1 * zeta2
                 + 448 * H0m1 * zeta3 * x + 512 * H0m10m10 * x
                 - 512 * H0m1m1 * zeta2 * x - 512 * H0m1m1m10 * x
                 - 256 * H0m1m10 + 4736. / 3 * H0m1m10 * x
                 + 2048. / 3 * H0m1m10 * x2 + 896 * H0m1m100 * x
                 + 256 * H0m1m101 * x + 12352. / 15 * H0m10
                 + 384 * H0m10 * zeta2 * x - 256. / 15 * H0m10 / x2
                 - 256. / 15 * H0m10 / x + 2944. / 45 * H0m10 * x
                 - 51392. / 45 * H0m10 * x2 - 1408. / 15 * H0m10 * x3
                 + 512 * H0m100 - 3904. / 3 * H0m100 * x - 1024 * H0m100 * x2
                 - 512. / 5 * H0m100 * x3 - 512 * H0m1000 * x + 256 * H0m101
                 - 512 * H0m101 * x - 1024. / 3 * H0m101 * x2
                 - 1024. / 5 * H0m101 * x3 - 128 * H0m1001 * x
                 + 352. / 15 * Hm1 * zeta2 / x2 + 896. / 45 * Hm1 * zeta2 / x
                 - 28256. / 15 * Hm1 * zeta2 * x
                 - 21664. / 45 * Hm1 * zeta2 * x2 + 704. / 5 * Hm1 * zeta2 * x3
                 - 18976. / 15 * Hm1 * zeta2
                 + 3072. / 5 * Hm1 * zeta2 * zeta2 * x
                 + 3072. / 5 * Hm1 * zeta2 * zeta2 * x2
                 + 128. / 5 * Hm1 * zeta3 / x2 + 288 * Hm1 * zeta3 * x
                 + 640 * Hm1 * zeta3 * x2 + 768. / 5 * Hm1 * zeta3 * x3
                 - 224 * Hm1 * zeta3 - 256 * Hm10m10 + 256 * Hm10m10 * x
                 + 512 * Hm10m10 * x2 - 512. / 15 * Hm1m1 * zeta2 / x2
                 - 1024. / 3 * Hm1m1 * zeta2 * x - 768 * Hm1m1 * zeta2 * x2
                 - 1024. / 5 * Hm1m1 * zeta2 * x3 + 256 * Hm1m1 * zeta2
                 + 256 * Hm1m1m10 - 256 * Hm1m1m10 * x - 512 * Hm1m1m10 * x2
                 - 2752. / 3 * Hm1m10 + 704. / 45 * Hm1m10 / x2
                 + 256. / 9 * Hm1m10 / x + 1984. / 9 * Hm1m10 * x
                 + 11200. / 9 * Hm1m10 * x2 + 1408. / 15 * Hm1m10 * x3
                 - 448 * Hm1m100 + 256. / 15 * Hm1m100 / x2
                 + 1472. / 3 * Hm1m100 * x + 1024 * Hm1m100 * x2
                 + 512. / 5 * Hm1m100 * x3 - 128 * Hm1m101
                 + 512. / 15 * Hm1m101 / x2 + 640. / 3 * Hm1m101 * x
                 + 512 * Hm1m101 * x2 + 1024. / 5 * Hm1m101 * x3
                 + 78944. / 45 * Hm10 + 128. / 5 * Hm10 * zeta2 / x2
                 + 64 * Hm10 * zeta2 * x + 384 * Hm10 * zeta2 * x2
                 + 768. / 5 * Hm10 * zeta2 * x3 - 192 * Hm10 * zeta2
                 + 768 * Hm10 * zeta3 * x + 768 * Hm10 * zeta3 * x2
                 - 8368. / 225 * Hm10 / x2 + 64. / 45 * Hm10 / x
                 + 107024. / 45 * Hm10 * x + 2176. / 5 * Hm10 * x2
                 - 17056. / 75 * Hm10 * x3 + 22528. / 15 * Hm100
                 - 768 * Hm100 * zeta2 * x - 768 * Hm100 * zeta2 * x2
                 - 1472. / 45 * Hm100 / x2 - 256. / 15 * Hm100 / x
                 + 66784. / 45 * Hm100 * x - 992. / 5 * Hm100 * x2
                 - 2944. / 15 * Hm100 * x3 + 256 * Hm1000
                 - 128. / 15 * Hm1000 / x2 - 256. / 3 * Hm1000 * x
                 - 384 * Hm1000 * x2 - 256. / 5 * Hm1000 * x3
                 + 4032. / 5 * Hm101 - 704. / 45 * Hm101 / x2
                 - 256. / 45 * Hm101 / x + 89728. / 45 * Hm101 * x
                 + 49664. / 45 * Hm101 * x2 - 1408. / 15 * Hm101 * x3
                 - 768 * Hm10100 * x - 768 * Hm10100 * x2 + 64 * Hm1001
                 - 256. / 15 * Hm1001 / x2 - 320. / 3 * Hm1001 * x
                 - 256 * Hm1001 * x2 - 512. / 5 * Hm1001 * x3
                 + 768 * Hm10001 * x + 768 * Hm10001 * x2 + 2296. / 225 * H0
                 + 64. / 3 * H0 * zeta2 / x + 8224. / 9 * H0 * zeta2 * x
                 + 30656. / 45 * H0 * zeta2 * x2 - 2816. / 15 * H0 * zeta2 * x3
                 + 64. / 15 * H0 * zeta2 + 384. / 5 * H0 * zeta2 * zeta2 * x
                 - 1696 * H0 * zeta3 * x - 2624. / 3 * H0 * zeta3 * x2
                 - 2944. / 5 * H0 * zeta3 * x3 + 2896. / 75 * H0 / x
                 + 294116. / 225 * H0 * x - 32768. / 25 * H0 * x2
                 - 500. / 9 * H00 - 32 * H00 * zeta2 * x
                 + 2240. / 3 * H00 * zeta2 * x2 + 128 * H00 * zeta2 * x3
                 - 256 * H00 * zeta3 * x + 704. / 45 * H00 / x
                 - 8900. / 9 * H00 * x - 48 * H00 * x2 + 17056. / 75 * H00 * x3
                 - 40 * H000 + 256. / 15 * H000 / x - 16864. / 45 * H000 * x
                 + 352 * H000 * x2 + 2944. / 15 * H000 * x3
                 + 496. / 3 * H0000 * x + 512. / 5 * H0000 * x3
                 + 6604. / 45 * H1 + 352. / 45 * H1 * zeta2 / x2
                 - 128. / 9 * H1 * zeta2 / x - 128. / 9 * H1 * zeta2 * x
                 + 3872. / 9 * H1 * zeta2 * x2 - 704. / 15 * H1 * zeta2 * x3
                 - 1088. / 3 * H1 * zeta2 + 1088. / 15 * H1 * zeta3 / x2
                 - 128 * H1 * zeta3 / x + 2144. / 3 * H1 * zeta3 * x
                 - 320 * H1 * zeta3 * x2 - 2176. / 5 * H1 * zeta3 * x3
                 + 96 * H1 * zeta3 - 2432. / 45 * H1 / x + 70036. / 45 * H1 * x
                 - 24736. / 15 * H1 * x2 + 256. / 15 * H10m10 / x2
                 + 1024. / 3 * H10m10 * x - 256 * H10m10 * x2
                 - 512. / 5 * H10m10 * x3 + 104 * H10
                 - 704. / 15 * H10 * zeta2 / x2 + 128 * H10 * zeta2 / x
                 - 1856. / 3 * H10 * zeta2 * x + 960 * H10 * zeta2 * x2
                 + 1408. / 5 * H10 * zeta2 * x3 - 704 * H10 * zeta2
                 - 56 * H10 * x - 48 * H10 * x2 - 1264. / 15 * H100
                 + 832. / 15 * H100 / x + 528. / 5 * H100 * x
                 - 384. / 5 * H100 * x2 + 256 * H1000 - 128. / 15 * H1000 / x2
                 + 256. / 3 * H1000 * x - 384 * H1000 * x2
                 + 256. / 5 * H1000 * x3 + 88 * H11 - 128 * H11 * zeta2 * x
                 + 256 * H11 * zeta2 * x2 - 128 * H11 * zeta2 + 40 * H11 * x
                 - 128 * H11 * x2 - 64 * H110 - 64 * H110 * x + 128 * H110 * x2
                 - 448 * H1100 - 832. / 15 * H1100 / x2 + 128 * H1100 / x
                 - 1600. / 3 * H1100 * x + 576 * H1100 * x2
                 + 1664. / 5 * H1100 * x3 - 80 * H111 - 80 * H111 * x
                 + 160 * H111 * x2 - 96 * H101 - 96 * H101 * x + 192 * H101 * x2
                 + 576 * H1001 + 832. / 15 * H1001 / x2 - 128 * H1001 / x
                 + 1984. / 3 * H1001 * x - 832 * H1001 * x2
                 - 1664. / 5 * H1001 * x3 - 3856. / 15 * H01
                 - 1792. / 3 * H01 * zeta2 * x + 1024. / 3 * H01 * zeta2 * x2
                 - 128 * H01 * zeta2 + 192 * H01 * zeta3 * x
                 + 768 * H01 * zeta3 * x2 + 2432. / 45 * H01 / x
                 + 50152. / 45 * H01 * x + 2944. / 15 * H01 * x2 - 48 * H010
                 - 1408 * H010 * zeta2 * x - 768 * H010 * zeta2 * x2
                 + 64 * H010 * x + 128 * H010 * x2 + 1952. / 3 * H0100 * x
                 + 576 * H0100 * x2 + 1664. / 5 * H0100 * x3 + 512 * H01000 * x
                 - 64 * H011 - 256 * H011 * zeta2 * x + 80 * H011 * x
                 + 160 * H011 * x2 - 128 * H0110 * x - 896 * H01100 * x
                 - 768 * H01100 * x2 - 160 * H0111 * x - 192 * H0101 * x
                 + 1152 * H01001 * x + 768 * H01001 * x2 - 64. / 15 * H001
                 - 256 * H001 * zeta2 * x - 192. / 5 * H001 / x
                 - 38176. / 45 * H001 * x - 30656. / 45 * H001 * x2
                 + 1408. / 15 * H001 * x3 - 224 * H0010 * x + 128 * H00100 * x
                 - 256 * H0011 * x - 160. / 3 * H0001 * x
                 - 2240. / 3 * H0001 * x2 - 1152. / 5 * H0001 * x3)
        + nf * CA * CA
              * (+6472. / 27 - 96 * zeta2 * zeta3 * x + 5248. / 45 * zeta2 / x
                 - 210824. / 45 * zeta2 * x + 319184. / 45 * zeta2 * x2
                 + 152. / 5 * zeta2 * x3 - 56. / 9 * zeta2
                 - 192. / 5 * zeta2 * zeta2 * x
                 - 7672. / 15 * zeta2 * zeta2 * x2
                 - 3456. / 25 * zeta2 * zeta2 * x3 + 64 * zeta3 / x
                 - 1720 * zeta3 * x + 146308. / 45 * zeta3 * x2
                 + 40 * zeta3 * x3 - 4928. / 15 * zeta3 - 960 * zeta5 * x
                 - 222964. / 405 / x + 301316. / 27 * x - 4393856. / 405 * x2
                 + 992. / 3 * H00m10 * x - 608. / 3 * H00m10 * x2
                 - 256. / 5 * H00m10 * x3 - 1648. / 3 * H0m1 * zeta2 * x
                 - 2080. / 3 * H0m1 * zeta2 * x2 + 512. / 5 * H0m1 * zeta2 * x3
                 + 160 * H0m1m10 * x - 1216. / 3 * H0m1m10 * x2
                 - 1264. / 15 * H0m10 + 64. / 5 * H0m10 / x
                 + 5504. / 15 * H0m10 * x + 76024. / 45 * H0m10 * x2
                 + 16 * H0m10 * x3 + 3008. / 3 * H0m100 * x + 448 * H0m100 * x2
                 - 256. / 5 * H0m100 * x3 + 1888. / 3 * H0m101 * x
                 + 1472. / 3 * H0m101 * x2 - 512. / 5 * H0m101 * x3
                 - 4 * Hm1 * zeta2 / x2 + 2768. / 45 * Hm1 * zeta2 / x
                 - 9416. / 5 * Hm1 * zeta2 * x - 71332. / 45 * Hm1 * zeta2 * x2
                 - 24 * Hm1 * zeta2 * x3 - 3848. / 15 * Hm1 * zeta2
                 + 768. / 5 * Hm1 * zeta2 * zeta2 * x
                 + 768. / 5 * Hm1 * zeta2 * zeta2 * x2
                 + 64. / 5 * Hm1 * zeta3 / x2 - 704 * Hm1 * zeta3 * x
                 - 640 * Hm1 * zeta3 * x2 + 384. / 5 * Hm1 * zeta3 * x3
                 - 192 * Hm10m10 * x - 192 * Hm10m10 * x2
                 - 256. / 15 * Hm1m1 * zeta2 / x2
                 + 2368. / 3 * Hm1m1 * zeta2 * x + 704 * Hm1m1 * zeta2 * x2
                 - 512. / 5 * Hm1m1 * zeta2 * x3 + 128 * Hm1m1m10 * x
                 + 128 * Hm1m1m10 * x2 - 272 * Hm1m10 - 8. / 3 * Hm1m10 / x2
                 - 32. / 9 * Hm1m10 / x - 752 * Hm1m10 * x
                 - 4472. / 9 * Hm1m10 * x2 - 16 * Hm1m10 * x3
                 + 128. / 15 * Hm1m100 / x2 - 2432. / 3 * Hm1m100 * x
                 - 768 * Hm1m100 * x2 + 256. / 5 * Hm1m100 * x3
                 + 256. / 15 * Hm1m101 / x2 - 2176. / 3 * Hm1m101 * x
                 - 640 * Hm1m101 * x2 + 512. / 5 * Hm1m101 * x3
                 + 19876. / 45 * Hm10 + 64. / 5 * Hm10 * zeta2 / x2
                 - 816 * Hm10 * zeta2 * x - 752 * Hm10 * zeta2 * x2
                 + 384. / 5 * Hm10 * zeta2 * x3 + 192 * Hm10 * zeta3 * x
                 + 192 * Hm10 * zeta3 * x2 + 76. / 15 * Hm10 / x2
                 + 1688. / 15 * Hm10 / x + 3536. / 45 * Hm10 * x
                 - 10136. / 45 * Hm10 * x2 + 152. / 5 * Hm10 * x3
                 + 64. / 15 * Hm100 - 192 * Hm100 * zeta2 * x
                 - 192 * Hm100 * zeta2 * x2 + 8. / 3 * Hm100 / x2
                 - 1408. / 15 * Hm100 / x + 34864. / 15 * Hm100 * x
                 + 33592. / 15 * Hm100 * x2 + 16 * Hm100 * x3
                 - 64. / 15 * Hm1000 / x2 + 1648. / 3 * Hm1000 * x
                 + 528 * Hm1000 * x2 - 128. / 5 * Hm1000 * x3
                 + 1808. / 15 * Hm101 + 8. / 3 * Hm101 / x2
                 - 2848. / 45 * Hm101 / x + 7536. / 5 * Hm101 * x
                 + 60152. / 45 * Hm101 * x2 + 16 * Hm101 * x3 + 128 * Hm1010 * x
                 + 128 * Hm1010 * x2 - 192 * Hm10100 * x - 192 * Hm10100 * x2
                 + 128 * Hm1011 * x + 128 * Hm1011 * x2
                 - 128. / 15 * Hm1001 / x2 + 2048. / 3 * Hm1001 * x
                 + 640 * Hm1001 * x2 - 256. / 5 * Hm1001 * x3
                 + 192 * Hm10001 * x + 192 * Hm10001 * x2 - 20944. / 45 * H0
                 + 128. / 5 * H0 * zeta2 / x - 49792. / 15 * H0 * zeta2 * x
                 + 99416. / 45 * H0 * zeta2 * x2 + 32 * H0 * zeta2 * x3
                 + 1544. / 15 * H0 * zeta2 - 5888. / 3 * H0 * zeta3 * x
                 - 752. / 3 * H0 * zeta3 * x2 - 256 * H0 * zeta3 * x3
                 - 10052. / 135 * H0 / x + 1012088. / 135 * H0 * x
                 - 24628. / 45 * H0 * x2 + 2764. / 15 * H00
                 - 5392. / 3 * H00 * zeta2 * x + 800. / 3 * H00 * zeta2 * x2
                 + 128. / 5 * H00 * zeta2 * x3 - 56. / 5 * H00 / x
                 + 197672. / 45 * H00 * x - 385064. / 45 * H00 * x2
                 - 152. / 5 * H00 * x3 - 160 * H000 + 128. / 15 * H000 / x
                 + 45296. / 15 * H000 * x - 704 * H000 * x2 - 16 * H000 * x3
                 + 4768. / 3 * H0000 * x + 256. / 5 * H0000 * x3 - 356 * H1
                 - 4. / 3 * H1 * zeta2 / x2 + 688. / 9 * H1 * zeta2 / x
                 - 1512 * H1 * zeta2 * x + 16100. / 9 * H1 * zeta2 * x2
                 + 8 * H1 * zeta2 * x3 - 360 * H1 * zeta2
                 + 448. / 15 * H1 * zeta3 / x2 - 32 * H1 * zeta3 / x
                 - 752. / 3 * H1 * zeta3 * x + 336 * H1 * zeta3 * x2
                 - 896. / 5 * H1 * zeta3 * x3 + 96 * H1 * zeta3
                 - 1432. / 15 * H1 / x + 122360. / 27 * H1 * x
                 - 550852. / 135 * H1 * x2 + 128. / 15 * H10m10 / x2
                 + 224. / 3 * H10m10 * x - 32 * H10m10 * x2
                 - 256. / 5 * H10m10 * x3 + 688. / 3 * H10
                 - 256. / 15 * H10 * zeta2 / x2 + 32 * H10 * zeta2 / x
                 - 1360. / 3 * H10 * zeta2 * x + 432 * H10 * zeta2 * x2
                 + 512. / 5 * H10 * zeta2 * x3 - 96 * H10 * zeta2
                 - 1712. / 9 * H10 / x + 40064. / 9 * H10 * x
                 - 13472. / 3 * H10 * x2 + 1112. / 3 * H100 - 64 * H100 / x
                 + 5008. / 3 * H100 * x - 1976 * H100 * x2
                 - 64. / 15 * H1000 / x2 + 464. / 3 * H1000 * x
                 - 176 * H1000 * x2 + 128. / 5 * H1000 * x3 + 896. / 3 * H11
                 - 384 * H11 * zeta2 * x + 384 * H11 * zeta2 * x2
                 - 1792. / 9 * H11 / x + 43040. / 9 * H11 * x
                 - 43936. / 9 * H11 * x2 + 224 * H110 - 224. / 3 * H110 / x
                 + 6368. / 3 * H110 * x - 2272 * H110 * x2 - 96 * H1100
                 - 64. / 3 * H1100 / x2 + 32 * H1100 / x + 928. / 3 * H1100 * x
                 - 352 * H1100 * x2 + 128 * H1100 * x3 + 192 * H111
                 - 64 * H111 / x + 5536. / 3 * H111 * x - 5920. / 3 * H111 * x2
                 + 448 * H1110 * x - 448 * H1110 * x2 + 384 * H1111 * x
                 - 384 * H1111 * x2 + 448 * H1101 * x - 448 * H1101 * x2
                 + 224 * H101 - 224. / 3 * H101 / x + 1888 * H101 * x
                 - 6112. / 3 * H101 * x2 + 384 * H1010 * x - 384 * H1010 * x2
                 + 448 * H1011 * x - 448 * H1011 * x2 + 96 * H1001
                 + 64. / 3 * H1001 / x2 - 32 * H1001 / x + 1760. / 3 * H1001 * x
                 - 544 * H1001 * x2 - 128 * H1001 * x3 + 56. / 9 * H01
                 - 1360 * H01 * zeta2 * x + 544. / 3 * H01 * zeta2 * x2
                 + 192 * H01 * zeta3 * x + 192 * H01 * zeta3 * x2
                 - 184. / 45 * H01 / x + 42872. / 9 * H01 * x
                 - 319184. / 45 * H01 * x2 + 64 * H010 - 192 * H010 * zeta2 * x
                 - 192 * H010 * zeta2 * x2 - 128. / 3 * H010 / x
                 + 3264 * H010 * x - 2048 * H010 * x2 + 4432. / 3 * H0100 * x
                 - 32 * H0100 * x2 + 128 * H0100 * x3 + 96 * H011
                 - 160. / 3 * H011 / x + 3680 * H011 * x - 6752. / 3 * H011 * x2
                 + 1280 * H0110 * x - 384 * H0110 * x2 - 192 * H01100 * x
                 - 192 * H01100 * x2 + 1152 * H0111 * x - 384 * H0111 * x2
                 + 1280 * H0101 * x - 384 * H0101 * x2 + 192 * H01001 * x
                 + 192 * H01001 * x2 - 1544. / 15 * H001 - 64. / 5 * H001 / x
                 + 18432. / 5 * H001 * x - 99416. / 45 * H001 * x2
                 - 16 * H001 * x3 + 1792 * H0010 * x - 128 * H0010 * x2
                 + 2048 * H0011 * x - 256 * H0011 * x2 + 2128 * H0001 * x
                 - 800. / 3 * H0001 * x2 - 384. / 5 * H0001 * x3)
        + nf * nf * CF
              * (+393496. / 675 - 448. / 45 * zeta2 / x
                 + 27152. / 135 * zeta2 * x - 1568. / 15 * zeta2 * x2
                 - 6848. / 225 * zeta2 * x3 - 6208. / 45 * zeta2
                 + 32. / 5 * zeta2 * zeta2 * x + 3712. / 9 * zeta3 * x
                 - 224. / 3 * zeta3 * x2 + 128. / 3 * zeta3 * x3 - 32 * zeta3
                 + 31888. / 2025 / x - 771296. / 675 * x + 1101512. / 2025 * x2
                 - 128. / 45 * H0m10 / x2 + 1216. / 9 * H0m10 * x
                 - 128. / 3 * H0m10 * x2 - 128. / 15 * H0m10 * x3
                 + 32. / 15 * Hm1 * zeta2 / x2 - 32. / 3 * Hm1 * zeta2 * x
                 + 64. / 5 * Hm1 * zeta2 * x3 + 64. / 45 * Hm1m10 / x2
                 - 64. / 9 * Hm1m10 * x + 128. / 15 * Hm1m10 * x3
                 - 2848. / 45 * Hm10 - 4544. / 675 * Hm10 / x2
                 - 128. / 15 * Hm10 / x + 10976. / 135 * Hm10 * x
                 + 5056. / 45 * Hm10 * x2 - 6848. / 225 * Hm10 * x3
                 - 64. / 15 * Hm100 / x2 + 64. / 3 * Hm100 * x
                 - 128. / 5 * Hm100 * x3 - 64. / 45 * Hm101 / x2
                 + 64. / 9 * Hm101 * x - 128. / 15 * Hm101 * x3
                 + 254648. / 675 * H0 + 3008. / 9 * H0 * zeta2 * x
                 + 416. / 3 * H0 * zeta2 * x2 - 192. / 5 * H0 * zeta2 * x3
                 - 64 * H0 * zeta2 + 128 * H0 * zeta3 * x + 3584. / 675 * H0 / x
                 - 105824. / 225 * H0 * x - 8512. / 225 * H0 * x2
                 + 2288. / 15 * H00 + 192 * H00 * zeta2 * x + 64. / 15 * H00 / x
                 - 57008. / 135 * H00 * x + 1024. / 15 * H00 * x2
                 + 6848. / 225 * H00 * x3 + 80 * H000 - 1648. / 3 * H000 * x
                 - 256 * H000 * x2 + 128. / 5 * H000 * x3 - 352 * H0000 * x
                 + 1464. / 5 * H1 + 32. / 45 * H1 * zeta2 / x2
                 + 128. / 9 * H1 * zeta2 * x + 32. / 3 * H1 * zeta2 * x2
                 - 128. / 5 * H1 * zeta2 * x3 - 352. / 135 * H1 / x
                 - 16696. / 45 * H1 * x + 10912. / 135 * H1 * x2
                 + 224. / 3 * H10 + 64. / 9 * H10 / x - 544. / 3 * H10 * x
                 + 896. / 9 * H10 * x2 + 224. / 3 * H11 + 64. / 9 * H11 / x
                 - 512. / 3 * H11 * x + 800. / 9 * H11 * x2 + 32. / 3 * H110 * x
                 + 32. / 3 * H110 * x2 - 64. / 3 * H110 * x3
                 - 32. / 3 * H101 * x - 32. / 3 * H101 * x2
                 + 64. / 3 * H101 * x3 + 6208. / 45 * H01 + 64. / 45 * H01 / x
                 - 5392. / 45 * H01 * x + 1568. / 15 * H01 * x2 + 32 * H010
                 - 32 * H010 * x2 - 64. / 3 * H010 * x3 + 32 * H011
                 - 32. / 3 * H011 * x - 128. / 3 * H011 * x2 + 64 * H001
                 - 1792. / 9 * H001 * x - 416. / 3 * H001 * x2
                 + 448. / 15 * H001 * x3 - 64 * H0010 * x - 64 * H0011 * x
                 - 192 * H0001 * x)
        + nf * nf * CA
              * (-6364. / 135 - 32. / 9 * zeta2 / x + 2128. / 9 * zeta2 * x
                 - 512. / 3 * zeta2 * x2 - 32. / 5 * zeta2 * x3
                 + 224. / 3 * zeta3 * x - 48 * zeta3 * x2 + 3896. / 135 / x
                 - 97376. / 135 * x + 99844. / 135 * x2 - 32. / 3 * H0m10 * x
                 - 128. / 3 * H0m10 * x2 + 32 * Hm1 * zeta2 * x
                 + 32 * Hm1 * zeta2 * x2 + 64. / 3 * Hm1m10 * x
                 + 64. / 3 * Hm1m10 * x2 - 16. / 3 * Hm10 - 16. / 15 * Hm10 / x2
                 - 32. / 9 * Hm10 / x - 512. / 9 * Hm10 * x
                 - 544. / 9 * Hm10 * x2 - 32. / 5 * Hm10 * x3 - 64 * Hm100 * x
                 - 64 * Hm100 * x2 - 64. / 3 * Hm101 * x - 64. / 3 * Hm101 * x2
                 - 424. / 45 * H0 + 96 * H0 * zeta2 * x
                 - 64. / 3 * H0 * zeta2 * x2 + 16. / 15 * H0 / x
                 - 86912. / 135 * H0 * x + 51584. / 135 * H0 * x2
                 + 16. / 3 * H00 - 3040. / 9 * H00 * x + 416. / 3 * H00 * x2
                 + 32. / 5 * H00 * x3 - 128 * H000 * x - 64. / 3 * H1
                 - 32. / 3 * H1 * zeta2 * x + 32. / 3 * H1 * zeta2 * x2
                 + 400. / 27 * H1 / x - 12032. / 27 * H1 * x
                 + 12208. / 27 * H1 * x2 - 32. / 3 * H10 + 32. / 9 * H10 / x
                 - 1184. / 9 * H10 * x + 416. / 3 * H10 * x2
                 - 64. / 3 * H100 * x + 64. / 3 * H100 * x2 - 32. / 3 * H11
                 + 32. / 9 * H11 / x - 1184. / 9 * H11 * x + 416. / 3 * H11 * x2
                 - 128. / 3 * H110 * x + 128. / 3 * H110 * x2
                 - 64. / 3 * H111 * x + 64. / 3 * H111 * x2 - 880. / 3 * H01 * x
                 + 512. / 3 * H01 * x2 - 64 * H010 * x + 64. / 3 * H010 * x2
                 - 64 * H011 * x + 64. / 3 * H011 * x2 - 320. / 3 * H001 * x
                 + 64. / 3 * H001 * x2);
}

//==========================================================================================//
//  Massless gluon coefficient functions for FL at
//  O(alpha_s^3) for mu=Q.
//  The term fl_ps_11 is put to zero for the reason
//  explained in page 15 of arXiv:1205.5727
//
//  Eq. (B.18) of Ref. [hep-ph/0504242v1]
//------------------------------------------------------------------------------------------//

double CL_ps3_massless(double x, int nf) {

    // the commented varibles are needed for the flav term

    double x2 = x * x;
    double x3 = x2 * x;

    // Allocate pointers for the harmonic polylogs
    double wx = x;
    int nw = 5;
    int n1 = -1;
    int n2 = 1;
    int sz = n2 - n1 + 1;
    double *Hr1 = new double[sz];
    double *Hr2 = new double[sz * sz];
    double *Hr3 = new double[sz * sz * sz];
    double *Hr4 = new double[sz * sz * sz * sz];
    double *Hr5 = new double[sz * sz * sz * sz * sz];

    // Call polylogs
    apf_hplog_(&wx, &nw, Hr1, Hr2, Hr3, Hr4, Hr5, &n1, &n2);

    // weight 1
    const double Hm1 = Hr1[0];
    const double H0 = Hr1[1];
    const double H1 = Hr1[2];

    // weight 2
    // const double Hm1m1 = Hr2[0];
    const double H0m1 = Hr2[1];
    const double Hm10 = Hr2[3];
    const double H00 = Hr2[4];
    const double H10 = Hr2[5];
    const double H01 = Hr2[7];
    const double H11 = Hr2[8];

    // weight 3
    const double Hm1m10 = Hr3[9];
    const double H0m10 = Hr3[10];
    const double Hm100 = Hr3[12];
    const double H000 = Hr3[13];
    const double H100 = Hr3[14];
    const double H010 = Hr3[16];
    const double H110 = Hr3[17];
    const double Hm101 = Hr3[21];
    const double H001 = Hr3[22];
    const double H101 = Hr3[23];
    const double H011 = Hr3[25];
    const double H111 = Hr3[26];

    // weight 4
    const double H0m1m10 = Hr4[28];
    const double H00m10 = Hr4[31];
    // const double H10m10 = Hr4[32];
    // const double Hm1m100 = Hr4[36];
    const double H0m100 = Hr4[37];
    // const double Hm1000 = Hr4[39];
    const double H0000 = Hr4[40];
    // const double H1000 = Hr4[41];
    const double H0100 = Hr4[43];
    const double H1100 = Hr4[44];
    const double H0010 = Hr4[49];
    const double H0110 = Hr4[52];
    // const double Hm1m101 = Hr4[63];
    const double H0m101 = Hr4[64];
    // const double Hm1001 = Hr4[66];
    const double H0001 = Hr4[67];
    const double H1001 = Hr4[68];
    const double H0101 = Hr4[70];
    const double H0011 = Hr4[76];
    const double H0111 = Hr4[79];

    //  weight 5
    // const double Hm10100 = Hr5[129];
    // const double H01100 = Hr5[133];
    // const double Hm10001 = Hr5[201];
    // const double H01001 = Hr5[205];

    delete[] Hr1;
    delete[] Hr2;
    delete[] Hr3;
    delete[] Hr4;
    delete[] Hr5;

    // double flav = 5. / 18 * fl11ps(nf) * nf;
    // this contribution is neglected
    // If you want to use it please check the 5/48 since
    // I'm not sure about it

    return /*flav
               * (+256 * zeta2 * zeta3 * x + 2048. / 3 * zeta2 / x
                  + 10112. / 15 * zeta2 * x - 512 * zeta2 * x2
                  - 3328. / 5 * zeta2 * x3 - 9088. / 15 * zeta2
                  + 576. / 5 * zeta2 * zeta2 * x
                  + 2048. / 5 * zeta2 * zeta2 * x2
                  - 3072. / 25 * zeta2 * zeta2 * x3 + 4096. / 15 * zeta3 / x
                  - 23488. / 15 * zeta3 * x + 896. / 5 * zeta3 * x2
                  - 1280 * zeta3 * x3 - 4864. / 15 * zeta3 + 2560 * zeta5 * x
                  + 512. / 5 / x - 1152. / 5 * x + 768. / 5 * x2 - 768. / 5
                  - 896. / 3 * H0m1 * zeta2 * x + 2048. / 5 * H0m1 * zeta2 * x3
                  + 256 * H0m1m10 * x + 256. / 3 * H0m10 * x - 768 * H0m10 * x2
                  - 512 * H0m10 * x3 + 256. / 3 * H0m100 * x
                  - 1024. / 5 * H0m100 * x3 + 1280. / 3 * H0m101 * x
                  - 2048. / 5 * H0m101 * x3 - 512 * Hm1 * zeta2 / x2
                  - 1216. / 15 * Hm1 * zeta2 / x + 3136. / 5 * Hm1 * zeta2 * x
                  + 128. / 5 * Hm1 * zeta2 * x2 + 768 * Hm1 * zeta2 * x3
                  - 3712. / 15 * Hm1 * zeta2
                  - 2048. / 5 * Hm1 * zeta2 * zeta2 * x
                  - 1024. / 5 * Hm1 * zeta3 / x2 - 512 * Hm1 * zeta3 * x
                  + 1536. / 5 * Hm1 * zeta3 * x3
                  + 4096. / 15 * Hm1m1 * zeta2 / x2
                  + 2048. / 3 * Hm1m1 * zeta2 * x
                  - 2048. / 5 * Hm1m1 * zeta2 * x3 - 1024. / 3 * Hm1m10 / x2
                  - 384 * Hm1m10 / x + 128. / 3 * Hm1m10 * x + 768 * Hm1m10 * x2
                  + 512 * Hm1m10 * x3 - 256 * Hm1m10 - 2048. / 15 * Hm1m100 / x2
                  - 1024. / 3 * Hm1m100 * x + 1024. / 5 * Hm1m100 * x3
                  - 4096. / 15 * Hm1m101 / x2 - 2048. / 3 * Hm1m101 * x
                  + 2048. / 5 * Hm1m101 * x3 - 2048. / 15 * Hm10 * zeta2 / x2
                  - 1408. / 3 * Hm10 * zeta2 * x + 1024. / 5 * Hm10 * zeta2 * x3
                  - 512 * Hm10 * zeta3 * x + 6656. / 15 * Hm10 / x2
                  + 1024. / 3 * Hm10 / x - 1664. / 3 * Hm10 * x
                  - 512 * Hm10 * x2 - 3328. / 5 * Hm10 * x3 + 256. / 3 * Hm10
                  + 512 * Hm100 * zeta2 * x + 1024. / 3 * Hm100 / x2
                  + 2048. / 15 * Hm100 / x - 4864. / 15 * Hm100 * x
                  - 1024. / 5 * Hm100 * x2 - 512 * Hm100 * x3
                  + 2816. / 15 * Hm100 + 128 * Hm1000 * x
                  + 1024. / 3 * Hm101 / x2 - 1664. / 15 * Hm101 / x
                  - 9088. / 15 * Hm101 * x + 1792. / 5 * Hm101 * x2
                  - 512 * Hm101 * x3 + 1792. / 15 * Hm101 + 512 * Hm10100 * x
                  + 2048. / 15 * Hm1001 / x2 + 1024. / 3 * Hm1001 * x
                  - 1024. / 5 * Hm1001 * x3 - 512 * Hm10001 * x
                  - 640. / 3 * H0 * zeta2 * x + 1792. / 5 * H0 * zeta2 * x2
                  - 1024 * H0 * zeta2 * x3 + 1792. / 15 * H0 * zeta2
                  + 512. / 3 * H0 * zeta3 * x + 512 * H0 * zeta3 * x2
                  - 2048. / 5 * H0 * zeta3 * x3 - 512. / 5 * H0 / x
                  - 896. / 5 * H0 * x + 3328. / 5 * H0 * x2 - 512. / 15 * H0
                  - 512 * H00 * zeta2 * x2 - 1024. / 3 * H00 / x
                  - 6272. / 15 * H00 * x + 512 * H00 * x2 + 3328. / 5 * H00 * x3
                  + 512. / 3 * H00 - 256. / 3 * H000 * x + 512 * H000 * x3
                  - 512. / 3 * H1 * zeta2 / x2 + 192 * H1 * zeta2 / x
                  - 64. / 3 * H1 * zeta2 * x + 384 * H1 * zeta2 * x2
                  - 256 * H1 * zeta2 * x3 - 128 * H1 * zeta2
                  - 1024. / 15 * H1 * zeta3 / x2 + 256 * H1 * zeta3 / x
                  + 896. / 3 * H1 * zeta3 * x + 512 * H1 * zeta3 * x2
                  - 512. / 5 * H1 * zeta3 * x3 - 512 * H1 * zeta3
                  + 1024. / 3 * H1 / x - 1152. / 5 * H1 * x + 512 * H1 * x2
                  - 9344. / 15 * H1 + 256 * H10m10 * x
                  + 2048. / 15 * H10 * zeta2 / x2 - 256 * H10 * zeta2 / x
                  + 128. / 3 * H10 * zeta2 * x - 512 * H10 * zeta2 * x2
                  + 1024. / 5 * H10 * zeta2 * x3 + 512 * H10 * zeta2
                  - 2048. / 15 * H100 / x + 128. / 5 * H100 * x
                  - 1024. / 5 * H100 * x2 + 2816. / 15 * H100 - 128 * H1000 * x
                  + 2048. / 15 * H1100 / x2 - 256 * H1100 / x
                  - 256. / 3 * H1100 * x - 512 * H1100 * x2
                  + 1024. / 5 * H1100 * x3 + 512 * H1100
                  - 2048. / 15 * H1001 / x2 + 256 * H1001 / x
                  + 256. / 3 * H1001 * x + 512 * H1001 * x2
                  - 1024. / 5 * H1001 * x3 - 512 * H1001 - 128 * H01 * zeta2 * x
                  - 512 * H01 * zeta3 * x - 1024. / 3 * H01 / x
                  - 6144. / 5 * H01 * x + 512 * H01 * x2 + 9088. / 15 * H01
                  + 512 * H010 * zeta2 * x - 256. / 3 * H0100 * x
                  - 512 * H0100 * x2 + 1024. / 5 * H0100 * x3 + 512 * H01100 * x
                  - 512 * H01001 * x + 896. / 3 * H001 * x
                  - 1792. / 5 * H001 * x2 + 512 * H001 * x3 - 1792. / 15 * H001
                  + 512 * H0001 * x2)*/
        +nf * CF * CA
            * (+6272. / 9 + 5296. / 45 * zeta2 / x - 1024. / 15 * zeta2 * x
               + 3408. / 5 * zeta2 * x2 + 32. / 15 * zeta2 * x3
               - 992. / 9 * zeta2 - 1384. / 15 * zeta2 * zeta2 * x
               + 512. / 15 * zeta2 * zeta2 * x2 - 512. / 25 * zeta2 * zeta2 * x3
               + 1856. / 15 * zeta3 / x - 4384. / 45 * zeta3 * x
               + 2208. / 5 * zeta3 * x2 + 64. / 3 * zeta3 * x3
               - 1824. / 5 * zeta3 - 71344. / 135 / x + 30944. / 27 * x
               - 59152. / 45 * x2 + 224. / 3 * H00m10 * x
               - 64 * H0m1 * zeta2 * x + 256. / 3 * H0m1m10 * x - 96 * H0m10
               + 64. / 3 * H0m10 / x - 1664. / 9 * H0m10 * x
               + 704. / 3 * H0m10 * x2 + 128. / 15 * H0m10 * x3
               + 608. / 3 * H0m100 * x + 320. / 3 * H0m101 * x
               + 128. / 15 * Hm1 * zeta2 / x2 + 80 * Hm1 * zeta2 / x
               - 464. / 3 * Hm1 * zeta2 * x - 160 * Hm1 * zeta2 * x2
               - 64. / 5 * Hm1 * zeta2 * x3 + 64 * Hm1 * zeta2
               - 256. / 3 * Hm1m10 + 256. / 45 * Hm1m10 / x2
               + 32. / 3 * Hm1m10 / x - 928. / 9 * Hm1m10 * x
               - 64. / 3 * Hm1m10 * x2 - 128. / 15 * Hm1m10 * x3
               + 9728. / 45 * Hm10 - 64. / 45 * Hm10 / x2
               + 4144. / 45 * Hm10 / x - 10192. / 45 * Hm10 * x
               - 15616. / 45 * Hm10 * x2 + 32. / 15 * Hm10 * x3
               - 608. / 3 * Hm100 - 256. / 45 * Hm100 / x2
               - 320. / 3 * Hm100 / x + 928. / 9 * Hm100 * x
               + 640. / 3 * Hm100 * x2 + 128. / 15 * Hm100 * x3
               - 320. / 3 * Hm101 - 256. / 45 * Hm101 / x2
               - 224. / 3 * Hm101 / x + 928. / 9 * Hm101 * x
               + 448. / 3 * Hm101 * x2 + 128. / 15 * Hm101 * x3
               - 6832. / 45 * H0 + 64. / 15 * H0 * zeta2 / x
               - 7376. / 45 * H0 * zeta2 * x + 512. / 5 * H0 * zeta2 * x2
               + 256. / 15 * H0 * zeta2 * x3 + 624. / 5 * H0 * zeta2
               - 1568. / 3 * H0 * zeta3 * x + 128. / 3 * H0 * zeta3 * x2
               - 128. / 5 * H0 * zeta3 * x3 - 11456. / 135 * H0 / x
               + 101824. / 135 * H0 * x + 161648. / 135 * H0 * x2
               + 15376. / 45 * H00 - 976. / 3 * H00 * zeta2 * x
               - 128. / 3 * H00 * zeta2 * x2 + 128. / 5 * H00 * zeta2 * x3
               + 256. / 45 * H00 / x + 8192. / 15 * H00 * x
               - 61472. / 45 * H00 * x2 - 32. / 15 * H00 * x3 - 160 * H000
               + 3536. / 9 * H000 * x - 128. / 15 * H000 * x3
               + 1216. / 3 * H0000 * x - 2528. / 9 * H1
               + 128. / 45 * H1 * zeta2 / x2 + 48 * H1 * zeta2 / x
               + 464. / 9 * H1 * zeta2 * x + 96 * H1 * zeta2 * x2
               + 64. / 15 * H1 * zeta2 * x3 - 608. / 3 * H1 * zeta2
               - 256. / 15 * H1 * zeta3 / x2 + 64. / 3 * H1 * zeta3 / x
               - 64. / 3 * H1 * zeta3 * x + 128. / 3 * H1 * zeta3 * x2
               - 128. / 5 * H1 * zeta3 * x3 + 8176. / 135 * H1 / x
               - 3424. / 9 * H1 * x + 81104. / 135 * H1 * x2 + 256. / 3 * H10
               + 256. / 15 * H10 * zeta2 / x2 - 64. / 3 * H10 * zeta2 / x
               + 64. / 3 * H10 * zeta2 * x - 128. / 3 * H10 * zeta2 * x2
               + 128. / 5 * H10 * zeta2 * x3 - 1264. / 9 * H10 / x
               + 752. / 3 * H10 * x - 1760. / 9 * H10 * x2 + 1424. / 5 * H100
               - 1216. / 15 * H100 / x - 752. / 15 * H100 * x
               - 768. / 5 * H100 * x2 + 224. / 3 * H11 - 1168. / 9 * H11 / x
               + 704. / 3 * H11 * x - 1616. / 9 * H11 * x2 + 160 * H110
               - 160. / 3 * H110 / x - 320. / 3 * H110 * x2
               + 256. / 15 * H1100 / x2 - 64. / 3 * H1100 / x
               + 64. / 3 * H1100 * x - 128. / 3 * H1100 * x2
               + 128. / 5 * H1100 * x3 + 160 * H111 - 160. / 3 * H111 / x
               - 320. / 3 * H111 * x2 + 160 * H101 - 160. / 3 * H101 / x
               - 320. / 3 * H101 * x2 - 256. / 15 * H1001 / x2
               + 64. / 3 * H1001 / x - 64. / 3 * H1001 * x
               + 128. / 3 * H1001 * x2 - 128. / 5 * H1001 * x3 + 992. / 9 * H01
               - 608. / 3 * H01 * zeta2 * x - 128. / 5 * H01 / x
               - 1424. / 9 * H01 * x - 3408. / 5 * H01 * x2
               - 128. / 3 * H010 / x + 32 * H010 * x - 448. / 3 * H010 * x2
               + 880. / 3 * H0100 * x - 128. / 3 * H0100 * x2
               + 128. / 5 * H0100 * x3 + 32 * H011 - 160. / 3 * H011 / x
               + 16. / 3 * H011 * x - 448. / 3 * H011 * x2 + 160 * H0110 * x
               + 160 * H0111 * x + 160 * H0101 * x - 624. / 5 * H001
               + 256. / 15 * H001 / x - 944. / 45 * H001 * x
               - 512. / 5 * H001 * x2 - 128. / 15 * H001 * x3 + 320 * H0010 * x
               + 352 * H0011 * x + 400 * H0001 * x + 128. / 3 * H0001 * x2
               - 128. / 5 * H0001 * x3)
        + nf * CF * CF
              * (+25648. / 225 - 352. / 9 * zeta2 / x - 80176. / 135 * zeta2 * x
                 + 1984. / 15 * zeta2 * x2 + 1088. / 225 * zeta2 * x3
                 - 11488. / 45 * zeta2 + 1424. / 15 * zeta2 * zeta2 * x
                 - 1024. / 15 * zeta2 * zeta2 * x2
                 + 1024. / 25 * zeta2 * zeta2 * x3 - 384. / 5 * zeta3 / x
                 - 23152. / 45 * zeta3 * x - 128. / 15 * zeta3 * x2
                 - 32. / 3 * zeta3 * x3 + 128. / 5 * zeta3 + 30464. / 675 / x
                 - 108544. / 675 * x + 1136. / 675 * x2 - 64. / 3 * H00m10 * x
                 - 96 * H0m1 * zeta2 * x - 704. / 3 * H0m1m10 * x
                 - 640. / 3 * H0m10 - 512. / 45 * H0m10 / x2
                 + 320. / 9 * H0m10 * x + 128. / 3 * H0m10 * x2
                 - 64. / 15 * H0m10 * x3 + 64 * H0m100 * x
                 - 64. / 3 * H0m101 * x - 64. / 15 * Hm1 * zeta2 / x2
                 - 32. / 3 * Hm1 * zeta2 / x + 352. / 3 * Hm1 * zeta2 * x
                 + 64. / 3 * Hm1 * zeta2 * x2 + 32. / 5 * Hm1 * zeta2 * x3
                 + 96 * Hm1 * zeta2 + 704. / 3 * Hm1m10
                 - 128. / 45 * Hm1m10 / x2 + 64. / 3 * Hm1m10 / x
                 + 1472. / 9 * Hm1m10 * x - 128. / 3 * Hm1m10 * x2
                 + 64. / 15 * Hm1m10 * x3 - 2176. / 5 * Hm10
                 - 4736. / 675 * Hm10 / x2 - 352. / 45 * Hm10 / x
                 - 56992. / 135 * Hm10 * x + 256. / 15 * Hm10 * x2
                 + 1088. / 225 * Hm10 * x3 - 64 * Hm100 - 128. / 15 * Hm100 / x2
                 - 256. / 3 * Hm100 * x + 64. / 5 * Hm100 * x3 + 64. / 3 * Hm101
                 + 128. / 45 * Hm101 / x2 + 64. / 3 * Hm101 / x
                 - 320. / 9 * Hm101 * x - 128. / 3 * Hm101 * x2
                 - 64. / 15 * Hm101 * x3 + 67952. / 675 * H0
                 + 512. / 15 * H0 * zeta2 / x - 5696. / 15 * H0 * zeta2 * x
                 + 1536. / 5 * H0 * zeta2 * x2 - 128. / 15 * H0 * zeta2 * x3
                 - 1408. / 5 * H0 * zeta2 + 736. / 3 * H0 * zeta3 * x
                 - 256. / 3 * H0 * zeta3 * x2 + 256. / 5 * H0 * zeta3 * x3
                 + 6656. / 675 * H0 / x - 11904. / 25 * H0 * x
                 - 12944. / 675 * H0 * x2 - 112. / 3 * H00
                 - 1024. / 3 * H00 * zeta2 * x + 256. / 3 * H00 * zeta2 * x2
                 - 256. / 5 * H00 * zeta2 * x3 + 128. / 15 * H00 / x
                 + 4544. / 27 * H00 * x - 64 * H00 * x2 - 1088. / 225 * H00 * x3
                 + 160 * H000 + 2000. / 9 * H000 * x - 256 * H000 * x2
                 - 64. / 5 * H000 * x3 + 64 * H0000 * x + 3344. / 5 * H1
                 - 64. / 45 * H1 * zeta2 / x2 + 32. / 3 * H1 * zeta2 / x
                 - 736. / 9 * H1 * zeta2 * x + 64. / 3 * H1 * zeta2 * x2
                 - 32. / 15 * H1 * zeta2 * x3 + 160. / 3 * H1 * zeta2
                 + 512. / 15 * H1 * zeta3 / x2 - 128. / 3 * H1 * zeta3 / x
                 + 128. / 3 * H1 * zeta3 * x - 256. / 3 * H1 * zeta3 * x2
                 + 256. / 5 * H1 * zeta3 * x3 - 19904. / 135 * H1 / x
                 - 13856. / 45 * H1 * x - 28816. / 135 * H1 * x2
                 + 544. / 3 * H10 - 512. / 15 * H10 * zeta2 / x2
                 + 128. / 3 * H10 * zeta2 / x - 128. / 3 * H10 * zeta2 * x
                 + 256. / 3 * H10 * zeta2 * x2 - 256. / 5 * H10 * zeta2 * x3
                 - 64. / 9 * H10 / x - 544. / 3 * H10 * x + 64. / 9 * H10 * x2
                 - 768. / 5 * H100 + 512. / 15 * H100 / x
                 + 1024. / 15 * H100 * x + 256. / 5 * H100 * x2 + 752. / 3 * H11
                 - 208. / 9 * H11 / x - 512. / 3 * H11 * x - 512. / 9 * H11 * x2
                 + 64 * H110 - 64. / 3 * H110 / x - 128. / 3 * H110 * x2
                 - 512. / 15 * H1100 / x2 + 128. / 3 * H1100 / x
                 - 128. / 3 * H1100 * x + 256. / 3 * H1100 * x2
                 - 256. / 5 * H1100 * x3 + 32 * H111 - 32. / 3 * H111 / x
                 - 64. / 3 * H111 * x2 + 64 * H101 - 64. / 3 * H101 / x
                 - 128. / 3 * H101 * x2 + 512. / 15 * H1001 / x2
                 - 128. / 3 * H1001 / x + 128. / 3 * H1001 * x
                 - 256. / 3 * H1001 * x2 + 256. / 5 * H1001 * x3
                 + 11488. / 45 * H01 + 160. / 3 * H01 * zeta2 * x
                 + 1408. / 45 * H01 / x + 2576. / 15 * H01 * x
                 - 1984. / 15 * H01 * x2 + 128 * H010 + 96 * H010 * x
                 - 128 * H010 * x2 - 512. / 3 * H0100 * x
                 + 256. / 3 * H0100 * x2 - 256. / 5 * H0100 * x3 + 128 * H011
                 + 208 * H011 * x - 320. / 3 * H011 * x2 + 64 * H0110 * x
                 + 32 * H0111 * x + 64 * H0101 * x + 1408. / 5 * H001
                 - 512. / 15 * H001 / x + 18688. / 45 * H001 * x
                 - 1536. / 5 * H001 * x2 + 64. / 15 * H001 * x3
                 + 128 * H0010 * x + 128 * H0011 * x + 320 * H0001 * x
                 - 256. / 3 * H0001 * x2 + 256. / 5 * H0001 * x3)
        + nf * nf * CF
              * (-5248. / 45 - 64. / 9 * zeta2 / x + 160. / 3 * zeta2 * x
                 + 64. / 9 * zeta2 * x2 + 64. / 15 * zeta2 * x3
                 - 32. / 3 * zeta3 * x + 15808. / 405 / x + 3584. / 135 * x
                 + 20672. / 405 * x2 - 128. / 45 * Hm10 / x2
                 - 64. / 9 * Hm10 / x + 128. / 9 * Hm10 * x
                 + 128. / 9 * Hm10 * x2 + 64. / 15 * Hm10 * x3 - 256. / 5 * H0
                 + 128. / 45 * H0 / x - 8992. / 135 * H0 * x
                 + 8224. / 135 * H0 * x2 - 64. / 3 * H00 - 128 * H00 * x
                 + 128. / 3 * H00 * x2 - 64. / 15 * H00 * x3 - 64 * H000 * x
                 - 224. / 9 * H1 + 160. / 27 * H1 / x - 32. / 9 * H1 * x
                 + 608. / 27 * H1 * x2 + 32. / 3 * H11 - 32. / 9 * H11 / x
                 - 64. / 9 * H11 * x2 - 352. / 9 * H01 * x - 64. / 9 * H01 * x2
                 + 32. / 3 * H011 * x);
}

//==========================================================================================//
//  Charge factors for the flavour topologies entering up to
//  three loops.
//
//  Table 2 of Ref. [arXiv:hep-ph/0504242v1]
//------------------------------------------------------------------------------------------//

// u, d, s, c, b, t
double charges[] = { 2. / 3., -1. / 3., -1. / 3., 2. / 3., -1. / 3., 2. / 3. };

double fl11g(int nf) {

    double eavg = 0., e2avg = 0.;

    for (int i = 0; i < nf; i++) {
        eavg += charges[i];
        e2avg += charges[i] * charges[i];
    }
    eavg /= nf;
    e2avg /= nf;

    return eavg * eavg / e2avg;
}

double fl11ps(int nf) {

    double eavg = 0., e2avg = 0.;

    for (int i = 0; i < nf; i++) {
        eavg += charges[i];
        e2avg += charges[i] * charges[i];
    }
    eavg /= nf;
    e2avg /= nf;

    return eavg * eavg / e2avg - 3 * eavg;
}

//==========================================================================================//
//  Massless gluon coefficient functions for F2 at O(alpha_s^2) for mu=Q.
//  Observe that this result is a parameterization of the exact (known but long)
//  result
//
//  Eq. (4.10) of Ref. [hep-ph/0504242v1]
//------------------------------------------------------------------------------------------//

double C2_g2_massless_param(double x, int nf) {

    double x1 = 1. - x;
    double L0 = log(x);
    double L1 = log(x1);

    return nf
           * (58. / 9. * L1 * L1 * L1 - 24 * L1 * L1 - 34.88 * L1 + 30.586
              - (25.08 + 760.3 * x + 29.65 * L1 * L1 * L1) * x1
              + 1204 * x * L0 * L0 + L0 * L1 * (293.8 + 711.2 * x + 1043 * L0)
              + 115.6 * L0 - 7.109 * L0 * L0 + 70. / 9. * L0 * L0 * L0
              + 11.9033 * x1 / x);
}

//==========================================================================================//
//  Massless quark coefficient functions for F2 at O(alpha_s^2) for mu=Q.
//  Observe that this result is a parameterization of the exact (known but long)
//  result
//
//  Eq. (4.9) of Ref. [hep-ph/0504242v1]
//------------------------------------------------------------------------------------------//

double C2_ps2_massless_param(double x, int nf) {

    double x1 = 1 - x;
    double L0 = log(x);
    double L02 = L0 * L0;
    double L03 = L02 * L0;

    double L1 = log(x1);
    double L12 = L1 * L1;
    double L13 = L12 * L1;

    return nf
           * ((8. / 3 * L12 - 32. / 3 * L1 + 9.8937) * x1
              + (9.57 - 13.41 * x + 0.08 * L13) * x1 * x1 + 5.667 * x * L03
              - L02 * L1 * (20.26 - 33.93 * x) + 43.36 * x1 * L0 - 1.053 * L02
              + 40. / 9 * L03 + 5.2903 / x * x1 * x1);
}

//==========================================================================================//
//  Massless gluon coefficient functions for FL at O(alpha_s^2) for mu=Q.
//  Observe that this result is a parameterization of the exact (known but long)
//  result
//
//  Eq. (6) of Ref. [arXiv:hep-ph/0411112v2]
//------------------------------------------------------------------------------------------//

double CL_g2_massless_param(double x, int nf) {

    double L0 = log(x);
    double L02 = L0 * L0;

    double x1 = 1. - x;
    double L1 = log(x1);
    double L12 = L1 * L1;

    return nf
           * ((94.74 - 49.2 * x) * x1 * L12 + 864.8 * x1 * L1
              + 1161 * x * L1 * L0 + 60.06 * x * L02 + 39.66 * x1 * L0
              - 5.333 * (1. / x - 1));
}

//==========================================================================================//
//  Massless gluon coefficient functions for FL at O(alpha_s^2) for mu=Q.
//  Observe that this result is a parameterization of the exact (known but long)
//  result
//
//  Eq. (5) of Ref. [arXiv:hep-ph/0411112v2]
//------------------------------------------------------------------------------------------//

double CL_ps2_massless_param(double x, int nf) {

    double L0 = log(x);
    double L02 = L0 * L0;

    double x1 = 1. - x;
    double L1 = log(x1);

    return nf
           * ((15.94 - 5.212 * x) * x1 * x1 * L1 + (0.421 + 1.520 * x) * L02
              + 28.09 * x1 * L0 - (2.370 / x - 19.27) * x1 * x1 * x1);
}

//==========================================================================================//
//  Massless gluon coefficient functions for F2 at O(alpha_s^3) for mu=Q.
//  Observe that this result is a parameterization of the exact (known but long)
//  result. The term fl_g_11 is put to zero for the reason explained in page 15
//  of arXiv:1205.5727
//
//  Eq. (4.13) of Ref. [hep-ph/0504242v1]
//------------------------------------------------------------------------------------------//

double C2_g3_massless_param(
    double x, int nf
) { // remember that there is a delta(x1) that has been omitted

    // double fl_g_11 = fl11g(nf) ;

    double x2 = x * x;
    // double x3 = x2 * x;

    double x1 = 1 - x;
    double L0 = log(x);
    double L1 = log(x1);

    double L02 = L0 * L0;
    double L03 = L02 * L0;
    double L04 = L03 * L0;
    double L05 = L04 * L0;

    double L12 = L1 * L1;
    double L13 = L12 * L1;
    double L14 = L13 * L1;
    double L15 = L14 * L1;

    double c_nf =
        (966. / 81. * L15 - 1871. / 18. * L14 + 89.31 * L13 + 979.2 * L12
         - 2405 * L1 + 1372 * x1 * L14 - 15729 - 310510 * x + 331570 * x2
         - 244150 * x * L02 - 253.3 * x * L05 + L0 * L1 * (138230 - 237010 * L0)
         - 11860 * L0 - 700.8 * L02 - 1440 * L03 + 4961. / 162. * L04
         - 134. / 9. * L05
         - (6362.54 + 932.089 * L0)
               / x // there's a typo in the paper: -932.089 -> +932.089
        );         // there is + 0.625*delta(x1)

    double c_nf2 =
        (131. / 81. * L14 - 14.72 * L13 + 3.607 * L12 - 226.1 * L1 + 4.762
         - 190 * x - 818.4 * x2 - 4019 * x * L02 - L0 * L1 * (791.5 + 4646 * L0)
         + 739.0 * L0 + 418.0 * L02 + 104.3 * L03 + 809. / 81. * L04
         + 12. / 9. * L05 + 84.423 / x);

    /*double c_nf_fl = (
        3.211 * L12 + 19.04 * x * L1 + 0.623 * x1 * L13 - 64.47
        * x + 121.6 * x2 - 45.82 * x3 - x * L0 * L1 * (31.68 +
        37.24 * L0) + 11.27 * x2 * L03 - 82.40 * x * L0 - 16.08
        * x * L02 + 520./81. * x * L03 + 20./27. * x * L04
    );*/

    return (
        c_nf * nf + c_nf2 * nf * nf
        //+ c_nf_fl * nf * nf * fl_g_11
    );
}

//==========================================================================================//
//  Massless quark coefficient functions for F2 at O(alpha_s^3) for mu=Q.
//  Observe that this result is a parameterization of the exact (known but long)
//  result. The term fl_ps_11 is put to zero for the reason explained in page 15
//  of arXiv:1205.5727
//
//  Eq. (4.12) of Ref. [hep-ph/0504242v1]
//------------------------------------------------------------------------------------------//

double C2_ps3_massless_param(
    double x, int nf
) { // remember that there is a delta(x1) that has been omitted

    // double fl_ps_11 = fl11ps(nf) ;

    double x2 = x * x;

    double x1 = 1. - x;
    double L0 = log(x);
    double L1 = log(x1);

    double L02 = L0 * L0;
    double L03 = L02 * L0;
    double L04 = L03 * L0;
    double L05 = L04 * L0;

    double L12 = L1 * L1;
    double L13 = L12 * L1;
    double L14 = L13 * L1;

    double c_nf =
        ((856. / 81 * L14 - 6032. / 81 * L13 + 130.57 * L12 - 542 * L1 + 8501
          - 4714 * x + 61.5 * x2)
             * x1
         + L0 * L1 * (8831 * L0 + 4162 * x1) - 15.44 * x * L05 + 3333 * x * L02
         + 1615 * L0 + 1208 * L02 - 333.73 * L03 + 4244. / 81 * L04
         - 40. / 9 * L05 - 1. / x * (2731.82 * x1 + 414.262 * L0));

    double c_nf2 =
        ((-64. / 81 * L13 + 208. / 81 * L12 + 23.09 * L1 - 220.27 + 59.80 * x
          - 177.6 * x2)
             * x1
         - L0 * L1 * (160.3 * L0 + 135.4 * x1) - 24.14 * x * L03
         - 215.4 * x * L02 - 209.8 * L0 - 90.38 * L02 - 3568. / 243 * L03
         - 184. / 81 * L04 + 40.2426 * x1 / x);

    // double c_fl_nf = (
    //     (126.42 - 50.29 * x - 50.15 * x2) * x1 - 26.717
    //     - 9.075 * x * x1 * L1 - x * L02 * (101.8 + 34.79 * L0 + 3.070 * L02)
    //     + 59.59 * L0 - 320./81 * L02 * (5 + L0)
    // ) * x ;

    return (
        c_nf * nf + c_nf2 * nf * nf
        //+ c_fl_nf * fl_ps_11 * nf
    );
}

//==========================================================================================//
//  Massless gluon coefficient functions for FL at O(alpha_s^3) for mu=Q.
//  Observe that this result is a parameterization of the exact (known but long)
//  result. The term fl_ps_11 is put to zero for the reason explained in page 15
//  of arXiv:1205.5727
//
//  Eq. (10) of Ref. [arXiv:hep-ph/0411112v2]
//------------------------------------------------------------------------------------------//

double CL_g3_massless_param(
    double x, int nf
) { // remember that there is a delta(x1) that has been omitted

    // double fl_g_11 = fl11g(nf) ;

    double x2 = x * x;

    double x1 = 1. - x;
    double L0 = log(x);
    double L1 = log(x1);

    double L02 = L0 * L0;
    double L03 = L02 * L0;
    // double L04 = L03 * L0;

    double L12 = L1 * L1;
    double L13 = L12 * L1;
    double L14 = L13 * L1;

    double c_nf =
        ((144. * L14 - 47024. / 27. * L13 + 6319. * L12 + 53160. * L1) * x1
         + 72549. * L0 * L1 + 88238. * L02 * L1
         + (3709. - 33514. * x - 9533. * x2) * x1 + 66773. * x * L02
         - 1117. * L0 + 45.37 * L02 - 5360. / 27. * L03
         - (2044.70 * x1 + 409.506 * L0) / x);

    double c_nf2 =
        ((32. / 3. * L13 - 1216. / 9. * L12 - 592.3 * L1 + 1511. * x * L1) * x1
         + 311.3 * L0 * L1 + 14.24 * L02 * L1 + (577.3 - 729. * x) * x1
         + 30.78 * x * L03 + 366. * L0 + 1000. / 9. * L02 + 160. / 9. * L03
         + 88.5037 / x * x1);

    /*double c_nf_fl = (
        (
            -0.0105 * L13 + 1.55 * L12 + 19.72 * x * L1 - 66.745 * x
            + 0.615 * x2
        ) * x1
        + 20./27. * x * L04 + (280./81. + 2.26 * x) * x * L03
        - (15.4 - 2.201 * x) * x * L02 - (71.66 - 0.121 * x) * x * L0
    ) ;*/

    return (
        c_nf * nf + c_nf2 * nf * nf
        //+ c_nf_fl * nf * nf * fl_g_11
    );
}

//==========================================================================================//
//  Massless gluon coefficient functions for FL at O(alpha_s^3) for mu=Q.
//  Observe that this result is a parameterization of the exact (known but long)
//  result. The term fl_ps_11 is put to zero for the reason explained in page 15
//  of arXiv:1205.5727
//
//  Eq. (9) of Ref. [arXiv:hep-ph/0411112v2]
//------------------------------------------------------------------------------------------//

double CL_ps3_massless_param(
    double x, int nf
) { // remember that there is a delta(x1) that has been omitted

    // double fl_ps_11 = fl11ps(nf) ;

    // double x2=x*x;

    double x1 = 1. - x;
    double L0 = log(x);
    double L1 = log(x1);

    double L02 = L0 * L0;
    double L03 = L02 * L0;

    double L12 = L1 * L1;
    double L13 = L12 * L1;

    double c_nf =
        ((1568. / 27 * L13 - 3968. / 9 * L12 + 5124 * L1) * x1 * x1
         + (2184 * L0 + 6059 * x1) * L0 * L1 - (795.6 + 1036 * x) * x1 * x1
         - 143.6 * x1 * L0 + 2848. / 9 * L02 - 1600. / 27 * L03
         - (885.53 * x1 + 182.00 * L0) / x * x1);

    double c_nf2 =
        ((-32. / 9 * L12 + 29.52 * L1) * x1 * x1
         + (35.18 * L0 + 73.06 * x1) * L0 * L1 - 35.24 * x * L02
         - (14.16 - 69.84 * x) * x1 * x1 - 69.41 * x1 * L0 - 128. / 9 * L02
         + 40.239 / x * x1 * x1);

    /*double c_fl_nf = (
        ((107.0 + 321.05 * x - 54.62 * x2) * x1 - 26.717 + 9.773
        * L0 + (363.8 + 68.32 * L0) * x * L0 - 320./81 * L02 * (2
        + L0)) * x
    );*/

    return (
        c_nf * nf + c_nf2 * nf * nf
        //+ c_fl_nf * fl_ps_11 * nf
    );
}
