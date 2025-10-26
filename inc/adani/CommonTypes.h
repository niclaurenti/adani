/*
 * =====================================================================================
 *
 *       Filename:  CommonTypes.h
 *
 *         Author:  Daniele Adani
 *
 *    Description:  Noi parliamo di idee che danno animo a un gioco che poi è
 *                  quello che amiamo, ovvero il calcio
 *
 * =====================================================================================
 */

#ifndef Common_h
#define Common_h

#include <cstdint>
#include <string>

//==========================================================================================//
//  method used for solving the double integral in the log(mu^2/m^2)^2 terms of
//  the exact coefficient function
//------------------------------------------------------------------------------------------//

enum class DoubleIntegralMethod : uint8_t {
    Analytical,
    DoubleNumerical,
    MonteCarlo
};

//==========================================================================================//
//  version of aQg30
//------------------------------------------------------------------------------------------//

enum class HighScaleVersion : uint8_t { Exact, GM, ABMP, KLMV };

std::string to_string(HighScaleVersion hs_version);

//==========================================================================================//
//  struct approximation_parameters: parameters of the damping functions
//------------------------------------------------------------------------------------------//

struct approximation_parameters {
        double A;
        double B;
        double C;
        double D;
};

//==========================================================================================//
//  struct variation_parameters: parameters of the variation of the damping
//  functions
//------------------------------------------------------------------------------------------//

struct variation_parameters {
        double var1;
        double var2;
};

//==========================================================================================//
//  struct klmv_params: parameters for klmv approximation
//------------------------------------------------------------------------------------------//

struct klmv_params {
        double eta_exponent;
        double shift;
        double log_coeff;
        double log_pow;
        double const_coeff;
};

#endif
