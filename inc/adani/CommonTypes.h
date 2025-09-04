/*
 * =====================================================================================
 *
 *       Filename:  CommonTypes.h
 *
 *         Author:  Daniele Adani
 *
 *    Description:  ???
 *
 * =====================================================================================
 */

#ifndef Common_h
#define Common_h

#include <cstdint>
#include <string>

enum class DoubleIntegralMethod : uint8_t {
    Analytical, DoubleNumerical, MonteCarlo
};

enum class HighScaleVersion : uint8_t {
    Exact, GM, ABMP, KLMV
};

std::string to_string(HighScaleVersion hs_version);

#endif
