#include "adani/CommonTypes.h"

#include <iostream>

std::ostream& operator<<(std::ostream& os, HighScaleVersion hs_version) {
    return os << to_string(hs_version);
}

std::string to_string(HighScaleVersion hs_version) {
    switch (hs_version) {
        case HighScaleVersion::Exact: return "Exact";
        case HighScaleVersion::GM: return "GM";
        case HighScaleVersion::ABMP: return "ABMP";
        case HighScaleVersion::KLMV: return "KLMV";
    }
}
