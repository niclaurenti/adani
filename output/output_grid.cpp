// See README.md for compilation suggestions

#include "adani/adani.h"

#include <cmath>
#include <ctime>
#include <fstream>
#include <iostream>
#include <sstream>
#include <string.h>
#include <vector>

using namespace std;

std::string print_time(time_t seconds);

int main(int argc, char **argv) {

    if (argc != 5) {
        cout << "ERROR!" << endl
             << "Usage: ./output_grid.exe mufrac = mu/Q m nf channel"
             << "Exiting..." << endl;
        return -1;
    }

    double mufrac = atof(argv[1]);
    double m = atof(argv[2]);

    int nf = atoi(argv[3]);

    char kind = argv[4][0];
    char channel = argv[4][1];

    ifstream inputQ;
    inputQ.open("Q.txt");

    ifstream inputx;
    inputx.open("x.txt");

    std::vector<double> Q, x;

    double xtmp;
    inputx >> xtmp;
    while (!inputx.eof()) {
        x.push_back(xtmp);
        inputx >> xtmp;
    }

    double Qtmp;
    inputQ >> Qtmp;
    while (!inputQ.eof()) {
        Q.push_back(Qtmp);
        inputQ >> Qtmp;
    }

    // clang-format off
    cout << "                                                                                      " << endl;
    cout << "                    888       888          888                                        " << endl;
    cout << "                    888   o   888          888                                        " << endl;
    cout << "                    888  d8b  888          888                                        " << endl;
    cout << "                    888 d888b 888  .d88b.  888  .d8888b .d88b.  88888b.d88b.   .d88b. " << endl;
    cout << "                    888d88888b888 d8P  Y8b 888 d88P'   d88''88b 888 '888 '88b d8P  Y8b" << endl;
    cout << "                    88888P Y88888 88888888 888 888     888  888 888  888  888 88888888" << endl;
    cout << "                    8888P   Y8888 Y8b.     888 Y88b.   Y88..88P 888  888  888 Y8b.    " << endl;
    cout << "                    888P     Y888  'Y8888  888  'Y8888P 'Y88P'  888  888  888  'Y8888 " << endl;
    cout << "                                                                                      " << endl;
    cout << "                                                                                      " << endl;
    cout << "                    888                       d8888      888                   d8b    " << endl;
    cout << "                    888                      d88888      888                   Y8P    " << endl;
    cout << "                    888                     d88P888      888                          " << endl;
    cout << "                    888888 .d88b.          d88P 888  .d88888  8888b.  88888b.  888    " << endl;
    cout << "                    888   d88''88b        d88P  888 d88' 888     '88b 888 '88b 888    " << endl;
    cout << "                    888   888  888       d88P   888 888  888 .d888888 888  888 888    " << endl;
    cout << "                    Y88b. Y88..88P      d8888888888 Y88b 888 888  888 888  888 888    " << endl;
    cout << "                     Y888   Y88P       d88P     888   Y88888  Y888888 888  888 888    " << endl;
    cout << "                                                                                      " << endl;
    cout << "                                                                                      " << endl;
    cout << "                                                                                                     " << endl;
    cout << "                                          -=+:  -%+**  :=+#---                                       " << endl;
    cout << "                                       +%*:=-=#@@@@@@@@@@@@@@@@@@@@%@=+-                             " << endl;
    cout << "                                    :: *@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@*                           " << endl;
    cout << "                                 .= @@@@@@@@@@@@@@@@@@%@@@@@@@@@@@@@@@@@%@@#:                        " << endl;
    cout << "                              .-#@@@@@@@@@@@*:                   +@@@@@@@@%@@#:                      " << endl;
    cout << "                           .-=@%@@@@@@=.                           :@@@@@@@@@@@=-                    " << endl;
    cout << "                         ::+%@@@@@@#-:.                              #@@@@@@@%@@+=                   " << endl;
    cout << "                        ==*@@@@@@@%-..                               .=@@@@##@@@@-                   " << endl;
    cout << "                       %*%@@@@@@@@+:                                  .#@@@@@@@@@@=                  " << endl;
    cout << "                     =%*@@@@@@@@@%=:                                  .-@@@@@@@@@@*                  " << endl;
    cout << "                     %%@@@@@@@@@@#+.                                  .-@@@@@@@@@@@*                 " << endl;
    cout << "                    *%@@@@@@@@@@@*=.                                    +@@@@@@@@@%@-                " << endl;
    cout << "                   :%@@@@@@@@@@@##@@@@@@@#**=             -==#@@@@@@#:. :#@@@@@@@@@@@=               " << endl;
    cout << "                   %@@@@@@@@@@@%#%@@=.    .==-=.          :-      .=#%-. :%@@@@@@@@@@@-              " << endl;
    cout << "                  .@@@@@@@@@@@@+=....%@#@@@  -:=         #:@%*+#%#:   .  .=%@@@@@@@@@@               " << endl;
    cout << "                 *@@@@@@@@@@@@=:.             .                          .:+@@@@@@@@@@=              " << endl;
    cout << "                =@@@@@@@@@@@@@-.                                          :+@@@@@@@@@@%              " << endl;
    cout << "                %@@@@@@@@@@@@@-                                          .:*@ %%#@@@@@@              " << endl;
    cout << "                #@@@@@@@@@@@@@*             :                            -*##   =-@@@@@              " << endl;
    cout << "                *@@@@@@@@@@@#@#.=       .   .++%%: .-@::       .        .=*%#   .@@@@@@@             " << endl;
    cout << "                :@@@@@@@@@@#@#@==+:   =@%@@%%*@@@#.:@#=-.   %.  :*     :  -%@   %@@@@@@%             " << endl;
    cout << "                 *#@@@@@@@@@@#@@@##:  @@@@@+*+##+=-.*       :%@#@@  =+::.  =@@@@@@@@@@@#=            " << endl;
    cout << "                 = @@@@@@@@@@@@@@%*#+  **:#@@-              *   +# :.#*:=-:#@@@@@@@@@@@%++           " << endl;
    cout << "               *@@@@@@@@@@@@@@@@@@@*@#.=+    .                  @#++-===-##@@@@@@@@@@@@@=            " << endl;
    cout << "               =@@@@@@@@@@@@@@@@@@@#@#@*+  .-=-               :*@@@#@#%#**@@@@@@@@@@@@@@=            " << endl;
    cout << "                 *@@@@@@@@@@@@@@@@@@@@@@*. .==+##:  :.:.. ...  =@@@#@@@%@@+@@@@@@@@@@@@@@            " << endl;
    cout << "                  @@@@@@@@@@@@@@@@@@@@@@%@#*  -+.   -        :-=@@@@@@@@# =@@@@@@@@@@@@@@+           " << endl;
    cout << "                   =%@@@@@@@@@@@@@@@@@@@@@#@:               =@@@@@@@%@@@ @@@@@@@@@@@@@@@@            " << endl;
    cout << "                   ..@:@@@@@@@@@@@@@@@@@@@@-               .#%@@@@@@@%%@@@@ @@@@@@@@@@@@@@=          " << endl;
    cout << "                     %@@@@@@@@@@@@@+%@@@@@@@*=  -   .  . %=@+@@@@@@@*+*@@@+ @@@@@@@@@:@@@@           " << endl;
    cout << "                 .=#@@@@@@@@@@@@@@@@==#@@@@@@@@@@#@%=#+#@@@@@@@@@@*-::+@@@@@@@@@@@@@@*=              " << endl;
    cout << "                     .+@@@@@@@@@@@@@@@#=*%%@@@@@@@@@@@@@@@@@%#**+*-:.-%@@@@@@@@@@@@@                 " << endl;
    cout << "                       @@@@@@@@@@@@@@@@+=-:=+%%@@@@@@@@@=:. .-=-+=. .+@@@@@@@@@@@@@@*                " << endl;
    cout << "                        @@@@@@@@@@@@@@@@+-::...-==+=++=.    -=----: .*@@@@@@@@@@@@@@%                " << endl;
    cout << "                         @@@@@@@@@@@@@@@#.         :.---  .::.:-:  .-*%@@@@@@@@@@@@+                 " << endl;
    cout << "                           @@@@@@@@@@@@@@@.               ..       .=-=*-@@@@@@@@@@@                 " << endl;
    cout << "                             @@@@@@@@@@@@@@@                       :::=- -@@@@@@@@@                  " << endl;
    cout << "                                @@@@@@@@@@@@@@%                    ::-.  %@@@@@@@                    " << endl;
    cout << "                                  @@@@@@@@@@@@@@                   .+-  =@@@@@@@                     " << endl;
    cout << endl;
    cout << endl;
    // clang-format on

    cout << "Computation of the grid for the coefficient function C" << channel
         << " for m = " << m << " GeV, nf = " << nf << " and µ/Q = " << mufrac
         << endl;

    cout << "Size of the grid (x,Q) = (" << x.size() << "," << Q.size() << ")"
         << endl;

    string filename = "results/C_";
    filename.append(argv[4]);
    filename.append("_nf" + to_string(nf));

    ofstream output_c;
    ofstream output_h;
    ofstream output_l;

    output_c.open(filename + "_central.dat");
    if (!output_c.is_open()) {
        cout << "Problems in opening " << filename + "_central.dat" << endl;
        exit(-1);
    } else {
        cout << "Saving central grid in " << filename + "_central.dat" << endl;
    }

    output_h.open(filename + "_higher.dat");
    if (!output_h.is_open()) {
        cout << "Problems in opening " << filename + "_higher.dat" << endl;
        exit(-1);
    } else {
        cout << "Saving higher grid in " << filename + "_higher.dat" << endl;
    }

    output_l.open(filename + "_lower.dat");
    if (!output_l.is_open()) {
        cout << "Problems in opening " << filename + "_lower.dat" << endl;
        exit(-1);
    } else {
        cout << "Saving lower grid in " << filename + "_lower.dat" << endl;
    }

    Value res = Value(0, 0, 0);

    time_t starting_time = time(NULL);

    string hs_version;

    if (channel == 'q')
        hs_version = "exact";
    else
        hs_version = "abmp";

    ApproximateCoefficientFunction Approx =
        ApproximateCoefficientFunction(3, kind, channel, true, hs_version);

    for (double Q_ : Q) {
        for (double x_ : x) {
            double m2Q2, mu, m2mu2;
            m2Q2 = pow(m / Q_, 2);
            mu = mufrac * Q_;
            m2mu2 = pow(m / mu, 2);
            res = Approx.fxBand(x_, m2Q2, m2mu2, nf);
            output_c << res.GetCentral() << "   ";
            output_h << res.GetHigher() << "   ";
            output_l << res.GetLower() << "   ";
        }
        output_c << endl;
        output_h << endl;
        output_l << endl;
    }

    time_t ending_time = time(NULL);

    cout << "Total running time is " << print_time(ending_time - starting_time)
         << endl;

    output_c.close();
    output_h.close();
    output_l.close();

    inputQ.close();
    inputx.close();

    return 0;
}

std::string print_time(time_t seconds) {

    stringstream ss;

    int hour = (int)seconds / 3600;
    int minute = (int)(seconds - hour * 3600) / 60;
    int second = seconds - hour * 3600 - minute * 60;

    ss << hour << "h:" << minute << "m:" << second << "s";

    return ss.str();
}
