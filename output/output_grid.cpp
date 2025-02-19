// See README.md for compilation suggestions

#include "adani/adani.h"

#include <cmath>
#include <ctime>
#include <fstream>
#include <iostream>
#include <sstream>
#include <string.h>
#include <thread>
#include <vector>

using namespace std;

std::string print_time(time_t seconds);

void loopFunction(
    int start, int end, vector<double> &xvec, vector<double> &Qvec,
    ApproximateCoefficientFunction &app, double **res_c, double **res_h,
    double **res_l, double m, int nf, double mufrac
);

int main(int argc, char **argv) {

    if (argc != 6) {
        cout << "ERROR!" << endl
             << "Usage: ./output_grid.exe mufrac = mu/Q m nf channel threads"
             << "Exiting..." << endl;
        return -1;
    }

    double mufrac = atof(argv[1]);
    double m = atof(argv[2]);

    int nf = atoi(argv[3]);

    char kind = argv[4][0];
    char channel = argv[4][1];

    int numThreads = atoi(argv[5]);

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

    int rows = x.size();
    int cols = Q.size();

    double **res_c = new double *[rows];

    for (int i = 0; i < rows; ++i) {
        res_c[i] = new double[cols];
    }

    double **res_h = new double *[rows];

    for (int i = 0; i < rows; ++i) {
        res_h[i] = new double[cols];
    }

    double **res_l = new double *[rows];

    for (int i = 0; i < rows; ++i) {
        res_l[i] = new double[cols];
    }

    string hs_version;

    if (channel == 'q')
        hs_version = "exact";
    else
        hs_version = "abmp";

    ApproximateCoefficientFunction Approx =
        ApproximateCoefficientFunction(3, kind, channel, true, hs_version);

    int iterationsPerThread = rows / numThreads;

    std::vector<std::thread> threads;

    time_t starting_time = time(NULL);

    for (int i = 0; i < numThreads; ++i) {
        int start = i * iterationsPerThread;
        int end = (i == numThreads - 1) ? rows : (i + 1) * iterationsPerThread;
        threads.emplace_back(
            loopFunction, start, end, ref(x), ref(Q), ref(Approx), res_c, res_h,
            res_l, m, nf, mufrac
        );
    }

    for (std::thread &thread : threads) {
        thread.join();
    }

    time_t ending_time = time(NULL);

    cout << "Total running time is " << print_time(ending_time - starting_time)
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

    for (int j = 0; j < cols; j++) {
        for (int i = 0; i < rows; i++) {
            output_c << res_c[i][j] << "   ";
            output_h << res_h[i][j] << "   ";
            output_l << res_l[i][j] << "   ";
        }
        output_c << endl;
        output_h << endl;
        output_l << endl;
    }

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

void loopFunction(
    int start, int end, vector<double> &xvec, vector<double> &Qvec,
    ApproximateCoefficientFunction &app, double **res_c, double **res_h,
    double **res_l, double m, int nf, double mufrac
) {
    double x, Q, m2Q2, mu, m2mu2;
    for (int i = start; i < end; i++) {
        for (int j = 0; j < int(Qvec.size()); j++) {
            x = xvec[i];
            Q = Qvec[j];
            mu = mufrac * Q;
            m2Q2 = m * m / (Q * Q);
            m2mu2 = m * m / (mu * mu);
            Value res = app.fxBand(x, m2Q2, m2mu2, nf);
            res_c[i][j] = res.GetCentral();
            res_h[i][j] = res.GetHigher();
            res_l[i][j] = res.GetLower();
            // std::cout << "Thread ID: " << std::this_thread::get_id() <<"
            // x="<< x << " Q=" << Q << " mu=" << mu << " nf" << nf << " res_c="
            // << res_c[i][j] << endl;
        }
    }
}
