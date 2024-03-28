// See README.md for compilation suggestions

#include <adani/adani.h>
#include <cmath>
#include <iostream>

using namespace std;

int main() {

    double x, logx, logx_min = -4, logx_max = 0;
    double Q, logQ, logQ_min = 1, logQ_max = 2;
    int N = 5;

    double dlogx = (logx_max - logx_min) / N;
    double dlogQ = (logQ_max - logQ_min) / N;

    double m = 4.92, m2Q2, m2mu2;
    int nf = 4;

    ApproximateCoefficientFunction F2g(
        3, '2', 'g', true, "abmp", 1e-3, 1e-3, 1000, false, 25000
    );
    ApproximateCoefficientFunction FLg(
        3, 'L', 'g', true, "abmp", 1e-3, 1e-3, 1000, false, 25000
    );
    ApproximateCoefficientFunction F2q(
        3, '2', 'q', true, "exact", 1e-3, 1e-3, 1000, false, 25000
    );
    ApproximateCoefficientFunction FLq(
        3, 'L', 'q', true, "exact", 1e-3, 1e-3, 1000, false, 25000
    );

    for (int i = 0; i < N; i++) {
        for (int j = 0; j < N; j++) {

            logx = logx_min + i * dlogx;
            x = pow(10, logx);
            logQ = logQ_min + j * dlogQ;
            Q = pow(10, logQ);

            m2Q2 = m * m / (Q * Q);
            m2mu2 = m2Q2;

            cout << "C2_g^(3)(x=" << x << ", Q=" << Q
                 << ") = " << F2g.fxBand(x, m2Q2, m2mu2, nf) << endl;
            cout << "C2_ps^(3)(x=" << x << ", Q=" << Q
                 << ") = " << F2q.fxBand(x, m2Q2, m2mu2, nf) << endl;
            cout << "CL_g^(3)(x=" << x << ", Q=" << Q
                 << ") = " << FLg.fxBand(x, m2Q2, m2mu2, nf) << endl;
            cout << "CL_ps^(3)(x=" << x << ", Q=" << Q
                 << ") = " << FLq.fxBand(x, m2Q2, m2mu2, nf) << endl;
        }
    }
    return 0;
}
