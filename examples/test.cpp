// On my machine it compiles with :
// g++ -o test.exe test.cpp -ladani -std=c++17 -stdlib=libc++
#include <adani/adani.h>
#include <iostream>
#include <cmath>

using namespace std ;

int main() {

    double x, logx, logx_min=-4, logx_max=0;
    double Q, logQ, logQ_min=1, logQ_max=2;
    int N=5;

    double dlogx=(logx_max - logx_min) / N;
    double dlogQ=(logQ_max - logQ_min) / N;

    double m = 4.92, mQ, mMu ;
    int nf = 4 ;

    for(int i=0; i<N; i++) {
        for(int j=0; j<N; j++) {

            logx = logx_min + i * dlogx;
            x = pow(10, logx);
            logQ = logQ_min + j * dlogQ;
            Q = pow(10, logQ);

            mQ = m * m / (Q * Q) ;
            mMu = mQ ;

            cout << "C2_g^(3)(x=" << x << ", Q="<< Q << ") = "<< C2m_g3_approximation(x, mQ, mMu, nf) << endl ;
            cout << "C2_ps^(3)(x=" << x << ", Q="<< Q << ") = "<< C2m_ps3_approximation(x, mQ, mMu, nf) << endl ;
            cout << "CL_g^(3)(x=" << x << ", Q="<< Q << ") = "<< CLm_g3_approximation(x, mQ, mMu, nf) << endl ;
            cout << "CL_ps^(3)(x=" << x << ", Q="<< Q << ") = "<< CLm_ps3_approximation(x, mQ, mMu, nf) << endl ;
        }
    }
    return 0;
}
