#include <iostream>
#include <fstream>
#include <cmath>

using namespace std;

int main(int argc, char** argv) {

    if(argc!=4) {
        cout << "ERROR\nUsage: ./grid_x.exe logx_min N_points filename\nExiting..." << endl;
        return -1;
    }

    double logx_min = atof(argv[1]);
    int N_points = atoi(argv[2]) ;

    ofstream output;
    output.open(argv[3]);

    double x, logx, logx_max=0;

    double dlog=(logx_max - logx_min)/N_points;

    for(int i=0; i <= N_points; i++) {

        logx=logx_min + i*dlog;

        x=pow(10, logx);

        output  << x << "   " ;
    }

    output << endl ;

    output.close();

    return 0;

}
