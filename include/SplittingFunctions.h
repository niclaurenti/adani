#ifndef Split
#define Split

double pgq(double x);
double pqg(double x);
double pggreg(double x);
double pggsing(double x);
double Pgq0(double x);
double Pqg0(double x, int nf);

double Pgg0reg(double x);
double Pgg0loc(int nf);
double Pgg0sing(double x);

double Pqq0reg(double x);
double Pqq0loc();
double Pqq0sing(double x);

double Pgq1(double x);

double Pgg1reg(double x, int nf);
double Pgg1sing (double x, int nf);
double Pgg1loc(int nf);

#endif