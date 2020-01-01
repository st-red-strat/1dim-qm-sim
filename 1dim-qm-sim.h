#include <iostream>
#include <complex>
#include <vector>

using namespace std;

void vInit (int imax, double h, vector<double>& v);

void phiInit (int imax, double h, vector<complex<double>>& phi);

void step (double dt, double imax, double h, vector<double>& v, vector<complex<double>>& phi, vector< vector<complex<double>> >& ws);

void printVandPhi (int int_t, int imax, double h, vector<double>& v, vector<complex<double>>& phi);
