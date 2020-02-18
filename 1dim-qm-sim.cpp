#define _USE_MATH_DEFINES

// potential shape
// 0: box, 1: harm-osci, else: flat
#define MODE 0

#include <cmath>
#include <complex>
#include <iostream>
#include <limits>
#include <vector>
#include "1dim-qm-sim.h"

using namespace std;

int main () {
    static const double h  = 0.004; // h: integral slice
    static const double dt = pow (h, 2.0) / 4.0; // time slice
    static const int wSize = 4; // workspace size
    static const int tmax  = 12800; // maximum time
    static const int imax  = 4001; // - h * imax / 2 < x < h * imax / 2

    vector<double> v;
    vector<complex<double>> phi;
    vector< vector<complex<double>> >
        workspace ( wSize, vector<complex<double>> (imax, 0.0) );

    v.reserve (imax); // potential
    vInit (imax, h, v); // initialize potential

    phi.reserve (imax); // wave function
    phiInit (imax, h, phi); // initialize wave function

    for (int int_t = 0; int_t < tmax; int_t++) {
        if ( int_t % 10 == 0 ) { // output the result for every 10 loop
            printVandPhi(int_t, imax, h, v, phi);
        }
        step (dt, imax, h, v, phi, workspace);
    }
    return 0;
}

void vInit (int imax, double h, vector<double>& v) {
    double v0 = - pow (70.7 * M_PI, 2); // box depth
    double vWidth = 2.0; // box width
    double w = 60 * M_PI; // harmonic oscillator freq

    for (int i = 0; i < imax; i++) {
        double x = ( i - (imax-1) / 2 ) * h;
#if MODE == 0 // finite depth box potential
        if (fabs(x) <= vWidth / 2) {
            v.emplace_back( v0 );
        } else {
            v.emplace_back( 0.0 );
        }
#elif MODE == 1 // harmonic oscillator
        double v_i = pow(w, 2) * pow(x, 2) / 4.0;
        // potentialが大きすぎると誤差がblow-upするので丸めておく
        double v_max = 3e5;
        if (v_i > v_max) {
            v_i = v_max;
        }
        v.emplace_back( v_i );
#else // flat
        v.emplace_back( 0.0 );
#endif
    }
}

void phiInit (int imax, double h, vector<complex<double>>& phi) {
    // parameters of the Gaussian wave packet
    double deltax = 0.035; // initial width
    double x0 = - 3; // initial position
    double p0 = M_PI * 40; // initial momentum

    phi.emplace_back (0.0);

    for (int j = 1; j < imax-1; j++) {
        double x = ( j - (imax-1) / 2 ) * h;
        complex<double> phi_i = 1.0 / sqrt( sqrt( 2.0 * M_PI * pow(deltax,2) ) )
            * exp( - pow(x-x0, 2) / (4.0 * pow(deltax, 2) ) + 1.0i * p0 * x );
        phi.emplace_back (phi_i);
    }

    phi.emplace_back (0.0);
}

void step(double dt, double imax, double h, vector<double>& v, vector<complex<double>>& phi, vector< vector<complex<double>> >& ws) {
    vector<double> factor = {0.5, 0.5, 1.0};

            // Note the fixed end boundary condition
            //  (j = 0, imax-1 components of ws[0,1,2,3] are not updated)
            // and set imax large enough to avoid unphysical reflection from end points
    // set ws[0]
    for (int j=1; j < imax-1; j++) {
        complex<double> ws0_j = 1.0i * dt
        * (
            1.0 / pow(h, 2)
            * (phi[j+1] - 2.0 * phi[j] + phi[j-1])
            - v[j] * phi[j]
        );
        ws[0][j] = ws0_j;
    }

    // set ws[1], ws[2], ws[3]
    for (int i=1; i <= 3; i++) {
        for (int j=1; j < imax-1; j++) {
            complex<double> ws_ij = 1.0i * dt
            * (
                1.0 / pow(h, 2)
                * (
                      ( phi[j+1] + ws[i-1][j+1] * factor[i-1] )
                    - 2.0 * ( phi[j] + ws[i-1][j] * factor[i-1] )
                    + ( phi[j-1] + ws[i-1][j-1] * factor[i-1] )
                )
                - v[j] * ( phi[j] + ws[i-1][j] * factor[i-1] )
            );
            ws[i][j] = ws_ij;
        }
    }

    // update phi using Runge-Kutta
    for (int j = 0; j < imax; j++) {
        phi[j] = phi[j] 
        + (ws[0][j] + 2.0 * ws[1][j] + 2.0 * ws[2][j] + ws[3][j] ) / 6.0;
    }
}

void printVandPhi (int int_t, int imax, double h, vector<double>& v, vector<complex<double>>& phi) {
    cout << int_t << endl;
    cout << endl << endl;

    // print wave function
    for (int j = 0; j < imax; j++) {
        double x = ( j - (imax-1) / 2 ) * h;
        printf ("%12.3f, %12.3f\n", x, abs(phi[j]) );
    }
    cout << endl << endl;

    // print potential
    for (int j = 0; j < imax; j++) {
        double x = ( j - (imax-1) / 2 ) * h;
        printf ("%12.3f, %9.1f\n", x, v[j] );
    }
    cout << endl << endl;
}
