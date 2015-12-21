#include "../src/kernels.h"
#include <iostream>

using namespace std;

int main(int argc, char **argv) {
    int steps = 1000;
    double h = 1;
    double cut_off = 2*h;
    double dx = cut_off/steps;
    double value = 0;

    for (double r = 0; r <= cut_off; r += dx) {
        value += 2*M_PI*r*kernel_wendland(r, h)*dx;
    }
    cout << "wendland: \t" << value << endl;

    cut_off = h;
    dx = cut_off/steps;
    value = 0;

    for (double r = 0; r <= cut_off; r += dx) {
        value += 2*M_PI*r*kernel_spiky(r, h)*dx;
    }
    cout << "spiky: \t" << value << endl;

    value = 0;

    for (double r = 0; r <= cut_off; r += dx) {
        value += 2*M_PI*r*kernel_dens(r, h)*dx;
    }
    cout << "poly6: \t" << value << endl;

    /*cout.setf(ios::fixed,ios::floatfield);
    cout.precision(10);

    steps = 100;
    dx = 2.0/steps;
    double dh = 0.00001;
    for (int i = 0; i < steps; i++) {
        double at = dx*i;
        number deriv = kernel_wendland_deriv(at, 1.0);
        number num_deriv = (kernel_wendland(at+dh, 1.0) - kernel_wendland(at, 1.0))/dh;
        cout << deriv << "\t  " <<  num_deriv  << "\t error: " << (deriv - num_deriv) << "\n";
    }*/



    return 0;
}
