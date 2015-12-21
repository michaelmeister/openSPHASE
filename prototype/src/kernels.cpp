#include "kernels.h"
#define _USE_MATH_DEFINES
#include "math.h"
#include <cmath>

/**
  * spiky pressure kernel
  */
inline double kernel_spiky(double r, double h) {
    if (r > h) {
        return 0.0;
    }
    double hl = (h - r);
    double h5 = h*h*h*h*h;
    double result = 10.0/(h5*M_PI)*(hl*hl*hl);
    assert(!std::isnan(result));
    return result;
}

/**
  * spiky pressure kernel
  */
inline double kernel_spiky_derivative(double r, double h) {
    if (r > h) {
        return 0.0;
    }
    double hl = (h - r);
    double h5 = h*h*h*h*h;
    double result = (-30.0/(M_PI*h5))*(hl*hl);
    assert(!std::isnan(result));
    return result;
}

/**
  * Poly 6 kernel
  */
inline double kernel_dens(double r, double h) {
    if (r > h) {
        return 0;
    }
    double r2 = r*r;
    double h2 = h*h;
    double h8 = h2*h2*h2*h2;
    double h2r2 = (h2 - r2);
    double result = 4.0/(M_PI*h8)*(h2r2*h2r2*h2r2);
    assert(!std::isnan(result));
    return result;
}

inline double kernel_dens_derivative(double r, double h) {
    if (r > h) {
        return 0;
    }
    double r2 = r*r;
    double h2 = h*h;
    double h8 = h2*h2*h2*h2;
    double h2r2 = (h2 - r2);
    double result = -24.0*r/(M_PI*h8)*(h2r2*h2r2);
    assert(!std::isnan(result));
    return result;
}

inline double kernel_dens_laplacian(double r, double h) {
    if (r > h) {
        return 0;
    }
    double r2 = r*r;
    double r4 = r2*r2;
    double h2 = h*h;
    double h4 = h2*h2;
    double h8 = h4*h4;

    //double result = 24.0 / (M_PI*h8) * (1.0 - x)*(5.0*x - 1);
    double result = -24.0 / (M_PI*h8) * (h4 - 6.0*h2*r2 + 5*r4);
    assert(!std::isnan(result));
    return result;
}

/**
  * -\frac{5 \left(h^4-4 h r^3+3 r^4\right)}{3 h^5 \pi  r^2}
  * -((5 (h^4-4 h r^3+3 r^4))/(3 h^5 \[Pi] r^2))
  */

inline double kernel_visc_laplacian(double r, double h) {
    if (r > h) {
        return 0;
    }
    /*double r2 = r*r;
    double r3 = r2*r;
    double r4 = r2*r2;*/
    double h2 = h*h;
    double h4 = h2*h2;
    double h5 = h4*h;
    double h6 = h5*h2;

    return 45/(M_PI*h6)*(h-r);
    //return -5.0*(h4-4*h*r3+3*r4)/(3*h5*M_PI*r2);
}

inline double kernel_wendland2(double r, double h) {                     /** Wendland C2 */
    if (r > h) {
        return 0;
    }
    double q = r/h;
    double oneq = 1.0 - q;
    double oneq2 = oneq*oneq;
    double oneq4 = oneq2*oneq2;
    double h2 = h*h;

    return 7.0/(h2*M_PI)*oneq4*(1+4.0*q);
}

inline double kernel_wendland2_deriv(double r, double h) {
    if (r > h) {
        return 0;
    }
    double q = r/h;
    double oneqh = (q - 1.0) / h;
    double oneqh2 = oneqh*oneqh;

    return 140.0 * oneqh2 * oneqh * M_1_PI * q;
}

inline double kernel_wendland4(double r, double h) {                     /** Wendland C4 */
    double q = r/h;
    if (q <= 1){
        double h2 = h*h;
        double q2 = q*q;

        double oneq = 1.0 - q;
        double oneq2 = oneq*oneq;
        double oneq4 = oneq2*oneq2;

        return 9.0/(h2*M_PI)*oneq4*oneq2*(1+6.0*q+35.0/3.0*q2);
    }
    else
        return 0;
}

inline double kernel_wendland4_deriv(double r, double h) {
    double q = r/h;
    if (q <= 1){
        double oneq = q - 1.0;
        double oneq2 = oneq*oneq;

        double oneqh = oneq / h;
        double oneqh3 = oneqh*oneqh*oneqh;

        return 168.0 / (M_PI)*q*oneqh3*oneq2*(1 + 5.0*q);
    }
    else
        return 0;
}

inline double kernel_wendland6(double r, double h) {                     /** Wendland C6 */
    double q = r/h;
    if (q <= 1){
        double q2 = q*q;
        double h2 = h*h;

        double oneq = 1.0 - q;
        double oneq2 = oneq*oneq;
        double oneq4 = oneq2*oneq2;
        double oneq8 = oneq4*oneq4;

        return 78.0 / (7.0*h2*M_PI)*oneq8*(1 + 8.0*q + 25.0*q2 + 32.0*q2*q);
    }
    else
        return 0;
}

inline double kernel_wendland6_deriv(double r, double h) {
    double q = r/h;
    if (q <= 1){
        double q2 = q*q;
        double oneq = q - 1.0;
        double oneq2 = oneq*oneq;
        double oneq4 = oneq2*oneq2;

        double oneqh = oneq / h;
        double oneqh3 = oneqh*oneqh*oneqh;

        return 1716.0 / (7.0*M_PI)*q*oneqh3*oneq4*(1 + 7.0*q + 16.0*q2);
    }
    else
        return 0;
}

inline double kernel_cubic(double r, double h) {                        /** Cubic spline */
    double q = r/h;
    double h2 = h*h;


    if (q <= 0.5){
        double q2 = q*q;
        double q3 = q2*q;

        return 40.0/(7.0*h2*M_PI)*(1.0-6.0*q2+6.0*q3);
    }
    else if ((q > 0.5) && (q <= 1)){
        double oneq = 1.0 - q;
        double oneq2 = oneq*oneq;

        return 80.0 / (7.0*h2*M_PI)*oneq*oneq2;
    }
    else
        return 0;
}

inline double kernel_cubic_deriv(double r, double h) {
    double q = r/h;
    double h3 = h*h*h;
    if (q <= 0.5)
        return 240.0/(7.0*h3*M_PI)*q*(3.0*q-2.0);
    else if ((q > 0.5) && (q <= 1)){
        double oneq = q - 1.0;

        return -240.0 / (7.0*h3*M_PI)*oneq*oneq;
    }
    else
        return 0;
}

inline double kernel_quintic(double r, double h) {                    /** Quintic spline */
    double q = r/h;
    double h2 = h*h;
    double oneq = 1.0 - q;
    double oneq2 = oneq*oneq;
    double oneq4 = oneq2*oneq2;

    double two3q = (2.0 - 3.0*q);
    double two3q2 = two3q*two3q;
    double two3q4 = two3q2*two3q2;

    if (q <= 1.0 / 3.0){
        double one3q = (1.0 - 3.0*q);
        double one3q2 = one3q*one3q;
        double one3q4 = one3q2*one3q2;

        return 63.0 / (478.0*h2*M_PI)*(243.0*oneq4*oneq - 6.0*two3q4*two3q + 15.0*one3q4*one3q);
    }
    else if ((q > 1.0 / 3.0) && (q <= 2.0 / 3.0)){
        return 63.0 / (478.0*h2*M_PI)*(243.0*oneq4*oneq - 6.0*two3q4*two3q);
    }
    else if ((q > 2.0 / 3.0) && (q <= 1)){
        return 15309.0 / (478.0*h2*M_PI)*oneq*oneq4;
    }
    else
        return 0;
}

inline double kernel_quintic_deriv(double r, double h) {
    double q = r / h;
    double q2 = q*q;
    double q4 = q2*q2;

    double h3 = h*h*h;

    if (q <= 1.0/3.0)
        return -8505.0/(239.0*h3*M_PI)*q*(4.0-36.0*q2+45*q2*q);
    else if ((q > 1.0/3.0) && (q <= 2.0/3.0))
        return 2835.0/(478.0*h3*M_PI)*(5.0-84.0*q+270.0*q2-324.0*q2*q+135*q4);
    else if ((q > 2.0 / 3.0) && (q <= 1)){
        double oneq = q - 1.0;
        double oneq2 = oneq * oneq;
        double oneq4 = oneq2 * oneq2;

        return -76545.0 / (478.0*h3*M_PI)*oneq4;
    }
    else
        return 0;
}

struct Cubic : public Kernel {
    Cubic(double h) : h(h) {
    }

    double f(double r) const {
        return kernel_cubic(r, h);
    }

    double df(double r) const {
        return kernel_cubic_deriv(r, h);
    }

    double h;
};

struct Quintic : public Kernel {
    Quintic(double h) : h(h) {
    }

    inline
    double f(double r) const {
        return kernel_quintic(r, h);
    }

    inline
    double df(double r) const {
        return kernel_quintic_deriv(r, h);
    }

    double h;
};

struct Wendland2 : public Kernel {
    Wendland2(double h) : h(h) {
    }

    inline
    double f(double r) const {
        return kernel_wendland2(r, h);
    }

    inline
    double df(double r) const {
        return kernel_wendland2_deriv(r, h);
    }

    double h;
};

struct Wendland4 : public Kernel {
    Wendland4(double h) : h(h) {
    }

    double f(double r) const {
        return kernel_wendland4(r, h);
    }

    double df(double r) const {
        return kernel_wendland4_deriv(r, h);
    }

    double h;
};

struct Wendland6 : public Kernel {
    Wendland6(double h) : h(h) {
    }

    double f(double r) const {
        return kernel_wendland4(r, h);
    }

    double df(double r) const {
        return kernel_wendland6_deriv(r, h);
    }

    double h;
};

Kernel *Kernel::get(const std::string &name, double h) {
    if (name == "quintic") {
        return new Quintic(h);
    }
    else if (name == "cubic") {
        return new Cubic(h);
    }
    else if (name == "wendland2") {
        return new Wendland2(h);
    }
    else if (name == "wendland4") {
        return new Wendland4(h);
    }
    else if (name == "wendland6") {
        return new Wendland6(h);
    }
    else return NULL;
}
