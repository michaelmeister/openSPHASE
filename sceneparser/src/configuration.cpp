#include "configuration.h"

namespace Sceneparser {

template <typename number>
Sph<number>::Stirrer::Stirrer()
    : is_active(false)
    , origin_x(0)
    , origin_y(0)
    , velocity(0)
    , length_x(0)
    , length_y(0)
    , type(ellipsoidal_rotator)
{
}

template <typename number>
Sph<number>::Sph()
    : volume(0)
    , rho0_1(1000)
    , rho0_2(0)
    , nu1(1.00e-6)
    , nu2(1.28e-4)
    , sd(1)
    , c1(30)
    , c2(30)
    , gamma1(7)
    , gamma2(7)
    , alpha(0)
    , epsilon_xsph(0)
    , epsilon_repulsion(0.08)
    , D_air(2.14e-5)
    , D_water(10.0e-3)
    , originx(0)
    , originy(0)
    , no_slip(true)
    , shepard(true)
    , moving_least_squares(false)
    , t_damp(0.0)
    , hydrostatic_initialization(false)
{
}

template <typename number>
Asm<number>::Asm()
    : sd(1.0)
    , width(0.0)
    , height(0.0)
    , ya(0.24)
    , yh(0.67)
    , fp(0.08)
    , ixb(0.086)
    , ixp(0.086)
    , ixe(0.06)
    , mh(6.0)
    , Ks(20.0)
    , Koh(0.2)
    , Kno(0.5)
    , bh(0.62)
    , ba(0.62)
    , eg(0.8)
    , eh(0.4)
    , kh(3.0)
    , rho7(3.0)
    , Kx(0.03)
    , ma(0.8)
    , Knh(1.0)
    , Koa(0.4)
    , ka(0.08)
    , SRT(12.0)
    , time_scaling(1.0)
{
}

template <typename number>
Configuration<number>::Configuration()
{

}

template <typename number>
Configuration<number>::~Configuration()
{

}

template class Configuration<float>;
template class Configuration<double>;
template struct Sph<float>;
template struct Sph<double>;
template struct Asm<float>;
template struct Asm<double>;

}
