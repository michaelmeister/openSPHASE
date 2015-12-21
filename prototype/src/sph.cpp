#include "sph.h"

#define _USE_MATH_DEFINES

#include "kernels.h"
#include "scene2.h"
#include <stdio.h>
#include <limits>
#include <math.h>
#include <array>

#include <boost/foreach.hpp>
#include <particle2.h>


void SPH::init() {
    has_second_phase = false;
#pragma omp parallel for
    for (long long i = 0; i < (long long)particles->size(); i++) {
        Particle2 *p = particles->at(i);
        if (!((p->is_inflow()) || (p->is_inflownew()))) {
            if (!(p->is_second_phase()))
                p->rho = configuration.rho0_1;
            else {
                p->rho = configuration.rho0_2;
                has_second_phase = true;
            }
        }

        p->volume = configuration.volume;
        p->mass = configuration.volume*p->rho;
    }
}

Sceneparser::Sph<number> SPH::getConfiguration() {
    return configuration;
}

void SPH::setConfiguration(Sceneparser::Sph<number> conf) {
    configuration = conf;
    init();
}

void SPH::update_kernels(const Kernel *kernel) {
#pragma omp parallel for
    for (long long i = 0; i < (long long)particles->size(); i++) {
        Particle2 *p_i = particles->at(i);
        BOOST_FOREACH(NeighbourData &nd, neighbours->at(i)) {
            Particle2 *p_j = particles->at(nd.j);
            nd.update_kernel(kernel, p_i->pos - p_j->pos - nd.displ);
        }
    }
}

void SPH::update_kernels_f_only(const Kernel *kernel) {
#pragma omp parallel for
    for (long long i = 0; i < (long long)particles->size(); i++) {
        Particle2 *p_i = particles->at(i);
        BOOST_FOREACH(NeighbourData &nd, neighbours->at(i)) {
            Particle2 *p_j = particles->at(nd.j);
            nd.update_kernel_f_only(kernel, p_i->pos - p_j->pos - nd.displ);
        }
    }
}

void SPH::update_kernels_df_only(const Kernel *kernel) {
#pragma omp parallel for
    for (long long i = 0; i < (long long)particles->size(); i++) {
        Particle2 *p_i = particles->at(i);
        BOOST_FOREACH(NeighbourData &nd, neighbours->at(i)) {
            Particle2 *p_j = particles->at(nd.j);
            nd.update_kernel_df_only(kernel, p_i->pos - p_j->pos - nd.displ);
        }
    }
}

void SPH::move_particles(const Scene2 &scene)
{
#pragma omp parallel for
    for (long long i = 0; i < (long long)particles->size(); i++) {

        Particle2 *p_i = particles->at(i);
        if (p_i->is_reactor2()) {
            p_i->v.clear();
            continue;
        }
        if (!p_i->is_wall()) {
            if ((!p_i->is_inflownew() && !p_i->is_inflow()) || !scene.isInflowActive()) {
                p_i->pos += dt/number(2.0)*(p_i->v+p_i->v_avg);

                if (scene.outOfScene(p_i->pos)) p_i->set(Particlegenerator::Out);
            }
        } else {
            p_i->v = p_i->v_wall * damp;
            if (p_i->is_moving_wall()) {
                p_i->pos += dt/number(2.0)*(p_i->v_wall*damp);

                if (scene.outOfScene(p_i->pos)) p_i->set(Particlegenerator::Out);
            }
        }
    }
}

// stirrer uses a damped 'ghost force' based on the velocity difference to the expected velocity of the stirrer
void SPH::handle_stirrer(long long i)
{
    Particle2 *p_i = particles->at(i);
    const number proportional_gain = 20.0;    // stirrer acts as a P-controller

    if (configuration.stirrer.type == configuration.stirrer.ellipsoidal_rotator) {
        // use polar coordinates
        number radius = sqrt((p_i->pos.x-configuration.stirrer.origin_x)*(p_i->pos.x-configuration.stirrer.origin_x)+(p_i->pos.y-configuration.stirrer.origin_y)*(p_i->pos.y-configuration.stirrer.origin_y));
        number phi = atan2(p_i->pos.y-configuration.stirrer.origin_y,p_i->pos.x-configuration.stirrer.origin_x);
        number r_x = configuration.stirrer.length_x/2.0;    // semi-axis in x direction of ellipse
        number r_y = configuration.stirrer.length_y/2.0;
        // ellipse: x^2/a^2+y^2/b^2=1 -> check later if LHS<1
        number is_in_ellipse = (p_i->pos.x-configuration.stirrer.origin_x)*(p_i->pos.x-configuration.stirrer.origin_x)/(r_x*r_x)+(p_i->pos.y-configuration.stirrer.origin_y)*(p_i->pos.y-configuration.stirrer.origin_y)/(r_y*r_y);

        // check if inside ellipse but exclude small circle in the middle (factor 1/5.0 is heuristic)
        if( is_in_ellipse < 1.0 && r_y/5.0 ) {
            number velocity_stirrer = configuration.stirrer.velocity*(radius-0.5);
            number velocity_difference_x = + velocity_stirrer*sin(phi) - p_i->v.x;
            number velocity_difference_y = - velocity_stirrer*cos(phi) - p_i->v.y;
            p_i->f.x += damp * velocity_difference_x * proportional_gain;
            p_i->f.y += damp * velocity_difference_y * proportional_gain;
        }
    }

    else if (configuration.stirrer.type == configuration.stirrer.rectangular_fan) {
        // check if inside rectangular fan interaction zone
        number dist_x = abs(p_i->pos.x-configuration.stirrer.origin_x);
        number dist_y = abs(p_i->pos.y-configuration.stirrer.origin_y);
        if( dist_x < configuration.stirrer.length_x  &&  dist_y < configuration.stirrer.length_y ) {
            // vary target velocity from 1 to 0 with cosine in y direction (maximum in centre)
            number velocity_difference_y = configuration.stirrer.velocity*cos(dist_y/configuration.stirrer.length_y*M_PI/2.0) - p_i->v.y;
            // vary velocity difference as a cosine (from 0 to 1 in the argument) in x-direction
            p_i->f.y += damp * velocity_difference_y * proportional_gain;
        }
    }
}

void SPH::handle_stirrer2(long long i)
{
    Particle2 *p_i = particles->at(i);
    const number proportional_gain = 20.0;    // stirrer acts as a P-controller

    if (configuration.stirrer.type == configuration.stirrer.ellipsoidal_rotator) {
        // use polar coordinates
        number radius = sqrt((p_i->pos.x-configuration.stirrer2.origin_x)*(p_i->pos.x-configuration.stirrer2.origin_x)+(p_i->pos.y-configuration.stirrer2.origin_y)*(p_i->pos.y-configuration.stirrer2.origin_y));
        number phi = atan2(p_i->pos.y-configuration.stirrer2.origin_y,p_i->pos.x-configuration.stirrer2.origin_x);
        number r_x = configuration.stirrer2.length_x/2.0;    // semi-axis in x direction of ellipse
        number r_y = configuration.stirrer2.length_y/2.0;
        // ellipse: x^2/a^2+y^2/b^2=1 -> check later if LHS<1
        number is_in_ellipse = (p_i->pos.x-configuration.stirrer2.origin_x)*(p_i->pos.x-configuration.stirrer2.origin_x)/(r_x*r_x)+(p_i->pos.y-configuration.stirrer2.origin_y)*(p_i->pos.y-configuration.stirrer2.origin_y)/(r_y*r_y);

        // check if inside ellipse but exclude small circle in the middle (factor 1/5.0 is heuristic)
        if( is_in_ellipse < 1.0 && r_y/5.0 ) {
            number velocity_stirrer = configuration.stirrer2.velocity*(radius-0.5);
            number velocity_difference_x = + velocity_stirrer*sin(phi) - p_i->v.x;
            number velocity_difference_y = - velocity_stirrer*cos(phi) - p_i->v.y;
            p_i->f.x += damp * velocity_difference_x * proportional_gain;
            p_i->f.y += damp * velocity_difference_y * proportional_gain;
        }
    }

    else if (configuration.stirrer.type == configuration.stirrer.rectangular_fan) {
        // check if inside rectangular fan interaction zone
        number dist_x = abs(p_i->pos.x-configuration.stirrer2.origin_x);
        number dist_y = abs(p_i->pos.y-configuration.stirrer2.origin_y);
        if( dist_x < configuration.stirrer2.length_x  &&  dist_y < configuration.stirrer2.length_y ) {
            // vary target velocity from 1 to 0 with cosine in y direction (maximum in centre)
            number velocity_difference_y = configuration.stirrer2.velocity*cos(dist_y/configuration.stirrer2.length_y*M_PI/2.0) - p_i->v.y;
            // vary velocity difference as a cosine (from 0 to 1 in the argument) in x-direction
            p_i->f.y += damp * velocity_difference_y * proportional_gain;
        }
    }
}


/**
  * moves the simulation one step ahead
  */

/* CFL after adams2012 */
number SPH::step(const Scene2 &scene) {
#define CFL

#ifdef CFL
    number max_dt = 0.005;
    dt = 0.25*configuration.sd/(max_c+max_v);
    dt = (std::min)(dt, max_dt);
    number bfc = 0.25*sqrt(configuration.sd/g.length());
    dt = (std::min)(dt, bfc);
    number sd2 = configuration.sd*configuration.sd;

    number nu;
    if (has_second_phase)
        nu = (std::max)(configuration.nu1,configuration.nu2);
    else
        nu = configuration.nu1;

    number vc = 0.125*sd2/nu;
    dt = (std::min)(dt, vc);
    if (dt == 0) (dt = 0.00001);
#else
    dt = 0.00001;   //SPH unit is (s)
#endif

#pragma omp parallel for
    for (long long i = 0; i < (long long)particles->size(); i++) {
        Particle2 *p = particles->at(i);
        if (p->is_fluid())
            p->v += dt / number(2.0)*p->f;
    }

    XSPH();

    move_particles(scene);

    update_kernels(scene.kernel);

    continuity();

    XSPH();

    move_particles(scene);

    update_kernels(scene.kernel);

    extrapolation();

    momentum();

    static int filter_step = 9;
    if (filter_step++ == 10) {
        if (configuration.moving_least_squares)
            moving_least_squares();
        else if (configuration.shepard)
            shepard();
        filter_step = 0;
    }

    update_wall_particles();

    number max_c, max_v;
    max_c = max_v = 0.0;

#if _OPENMP >= 200205
#pragma omp parallel for reduction(max: max_v, max_c)
#endif
    for (long long i = 0; i < (long long)particles->size(); i++) {
        Particle2 *p = particles->at(i);
        if (!(p->is_fluid()))
            continue;
        p->v += dt/number(2.0)*p->f;
        max_v = (std::max)(max_v, p->v.length());
        max_c = (std::max)(max_c, p->c);
    }
    this->max_c = max_c;
    this->max_v = max_v;

    return dt;
}

void SPH::extrapolation()
{
    const number wall_dens_const1 = configuration.gamma1 / (configuration.rho0_1*configuration.c1*configuration.c1);

    /* extrapolation of wall particle velocity and pressure */
#pragma omp parallel for
    for (long long w_i = 0; w_i < (long long)particles->size(); w_i++) {
        Particle2 *p_i = particles->at(w_i);

        if (!(p_i->is_wall()))
            continue;
        Vector2 v_sum, rho_r_wf_Wij;
        number volume = 0.0;
        number p_Wij = 0.0;
        BOOST_FOREACH(const NeighbourData &nd, neighbours->at(w_i)) {
            Particle2 *p_j = particles->at(nd.j);

            if (p_j->is_wall())
                continue;

            v_sum += p_j->v * nd.kernel;
            p_Wij += p_j->p * nd.kernel;
            rho_r_wf_Wij += p_j->rho * nd.dir*nd.kernel;
            volume += nd.kernel;
        }
        if (volume > 0.0) {
            p_i->v = p_i->v_wall * 2 * damp - v_sum / volume;
            p_i->p = (p_Wij + g*damp*rho_r_wf_Wij) / volume;
        }
        else {
            p_i->v = p_i->v_wall * 2 * damp;
        }

        p_i->rho = configuration.rho0_1 * std::pow(1 + p_i->p * wall_dens_const1, 1 / configuration.gamma1);
    }
}

/* XSPH Numerical simulation of interfacial flows by smoothed particle hydrodynamics */
inline
void SPH::XSPH()
{
    if (configuration.epsilon_xsph == 0.0)
        return;

#pragma omp parallel for
    for (long long i = 0; i < (long long)particles->size(); i++) {
        Particle2 *p_i = particles->at(i);

        Vector2 v_avg_tmp;
        if (!(p_i->is_fluid())) {
            p_i->v_avg.clear();
        }
        else {
            BOOST_FOREACH(const NeighbourData &nd, neighbours->at(i)) {
                Particle2 *p_j = particles->at(nd.j);

                if ((p_j->is_wall()) || (i == nd.j)
                        || ((p_i->is_second_phase()) != (p_j->is_second_phase()))) continue;
                Vector2 v_ji = p_j->v - p_i->v;

                number rhoij_bar = p_j->rho + p_i->rho;
                v_avg_tmp += p_j->mass / rhoij_bar * v_ji * nd.kernel;
            }
            p_i->v_avg = configuration.epsilon_xsph * v_avg_tmp;
        }
    }
}


inline
number SPH::continuity()
{
#pragma omp parallel for
    for (long long i = 0; i < (long long)particles->size(); i++) {
        Particle2 *p_i = particles->at(i);

        if (!(p_i->is_fluid())) continue;
        number dsum = 0.0;
        BOOST_FOREACH(const NeighbourData &nd, neighbours->at(i)) {
            if (i == nd.j)
                continue;

            Particle2 *p_j = particles->at(nd.j);

            const Vector2 v_ij = p_i->v + p_i->v_avg - p_j->v - p_j->v_avg;
            dsum += p_j->volume * v_ij * nd.kernel_gradient;
        }
        assert(!std::isnan(p_i->rho));
        p_i->rho *= 1 + dsum * dt;               /** Continuity equation **/
        assert(!std::isnan(p_i->rho));

        if (p_i->is_second_phase()) {
            const number rho_norm = p_i->rho / configuration.rho0_2;
            p_i->c = configuration.c2 * rho_norm * rho_norm * rho_norm;  /** d p2/d rho2 **/
            p_i->p = eos_complex(p_i->rho, configuration.rho0_2, configuration.c2, configuration.gamma2);             /** Equation of state phase 2**/
        }
        else {
            const number rho_norm = p_i->rho / configuration.rho0_1;
            p_i->c = configuration.c1 * rho_norm * rho_norm * rho_norm;  /** d p1/d rho1 **/
            p_i->p = eos_complex(p_i->rho, configuration.rho0_1, configuration.c1, configuration.gamma1);             /** Equation of state phase 1**/
        }
    }
#pragma omp parallel for
    for (long long i = 0; i < (long long)particles->size(); i++) {
        Particle2 *p_i = particles->at(i);

        if (!(p_i->is_wall()))
            p_i->volume = p_i->mass / p_i->rho;
    }
}
#define sq(x) ((x)*(x))

/**
  * more complex equation of state after (lopez2010smoothed).
  *
  */
number SPH::eos_complex(number rho, number rho0, number c, number gamma) {
    return (rho0*sq(c)) / gamma * (std::pow(rho / rho0, gamma) - 1);
}

void SPH::getVelocityRange(double &min_v, double &max_v) const {
    double _max_v = 0;
    double _min_v = 1 << 20;
#if _OPENMP >= 200205
#pragma omp parallel for reduction(min: _min_v) reduction(max: _max_v)
#endif
    for (long long i = 0; i < (long long) particles->size(); i++) {
        Particle2 *p = particles->at(i);
        double l = p->v.length_squared();
        _max_v = (std::max)(_max_v, l);
        _min_v = (std::min)(_min_v, l);
    }
    min_v = std::sqrt(_min_v);
    max_v = std::sqrt(_max_v);
}

void SPH::getDensityRange(double &min_rho, double &max_rho) const {
    double _max_rho = 0;
    double _min_rho = 1 << 20;
#if _OPENMP >= 200205
#pragma omp parallel for reduction(min: _min_rho) reduction(max: _max_rho)
#endif
    for (long long i = 0; i < (long long) particles->size(); i++) {
        Particle2 *p = particles->at(i);
        double l = p->rho;
        _max_rho = (std::max)(_max_rho, l);
        _min_rho = (std::min)(_min_rho, l);
    }
    min_rho = _min_rho;
    max_rho = _max_rho;
}

void SPH::getPressureRange(double &min_p, double &max_p) const {
    double _max_p = 0;
    double _min_p = 1 << 20;
#if _OPENMP >= 200205
#pragma omp parallel for reduction(min: _min_p) reduction(max: _max_p)
#endif
    for (long long i = 0; i < (long long) particles->size(); i++) {
        Particle2 *p = particles->at(i);
        double l = p->p;
        _max_p = (std::max)(_max_p, l);
        _min_p = (std::min)(_min_p, l);
    }
    min_p = _min_p;
    max_p = _max_p;
}

/**
  * calculate the pressure gradient, viscous forces and adds the resulting
  * acceleration to the particle
  */
inline
void SPH::momentum()
{
    if (configuration.rho0_2 == 0)
        momentum_single_fluid();
    else
        momentum_multi_fluid();
}



inline
void SPH::momentum_single_fluid()
{
    const Vector2 g_damp = g*damp;

#pragma omp parallel for
    for (long long i = 0; i < (long long)particles->size(); i++) {
        Particle2 *p_i = particles->at(i);

        if (!(p_i->is_fluid()))
            continue;

        p_i->f = g_damp;

        if (configuration.stirrer.is_active) {
            handle_stirrer(i);
            handle_stirrer2(i);
        }

        Vector2 f_p, f_v;

        number vol_i2 = p_i->volume * p_i->volume;
        BOOST_FOREACH(const NeighbourData &nd, neighbours->at(i))
        {
            const long long j = nd.j;
            if (i == j)
                continue;

            Particle2 *p_j = particles->at(nd.j);

            const Vector2 &Wij = nd.kernel_gradient;
            assert(!(std::isnan(Wij.x) || std::isnan(Wij.y)));

            /* artificial viscosity after Monaghan and Gingolds */
            number cij_bar = (p_i->c + p_j->c);
            number nue_2 = 0.01 * configuration.sd * configuration.sd;
            const Vector2 v_ij = p_i->v - p_j->v;
            const Vector2 r_ij = nd.dir;
            number v_ij_r_ij = v_ij*r_ij;
            number PI_ij = 0.0;

            number rho_sum = p_i->rho + p_j->rho;

            if (!(p_j->is_wall()) && (v_ij_r_ij < 0)) {
                number mue_ij = configuration.sd * v_ij_r_ij / (r_ij * r_ij + nue_2);
                PI_ij = -configuration.alpha * cij_bar* mue_ij / rho_sum;
            }

            number vol_j2 = p_j->volume * p_j->volume;

            /* pressure gradient after adami2012 */
            f_p -= ((vol_i2 + vol_j2) / p_i->mass *
                    (p_j->rho * p_i->p + p_i->rho * p_j->p) / rho_sum + p_j->rho * p_j->volume * PI_ij) * Wij;

            /* shear viscous forces after adami2012 */
            if ((configuration.no_slip || !(p_j->is_wall()))) {
                number mu_ij = 2 * configuration.nu1*(p_i->rho * p_j->rho / rho_sum)*(vol_i2 + vol_j2);
                number diff_Wij = 1.0/nd.dist * nd.kernel_deriv;
                Vector2 visc_factor = mu_ij * v_ij / p_i->mass;
                f_v += visc_factor * diff_Wij;
            }
        }
        p_i->f += f_v + f_p;
    }
}

inline
void SPH::momentum_multi_fluid()
{
    const number wall_dens_const2 = configuration.gamma2 / (configuration.rho0_2*configuration.c2*configuration.c2);
    const number rho_diff_rel = fabs(configuration.rho0_1 - configuration.rho0_2) / (configuration.rho0_1 + configuration.rho0_2);
    const Vector2 g_damp = g*damp;

#pragma omp parallel for
    for (long long i = 0; i < (long long)particles->size(); i++) {
        Particle2 *p_i = particles->at(i);

        if (!(p_i->is_fluid()))
            continue;

        p_i->f = g_damp;

        Vector2 f_p;

        number nu_i = (p_i->is_second_phase()) ? configuration.nu2 : configuration.nu1;

        BOOST_FOREACH(const NeighbourData &nd, neighbours->at(i)) {
            if (i == nd.j)
                continue;

            Particle2 *p_j = particles->at(nd.j);

            const Vector2 &Wij = nd.kernel_gradient;
            assert(!(std::isnan(Wij.x) || std::isnan(Wij.y)));

            /* artificial viscosity after Monaghan (2005 and 2013),
             * requires WendlandC2 for physical interpretation of viscosity */
            const Vector2 v_ij = p_i->v - p_j->v;
            number PI_ij = 0.0;
            number nu_j = 0.0;
            if( p_j->is_fluid() && p_j->is_second_phase() ) nu_j = configuration.nu2;
            else nu_j = configuration.nu1;
            if (!(p_j->is_wall()))
                PI_ij = - number(16.0) * nu_i * nu_j * v_ij * nd.dir / ( nd.dist * configuration.sd * (nu_i * p_i->rho+nu_j * p_j->rho) );

            /* density calculation of wall particles for first phase already done in extrapolation method
               -> recalculate wall velocity only if interacting particle is second phase */
            number rho_j = p_j->rho;
            if (p_j->is_wall() && p_i->is_second_phase())
                rho_j = configuration.rho0_2 * std::pow(1 + p_j->p * wall_dens_const2, 1 / configuration.gamma2);

            /* repulsive force to cure instability at the phase interface (after monaghan2013) */
            number p_div_rho = (p_i->p + p_j->p) / (p_i->rho * rho_j);

            number R_ij = 0.0;
            if ((p_j->is_fluid()) && ((p_i->is_second_phase()) != (p_j->is_second_phase())))
                R_ij = configuration.epsilon_repulsion * rho_diff_rel * fabs(p_div_rho);

            /* pressure gradient after monaghan2013 */
            f_p -= rho_j* p_j->volume * (p_div_rho + R_ij + PI_ij)*Wij;
        }
        p_i->f += f_p;
    }
}

void SPH::update_wall_particles()
{
    #pragma omp parallel for
    for (long long i = 0; i < (long long)particles->size(); i++) {
        Particle2 *p_i = particles->at(i);

        if ((p_i->is_wall()))
            p_i->p = 0;
    }
}

void SPH::shepard()
{
#pragma omp parallel for
    for (long long i = 0; i < (long long)particles->size(); i++) {
        Particle2 *p_i = particles->at(i);

        if (p_i->is_fluid()) {
            number tmp_volume = 0.0;

            BOOST_FOREACH(const NeighbourData &nd, neighbours->at(i)) {
                Particle2 *p_j = particles->at(nd.j);
                if (!(p_j->is_fluid()) || ((p_i->is_second_phase()) != (p_j->is_second_phase())))
                    continue;
                tmp_volume += p_j->volume * nd.kernel;
            }

            number rho_i_new = 0.0;

            BOOST_FOREACH(const NeighbourData &nd, neighbours->at(i)) {
                Particle2 *p_j = particles->at(nd.j);
                if (!(p_j->is_fluid()) || ((p_i->is_second_phase()) != (p_j->is_second_phase())))
                    continue;
                rho_i_new += p_j->volume * p_j->rho * nd.kernel;
            }
            p_i->rho = rho_i_new / tmp_volume;
            p_i->volume = p_i->mass/p_i->rho;
        }
    }

#pragma omp parallel for
    for (long long i = 0; i < (long long)particles->size(); i++) {
        Particle2 *p_i = particles->at(i);
        p_i->volume = p_i->mass/p_i->rho;
    }
}


number SPH::shepard_for_particle_i(long long i)
{
    Particle2 *p_i = particles->at(i);

    if (p_i->is_fluid()) {
        number tmp_volume = 0.0;

        BOOST_FOREACH(const NeighbourData &nd, neighbours->at(i)) {
            Particle2 *p_j = particles->at(nd.j);
            if (!(p_j->is_fluid()) || ((p_i->is_second_phase()) != (p_j->is_second_phase())))
                continue;
            tmp_volume += p_j->volume * nd.kernel;
        }

        number rho_i_new = 0.0;

        BOOST_FOREACH(const NeighbourData &nd, neighbours->at(i)) {
            Particle2 *p_j = particles->at(nd.j);
            if (!(p_j->is_fluid()) || ((p_i->is_second_phase()) != (p_j->is_second_phase())))
                continue;
            rho_i_new += p_j->volume * p_j->rho * nd.kernel;
        }
        return /*p_i->rho = */rho_i_new / tmp_volume;
    }
}

// density correction algorithm of first order after "State-of-the-art of classical SPH for free-surface flows"
void SPH::moving_least_squares()
{
    std::vector<number> rho_i_new(particles->size());

    #pragma omp parallel for
    for (long long i = 0; i < (long long)particles->size(); i++) {
        Particle2 *p_i = particles->at(i);

        if (p_i->is_fluid()) {
            number _rho_i_new = 0.0;
            bool MLS_failed = false;
            BOOST_FOREACH(const NeighbourData &nd, neighbours->at(i)) {
                Particle2 *p_j = particles->at(nd.j);

                if (!(p_j->is_fluid()) || ((p_i->is_second_phase()) != (p_j->is_second_phase())))
                   continue;
                _rho_i_new += p_j->mass * compute_W_MLS(i, nd, MLS_failed);
            }
            if( MLS_failed == false)
                rho_i_new[i] = _rho_i_new;
            else    // variable-rank algorithm (Tilts 1999)
                rho_i_new[i] = shepard_for_particle_i(i);
        }
    }

    // assign new values now
    #pragma omp parallel for
    for (long long i = 0; i < (long long)particles->size(); i++) {
        Particle2 *p_i = particles->at(i);
                if (p_i->is_fluid()) {
                    p_i->rho = rho_i_new[i];
                    p_i->volume = p_i->mass/p_i->rho;
                }
    }
}

inline
number SPH::compute_W_MLS(long long i, const NeighbourData &nd, bool &MLS_failed)
{
    Particle2 *p_i = particles->at(i);
    Particle2 *p_j = particles->at(nd.j);
    Vector2 x_i = p_i->pos;
    Vector2 x_j = p_j->pos;
    std::array<number,3> beta;

    compute_beta(beta, i, MLS_failed);

    number W_MLS = ( beta[0] + beta[1]*(x_i.x - x_j.x) + beta[2]*(x_i.y - x_j.y) ) * nd.kernel;

    return W_MLS;
}

inline
void SPH::compute_beta(std::array<number,3> &beta, long long i, bool &MLS_failed)
{
    std::array<number,6> A;   // the full matrix is not stored due to symmetry -> see definition in function 'compute_A_tilde'
    for (int j=0;j<6;j++) A[j] = 0.0;

    compute_A(A,i);

    // hard-coded matrix inversion using the analytic formula
    number determinantOfA = A[0]*( A[3]*A[5] - A[4]*A[4] ) - A[1]*( A[5]*A[1] - A[4]*A[2] ) + A[2]*( A[1]*A[4] - A[3]*A[2] );
    beta[0] = + ( A[3]*A[5] - A[4]*A[4] ) / determinantOfA;
    beta[1] = - ( A[5]*A[1] - A[4]*A[2] ) / determinantOfA;
    beta[2] = + ( A[1]*A[4] - A[3]*A[2] ) / determinantOfA;

    if( std::abs(determinantOfA) < 1E-16 )
         MLS_failed = true;
}

inline
void  SPH::compute_A(std::array<number,6> &A, long long i)
{
    Particle2 *p_i = particles->at(i);
    std::array<number,6> A_tilde;

    BOOST_FOREACH(const NeighbourData &nd, neighbours->at(i)) {
        Particle2 *p_j = particles->at(nd.j);
        if (!(p_j->is_fluid()) || ((p_i->is_second_phase()) != (p_j->is_second_phase())))
            continue;
        compute_A_tilde(A_tilde, p_i->pos,p_j->pos);

        for(int i=0; i<6; i++)
            A[i] += nd.kernel * p_j->volume * A_tilde[i];
    }
}

inline
void SPH::compute_A_tilde(std::array<number,6> &A_tilde, Vector2 x_i, Vector2 x_j)
{
    // Let X be the A_tilde matrix after the definition in the paper.
    // Then due to symmetry store only A_tilde as [X(0,0),X(0,1),X(0,2),X(1,1),X(1,2),X(2,2)]

    A_tilde[0] = 1.0;
    A_tilde[1] = (x_i.x - x_j.x);
    A_tilde[2] = (x_i.y - x_j.y);
    A_tilde[3] = (x_i.x - x_j.x) * (x_i.x - x_j.x);
    A_tilde[4] = (x_i.x - x_j.x) * (x_i.y - x_j.y);
    A_tilde[5] = (x_i.y - x_j.y) * (x_i.y - x_j.y);
}
