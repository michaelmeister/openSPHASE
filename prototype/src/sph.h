#ifndef SPH_H
#define SPH_H

#include "stable.h"
#include <stddef.h>
#include <array>
#include "kernels.h"

#include "neighbourdata.h"

class Scene2;

class SPH {
public:
    const Sceneparser::Sph<number>& getParameters() const {
        return configuration;
    }

    SPH(const Sceneparser::Sph<number>& configuration, Particles* particles, Neighbours *neighbours)
        : configuration(configuration)
        , g(configuration.g)
        , damp(1.0)
        , max_v(0.0)
        , max_c(0.0)
        , particles(particles)
        , neighbours(neighbours)
        , k_simple(14)
        , mod_time(0)
    {
        init();
    }

    number step(const Scene2 &scene);

    void XSPH();

    void move_particles(const Scene2 &scene);

    void handle_stirrer(long long i);
    void handle_stirrer2(long long i);

    number continuity();

    void momentum();

    void shepard();

    number shepard_for_particle_i(long long i);
    void moving_least_squares();
    inline number compute_W_MLS(long long i, const NeighbourData &nd, bool &MLS_failed);
    inline void compute_beta(std::array<number,3> &beta, long long i, bool &MLS_failed);
    inline void  compute_A(std::array<number,6> &A, long long i);
    inline void compute_A_tilde(std::array<number,6> &A_tilde, Vector2 x_i, Vector2 x_j);

    void update_kernels(const Kernel *kernel);

    void update_kernels_f_only(const Kernel *kernel);

    void update_kernels_df_only(const Kernel *kernel);

    void extrapolation();

    number eos_complex(number rho, number rho0, number c, number gamma);

    void getVelocityRange(double &min, double &max) const;
    void getDensityRange(double &min, double &max) const;
    void getPressureRange(double &min, double &max) const;

    Sceneparser::Sph<number> getConfiguration();
    void setConfiguration(Sceneparser::Sph<number> conf);

    Vector2 g;
    number damp, max_v, max_c;
    number dt;
private:

    void init();

    void momentum_single_fluid();

    void momentum_multi_fluid();

    void update_wall_particles();

    Sceneparser::Sph<number> configuration;

    Particles * const particles;
    Neighbours * const neighbours;

    number k_simple;
    number mod_time;
    bool has_second_phase;
};

#endif // SPH_H
