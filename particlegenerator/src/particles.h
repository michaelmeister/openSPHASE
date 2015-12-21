#ifndef PARTICLES_H
#define PARTICLES_H

#include "particle2.h"
#include "configuration.h"
#include <functional>

namespace Particlegenerator {

template<typename number>
using pair = std::pair<number,number>;
template<typename number>
using Particles2 = std::vector<Particle2<number>*>;
using Sceneparser::Configuration;
using Sceneparser::LineSegment;

template<typename number>
class ParticleHelper
{
public:
    explicit ParticleHelper(const Configuration<number>& configuration);
    ~ParticleHelper();

    Particles2<number> *getParticles();

    static unsigned int generateParticleID();

private:
    void initFromConfiguration(); /// called from constructor to build up particles from configuration
    void add_inflow_particles();
    void add_boundary_particles();

    void addParticle(Particle2<number>* p);
    void addRigidParticle(const Vector2<number>& pos, const Vector2<number>& vel);
    void addStirrerParticle(const Vector2<number>& pos);
    void popParticle();
    void addParticles(pair<number> pos, pair<number> dim, std::function<void(Particle2<number>&)> modifier);
    void addSquareFluid1(pair<number> pos, pair<number> dim);
    void addSquareFluid2(pair<number> pos,pair<number> dim);
    void addSquareFluid_Reactor1(pair<number> pos, pair<number> dim);
    void addSquareFluid_Reactor2(pair<number> pos,pair<number> dim);


    Configuration<number>   configuration;
    Particles2<number>      particles;

    bool    inflow;
    number  sample_dist, h;
    number  water_level;
    bool    hydro_pressure;

    std::function<void(Vector2<number>&)> align_to_grid;
};

}

#endif // PARTICLES_H
