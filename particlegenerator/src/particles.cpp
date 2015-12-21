#include "particles.h"

#include <algorithm>
#include <cassert>

#ifdef _MSC_VER
#define NOMINMAX
#include <Windows.h>
#endif

using namespace Sceneparser;

namespace Particlegenerator {

template<typename number>
ParticleHelper<number>::ParticleHelper(const Configuration<number>& configuration)
    : configuration(configuration)
    , sample_dist(configuration.scene.sample_dist)
    , h(configuration.scene.neighbours * configuration.scene.sample_dist)
    , water_level(0.f)
    , hydro_pressure(configuration.scene.hydrostatic_pressure)
{
    const double sdhalf = double(sample_dist) / 2.0;
    align_to_grid = [=](Vector2<number>& v)
    {
        v.x = round(double(v.x) / sdhalf) * sdhalf;
        v.y = round(double(v.y) / sdhalf) * sdhalf;
    };
}

template<typename number>
ParticleHelper<number>::~ParticleHelper()
{
}

template<typename number>
Particles2<number>* ParticleHelper<number>::getParticles()
{
    initFromConfiguration();
    Particles2<number>* ret = new Particles2<number>();
    std::swap(*ret, particles);
    return ret;
}

template<typename number>
unsigned int ParticleHelper<number>::generateParticleID()
{
    static unsigned int id = 0;
#ifdef __GNUC__
    return __sync_fetch_and_add(&id, 1);
#elif _MSC_VER
    return InterlockedIncrement(&id);
#endif
}

template<typename number>
void ParticleHelper<number>::initFromConfiguration()
{
    std::for_each(configuration.fluid1.cbegin(), configuration.fluid1.cend(), [this] (const AreaStruct<number>& ar) {
        addSquareFluid1(pair<number>(ar.x, ar.y), pair<number>(ar.width, ar.height));
    });
    std::for_each(configuration.fluid2.cbegin(), configuration.fluid2.cend(), [this] (const AreaStruct<number>& ar) {
        addSquareFluid2(pair<number>(ar.x, ar.y), pair<number>(ar.width, ar.height));
    });
    std::for_each(configuration.fluid_reactor1.cbegin(), configuration.fluid_reactor1.cend(), [this] (const AreaStruct<number>& ar) {
        addSquareFluid_Reactor1(pair<number>(ar.x, ar.y), pair<number>(ar.width, ar.height));
    });
    std::for_each(configuration.fluid_reactor2.cbegin(), configuration.fluid_reactor2.cend(), [this] (const AreaStruct<number>& ar) {
        addSquareFluid_Reactor2(pair<number>(ar.x, ar.y), pair<number>(ar.width, ar.height));
    });
    // TODO check rigid particles and stirrer particles
    //    std::for_each(configuration.boundaries.cbegin(), configuration.boundaries.cend(), [this] (const Boundary<number>& b) {
    //        std::unique_ptr<Vector2<number>> previous(nullptr);
    //        std::for_each(b.lines.cbegin(), b.lines.cend(), [this, &previous, &b] (const Vector2<number>& current) {
    //            if (previous) {
    //                addSolidBoundary(LineSegment<number>(*previous, current, b.velocity, b.is_moving));
    //            }
    //            previous.reset(new Vector2<number>(current));
    //        });
    //    });
    add_boundary_particles();
    add_inflow_particles();
}


template<typename number>
void ParticleHelper<number>::addParticle(Particle2<number>* p) {
    assert(p->pos.x >= 0 - configuration.sph.originx && p->pos.x < configuration.scene.width - configuration.sph.originx);
    assert(p->pos.y >= 0 - configuration.sph.originy && p->pos.y < configuration.scene.height - configuration.sph.originy);
    particles.push_back(p);
}

template<typename number>
void ParticleHelper<number>::addRigidParticle(const Vector2<number>& pos, const Vector2<number>& vel) {
    Particle2<number> *p_i = new Particle2<number>(generateParticleID(), pos, true, false, true);
    p_i->v_wall = vel;
    p_i->set(Rigid);
    addParticle(p_i);
}

template<typename number>
void ParticleHelper<number>::addStirrerParticle(const Vector2<number>& pos) {
    Particle2<number> *p_i = new Particle2<number>(generateParticleID(), pos, true, false, false);
    p_i->set(Stirrer);
    addParticle(p_i);
}

template<typename number>
void ParticleHelper<number>::popParticle() {
    particles.pop_back();
}

template<typename number>
void ParticleHelper<number>::addParticles(pair<number> pos, pair<number> dim, std::function<void(Particle2<number>&)> modifier) {
    const size_t xsamples = lround(abs(dim.first / sample_dist));
    const size_t ysamples = lround(abs(dim.second / sample_dist));
    number sd2 = sample_dist / static_cast<number>(2.0);
    const number xsgn = (dim.first < 0 ? -1.0 : 1.0);
    const number ysgn = (dim.second < 0 ? -1.0 : 1.0);


    for (size_t yit = 0; yit < ysamples; ++yit) {
        for (size_t xit = 0; xit < xsamples; ++xit) {
            Vector2<number> p(pos.first + xsgn * (xit * sample_dist + sd2)
                , pos.second + ysgn * (yit * sample_dist + sd2));
            align_to_grid(p);
            Particle2<number> *p_i = new Particle2<number>(generateParticleID(), p, false, false);
            modifier(*p_i);
            addParticle(p_i);
        }
    }
}

template<typename number>
void ParticleHelper<number>::addSquareFluid1(pair<number> pos, pair<number> dim) {
    addParticles(pos, dim, [](Particle2<number>& p)
    {
        p.set(Fluid);
    });
}

template<typename number>
void ParticleHelper<number>::addSquareFluid2(pair<number> pos, pair<number> dim) {
    addParticles(pos, dim, [](Particle2<number>& p)
    {
        p.set(Fluid);
        p.set(SecondPhase);
    });
}

template<typename number>
void ParticleHelper<number>::addSquareFluid_Reactor1(pair<number> pos, pair<number> dim) {
    addParticles(pos, dim, [](Particle2<number>& p)
    {
        p.set(Fluid);
        p.set(Reactor1);
    });
}

template<typename number>
void ParticleHelper<number>::addSquareFluid_Reactor2(pair<number> pos, pair<number> dim) {
    addParticles(pos, dim, [](Particle2<number>& p)
    {
        p.set(Fluid);
        p.set(Reactor2);
    });
}

template<typename number>
void ParticleHelper<number>::add_inflow_particles() {

    const number rho0_1 = configuration.sph.rho0_1;
    const number rho0_2 = configuration.sph.rho0_2;
    const number c1 = configuration.sph.c1;
    const number c2 = configuration.sph.c2;
    const number rho_part1 = configuration.sph.gamma1 / (rho0_1*c1*c1);
    const number rho_part2 = configuration.sph.gamma2 / (rho0_2*c2*c2);

    for (size_t b = 0; b < configuration.horizontal_in.size(); b++) {
        LineSegment<number> &in = configuration.horizontal_in[b];
        int layer = -1;
        const Vector2<number> in_m_normalized = in.m.normalized();

        while ((++layer)*sample_dist <= h)   {
            Vector2<number> dxv = in_m_normalized*sample_dist;
            Vector2<number> base = in.b + (in.wall_velocity.normalized()*sample_dist*layer);

            while ((base - in.b)*in_m_normalized < (in.m.length() + dxv.length() / 2)) {
                Particle2<number> *p = new Particle2<number>(generateParticleID(), base, false, false);

                if (!in.is_moving) {
                    p->p = hydro_pressure*configuration.sph.rho0_1*abs((configuration.sph.g*in.m.normalized())*((p->pos - in.b)*in.m.normalized()));
                    p->rho = rho0_1*std::pow(p->p*rho_part1 + 1, 1 / configuration.sph.gamma1);
                }
                else {
                    p->p = hydro_pressure*configuration.sph.rho0_2*abs((configuration.sph.g*in_m_normalized)*((p->pos - in.b)*in_m_normalized));
                    p->rho = rho0_2*std::pow(p->p*rho_part2 + 1, 1 / configuration.sph.gamma2);
                    p->set(SecondPhase);
                }
                p->v = in.wall_velocity;
                if (layer == 0)
                    p->set(InFlowNew);
                else
                    p->set(InFlow);
                particles.push_back(p);
                base += dxv;
            }
        }
    }

    for (size_t b = 0; b < configuration.bottom_in.size(); b++) {
        LineSegment<number> &in = configuration.bottom_in[b];
        int layer = -1;
        const Vector2<number> in_m_normalized = in.m.normalized();

        while ((++layer)*sample_dist <= h) {
            Vector2<number> dxv = in_m_normalized*sample_dist;
            Vector2<number> base = in.b + (in.wall_velocity.normalized()*sample_dist*layer);

            while ((base - in.b)*in_m_normalized < (in.m.length() + dxv.length() / 2)) {
                Particle2<number> *p = new Particle2<number>(generateParticleID(), base, false, false);
                p->p = 0;
                if (in.is_moving) {
                    p->rho = rho0_2;
                    p->set(SecondPhase);
                }
                else {
                    p->rho = rho0_1;
                }

                if (layer == 0)
                    p->set(InFlowNew);
                else
                    p->set(InFlow);
                particles.push_back(p);
                base += dxv;
            }
        }
        for (size_t i = 0; i < particles.size(); i++) {
            Particle2<number> *p_i = particles.at(i);

            if (!(p_i->is_wall()))
                water_level = std::max(water_level, in.normal_distance(p_i->pos).length());
        }
    }
}

template<typename number>
void ParticleHelper<number>::add_boundary_particles() {
    const int wall_width = floor(abs(h / sample_dist)) + 1;


    for (size_t i = 0; i < configuration.solid_boundaries.size(); i++) {
        LineSegment<number> &ls = configuration.solid_boundaries[i];

        // check if in continuous line
        bool is_line_strip = false;
        if (i > 0)
        {
            auto&& prev = configuration.solid_boundaries[i - 1];
            Vector2<number> diff(prev.b + prev.m - ls.b);
            if (configuration.scene.connect_boundaries && diff.length() < (sample_dist/10.0))
                is_line_strip = true;
        }

        for (size_t wall_depth = 0; wall_depth < wall_width; ++wall_depth)
        {
            int samples = round(ls.m.length() / sample_dist);
            Vector2<number> dxv = ls.m.normalized()*sample_dist; // sampling step size
            Vector2<number> base = ls.b;
            base += ls.normal()*sample_dist*(number(wall_depth) + 0.5); // into wall depth dir
            base += dxv / 2; // into line dir
            if (is_line_strip)
            {
                base -= dxv * wall_width;
                samples += wall_width;
            }

            for (int s = 0; s < samples; ++s)
            {
                align_to_grid(base);
                Particle2<number> *p = new Particle2<number>(generateParticleID(), base, true, false, ls.is_moving);
                p->v_wall = ls.wall_velocity;
                particles.push_back(p);
                base += dxv;
            }
        }
    }
}

template class ParticleHelper<float>;
template class ParticleHelper<double>;

}
