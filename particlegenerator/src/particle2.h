#ifndef PARTICLE2_H
#define PARTICLE2_H

#include <queue>
#include "vector2.h"

namespace Particlegenerator {

using Sceneparser::Vector2;

enum ParticleType {
    // implicit type nothing <= not fluid not wall
    // Normal  = 1,     // no used
    Out     = 2,        // particle is out and will be removed
    Wall    = 4,        // wall particle
    // OutFlow = 8,     // particle is in outflow buffer zone
    InFlow = 16,        // different from fluid, no fluid calculations performed
    InFlowNew = 32,     // different from fluid, no fluid calculations performed
    Fluid = 64,         // fluid particle, phase 1
    MovingWall = 128,   // sub type of wall, moved with wall velocity
    SecondPhase = 256,  // fluid particle, phase 2
    Rigid = 512,        // sub type of wall, moved with velocity, mirrored when exiting scene (used for aeration)
    Reactor1 = 1024,    // only for ASM (marker), can be calculated with ASM calculation depending on position
    Reactor2 = 2048,    // only for ASM (marker), can be calculated with ASM calculation depending on position
    Stirrer = 4096      // sub type of wall, moving wall with velocity v(t,r)
};

inline void SetFlag(unsigned int& target, ParticleType flag){
    target = target | flag;
}

inline void UnsetFlag(unsigned int& target, ParticleType flag){
    target = target & ~flag;
}

template<typename number>
class Particle2 {
public:
    Particle2(unsigned int id, Vector2<number> position, bool is_wall = false, bool is_fluid = false, bool is_moving_wall = false)
        : id(id), pos(position), type(0x0) {
        if (is_wall) {
            set(Wall);
        }
        if (is_fluid) {
            set(Fluid);
        }
        if (is_moving_wall) {
            set(MovingWall);
        }
        v = v_avg = f = Vector2<number>();
        c = p = rho = mass = t_offset = 0.0;

        selected = false;
    }

    unsigned int id;
    Vector2<number> pos;
    Vector2<number> v;
    Vector2<number> v_avg;
    Vector2<number> v_wall;
    Vector2<number> f;
    number p;
    number rho, volume;
    number c;
    number mass;
    number t_offset;
    unsigned int type;

    bool selected;

    inline bool is(ParticleType t) const {
        return (type & t) != 0;
    }

    inline bool is_out() const {
        return is(Out);
    }

    inline bool is_wall() const {
        return is(Wall);
    }

    inline bool is_inflow() const {
        return is(InFlow);
    }

    inline bool is_inflownew() const {
        return is(InFlowNew);
    }

    inline bool is_fluid() const {
        return is(Fluid);
    }

    inline bool is_moving_wall() const {
        return is(MovingWall);
    }

    inline bool is_second_phase() const {
        return is(SecondPhase);
    }

    inline bool is_rigid() const {
        return is(Rigid);
    }

    inline bool is_reactor1() const {
        return is(Reactor1);
    }

    inline bool is_reactor2() const {
        return is(Reactor2);
    }

    inline bool is_stirrer() const {
        return is(Stirrer);
    }

    inline void set(ParticleType t) {
        type = type | t;
    }

    inline void unset(ParticleType t) {
        type = type & ~t;
    }

    inline void set_selected() {
        selected = true;
    }

    inline void unset_selected() {
        selected = false;
    }
};

}

#endif // PARTICLE2_H
