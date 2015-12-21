#ifndef CONFIGURATION_H
#define CONFIGURATION_H

#include <vector>
#include "vector2.h"
#include "linesegment.h"

namespace Sceneparser {

template <typename number>
struct AreaStruct
{
    AreaStruct() : x(0.0), y(0.0), width(0.0), height(0.0) {}
    AreaStruct(number x, number y, number width, number height)
        : x(x), y(y), width(width), height(height) {}
    ~AreaStruct() { }

    number x;
    number y;
    number width;
    number height;
};

template <typename number>
struct Grid
{
    Grid() : origin_x(0.0), origin_y(0.0), size(25) { }
    ~Grid() { }

    number origin_x, origin_y;
    number neighbours;
    int size;
};

template <typename number>
struct Reactor
{
    Reactor() : oxygen(0.0), cstr(false), cstr_volume(0.0), inflow(0.0) {}
    ~Reactor() { }

    number oxygen;
    bool   cstr;
    number cstr_volume;
    number inflow;

    struct AreaStruct<number> area;
};

template <typename number>
struct Sph
{
    Sph();
    ~Sph() { }

    struct Stirrer {
        Stirrer();
        ~Stirrer() { }

        bool is_active;
        number origin_x, origin_y, velocity, length_x, length_y;
        enum Type {ellipsoidal_rotator=1, rectangular_fan=2};
        Type type;
    };

    number volume;
    number rho0_1;
    number rho0_2;
    number nu1;
    number nu2;
    number sd;
    number c1;
    number c2;
    number gamma1;
    number gamma2;
    number alpha;
    number epsilon_xsph;
    number epsilon_repulsion;
    number D_air;
    number D_water;
    number originx;
    number originy;
    bool no_slip;
    bool shepard;
    bool moving_least_squares;
    struct Stirrer stirrer, stirrer2;
    Vector2<number> g;
    number t_damp;
    bool hydrostatic_initialization;
};

template <typename number>
struct Asm
{
    Asm();
    ~Asm() { }

    number sd;
    number width;
    number height;
    number ya;
    number yh;
    number fp;
    number ixb;
    number ixp;
    number ixe;
    number mh;
    number Ks;
    number Koh;
    number Kno;
    number bh;
    number ba;
    number eg;
    number eh;
    number kh;
    number rho7;
    number Kx;
    number ma;
    number Knh;
    number Koa;
    number ka;
    number SRT;
    number time_scaling;

    struct Grid<number> grid;
    struct AreaStruct<number> inflow_area, outflow_area;
    struct Reactor<number> reactor1;
    struct Reactor<number> reactor2;
};

template <typename number>
struct Scene
{
    Scene() : sample_dist(0.0), width(0.0), height(0.0), neighbours(0.0), originx(0.0), originy(0.0), hydrostatic_pressure(false), connect_boundaries(false) { }
    ~Scene() { }

    number		sample_dist;
    number		width;
    number		height;
    number		neighbours;
    number      originx;
    number      originy;

    bool        hydrostatic_pressure;
    bool        connect_boundaries;
};

template <typename number>
struct Boundary
{
    Boundary() : is_moving(false) { }
    ~Boundary() { }

    bool is_moving;
    Vector2<number> velocity;
    std::vector<Vector2<number> >   lines;
};

template <typename number>
class Configuration
{
public:
    Configuration();
    ~Configuration();

    Scene<number> scene;
    Sph<number> sph;
    Asm<number> sm;

    std::string json_file;
    std::string import_file;
    std::string export_file;
    std::string inflow_file;
    bool hydrostatic_pressure;

    number rec_step;
    number max_time;

    std::vector<AreaStruct<number> > fluid1;
    std::vector<AreaStruct<number> > fluid2;
    std::vector<AreaStruct<number> > fluid_reactor1;
    std::vector<AreaStruct<number> > fluid_reactor2;
    std::vector<Boundary<number> > boundaries;

    std::vector<LineSegment<number> > solid_boundaries;
    std::vector<LineSegment<number> > periodic_in, periodic_out, rigid_periodic_in, rigid_periodic_out;
    std::vector<LineSegment<number> > horizontal_in, bottom_in, open_out;
    std::vector<LineSegment<number> > vprofile;
    std::vector<LineSegment<number> > qcounter;
    std::vector<LineSegment<number> > hcounter;
};

}

#endif // CONFIGURATION_H
