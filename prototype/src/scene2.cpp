#include "scene2.h"

#define _USE_MATH_DEFINES

#include <H5Part.h>
#include <boost/algorithm/string/replace.hpp>
#include <boost/uuid/uuid.hpp>
#include <boost/uuid/uuid_generators.hpp>
#include <boost/uuid/uuid_io.hpp>
#include <fstream>
#include <cmath>
#include "math.h"
#include <memory>

#include <configuration.h>


inline bool is_out(const Particle2 *const p) {
    if (p->is_out()) {
        //delete p;
        return true;
    }
    return false;
}

Scene2::Scene2(eCalculationMode mode, const Sceneparser::Configuration<number>& conf, Particles *particles, bool single_file)
    : mode(mode)
    , configuration(conf)
    , h(conf.scene.neighbours*conf.sph.sd)
    , particles(particles)
    , sph(SPH(conf.sph, this->particles, &(this->neighbours)))
    , dump_time(0)
    , single_file(single_file)
{
    bottom_in_state = false;
    measure = false;
    numImportedSteps = 0;
    frameDelay = 0;
    inflow = 0;
    active_inflow = 0;
    tmp_step = 0;
    has_asm_import = false;
    has_sph_import = false;
    record_step = 0;
    inflow_map.clear();

    if (!conf.inflow_file.empty())
        addInflowFile(configuration.inflow_file);

    if (mode ==  eCalculationMode::eCalcSPH) {
        kernel = Kernel::get("wendland2", h);

        // if we don't use xsph, the unigrid doesn't have to initialize f & df from any kernel, we update them anyway at sph.step()
        if (configuration.sph.epsilon_xsph == 0)
            grid = new UniGrid(number_pair(configuration.scene.width, configuration.scene.height), NULL, h);
        else
            grid = new UniGrid(number_pair(configuration.scene.width, configuration.scene.height), kernel, h);

        addSPHExportFile(configuration.export_file);

        init();
    }
}

Scene2::~Scene2()
{
    std::for_each(particles->begin(), particles->end(), [] (Particle2* p) {
        delete p;
    });
    delete particles;
    if (grid)
        delete grid;
}

void Scene2::addParticle(Particle2 *p)
{
    assert(p->pos.x >= 0 - sph.getParameters().originx && p->pos.x <  configuration.scene.width - sph.getParameters().originx);
    assert(p->pos.y >= 0 - sph.getParameters().originy && p->pos.y < configuration.scene.height - sph.getParameters().originy);
    particles->push_back(p);
}

void Scene2::check_overlapping_particles() {
    number sample_dist_limit2 = 0.03*configuration.sph.sd*configuration.sph.sd;
    for (int i = 0; i < particles->size(); i++) {
        Particle2 *p_i = particles->at(i);
        BOOST_FOREACH(const NeighbourData &nd, neighbours.at(i)) {
            Particle2 *p_j = particles->at(nd.j);
            if (((p_i->pos - p_j->pos).length_squared() < sample_dist_limit2) && (i != nd.j)) {
                // mark particle out if fluid / stirrer conflict
                if (p_i->is_fluid() && p_j->is_stirrer()) {
                    cout << "resolved fluid / stirrer conflict of overlapping particles [" << i << ',' << nd.j << "] \n" << endl;
                    p_i->set(Particlegenerator::Out);
                } else if (!(p_j->is_fluid() && p_i->is_stirrer())) {
                    cout << "overlapping particles [" << i << ',' << nd.j << "] \n" << endl;
                }

            }
        }
    }
}

/**
  * calculates the mass for the particles according to rho0 and the
  * particle distribution. also sets the  influence radius according to
  * particles volumes.
  */
void Scene2::init() {
    assert(configuration.periodic_in.size() == configuration.periodic_out.size());
    cout << "initing " << particles->size() << " particles" << endl;
    sph.max_v = 0;
    sph.damp = 0;
    measure = false;
    cache_neighbours();
    check_overlapping_particles();
    bottom_in_state = false;
}

Sceneparser::Sph<number> Scene2::getSPHParameters() {
    return sph.getConfiguration();
}

void Scene2::setParameters(Sceneparser::Configuration<number> conf) {
    configuration = conf;
    sph.setConfiguration(conf.sph);
}

inline number Scene2::calc_damping_factor(number time) {
    return 0.5*(sin((-0.5 + (time / configuration.sph.t_damp))*M_PI) + 1.0);
}

void Scene2::ExportConfigurationToH5Part(H5PartFile* part_file) {
    Sceneparser::ConfigurationParser<number> parser(configuration);
    std::string json = parser.getJson();
    int json_size = json.size();

    if (H5PartWriteFileAttribString(part_file, "json_blob", json.c_str()) < 0 || H5PartWriteFileAttribInt64(part_file, "json_blob_size", json_size) < 0)
        throw std::runtime_error("Could not write configuration to H5Part file!");
}

void Scene2::ImportConfigurationFromH5Part(H5PartFile* part_file) {
    h5part_int64_t json_size = -1;
    if (H5PartReadFileAttrib(part_file, "json_blob_size", &json_size) < 0 || json_size == -1)
        throw std::runtime_error("Could not read configuration from H5Part file!");

    char* json = new char[json_size+1];
    H5PartReadFileAttrib(part_file, "json_blob", json);
    json[json_size] = '\0';

    std::string fileName = "tmpConfigurationFile.json";

    std::ofstream outfile;
    outfile.open(fileName, ios::out);
    outfile.write(json, sizeof(char)*json_size);
    outfile.close();

    Sceneparser::JSONParser<number> parser(fileName, "", "", 0, 0, "");
    setParameters(parser.getConfiguration());

    remove(fileName.c_str());

    delete[] json;
}

void Scene2::ExportStepToH5Part(number time)
{
    switch (this->mode) {
    case eCalculationMode::eCalcSPH:
        ExportSphStepToH5Part(time);
        break;
    case eCalculationMode::eVisualize:
        // no-op
        break;
    default:
        break;
    }
}

void Scene2::ExportSphStepToH5Part(number time)
{
    // only export if rec step condition holds
    if (time < dump_time) {
        return;
    }
    dump_time = time + rec_step();

    static std::string filename = getSphPartFileName();
    if (!single_file && record_step >= 1000 && record_step%1000 == 0) {
        filename = getSphPartFileName().substr(0, getSphPartFileName().find_first_of('.')) + "_" + std::to_string(record_step/1000) + "." + getSphPartFileName().substr(getSphPartFileName().find_first_of('.') + 1);
    }

    H5PartFile* part_file = H5PartOpenFile(filename.c_str(), H5PART_APPEND);
    if (!part_file) {
        throw std::runtime_error("Could not open SPH export file!");
    }

    H5PartSetStep(part_file, record_step++);

    size_t npart = getParticles()->size();
    H5PartSetNumParticles(part_file, npart);

    std::unique_ptr<unsigned int[]> id(new unsigned int[npart]);
    std::unique_ptr<double[]> x(new double[npart]);
    std::unique_ptr<double[]> y(new double[npart]);
    std::unique_ptr<double[]> z(new double[npart]);
    std::unique_ptr<double[]> vx(new double[npart]);
    std::unique_ptr<double[]> vy(new double[npart]);
    std::unique_ptr<double[]> vz(new double[npart]);
    std::unique_ptr<double[]> density(new double[npart]);
    std::unique_ptr<double[]> pressure(new double[npart]);
    std::unique_ptr<unsigned int[]> type(new unsigned int[npart]);

    for (int i = 0; i < npart; i++) {
        Particle2 *p = getParticles()->at(i);
        id[i] = p->id;
        x[i] = p->pos.x;
        y[i] = p->pos.y;
        z[i] = 0.0;

        density[i] = p->rho;
        pressure[i] = p->p;
        vx[i] = p->v.x;
        vy[i] = p->v.y;
        vz[i] = 0.0;
        type[i] = p->type;
    }

    H5PartWriteStepAttribFloat64(part_file, "time", time);
    H5PartWriteDataInt32(part_file, "id", (int*)id.get());
    H5PartWriteDataFloat64(part_file, "x", x.get());
    H5PartWriteDataFloat64(part_file, "y", y.get());
    H5PartWriteDataFloat64(part_file, "z", z.get());
    H5PartWriteDataFloat64(part_file, "vx", vx.get());
    H5PartWriteDataFloat64(part_file, "vy", vy.get());
    H5PartWriteDataFloat64(part_file, "vz", vz.get());
    H5PartWriteDataFloat64(part_file, "density", density.get());
    H5PartWriteDataFloat64(part_file, "pressure", pressure.get());
    H5PartWriteDataInt32(part_file, "type", (int*)type.get());
    //	H5Block3dWrite3dVectorFieldFloat64(part_file, "vel", vx, vy, vz);

    H5PartCloseFile(part_file);
}

bool Scene2::ImportStepFromH5Part(size_t step, number& time)
{
    switch (this->mode) {
    case eCalculationMode::eCalcSPH:
        // no-op
        return true;
        break;
    case eCalculationMode::eVisualize:
        return  ImportH5PartSph(step, time);
        break;
    default:
        break;
    }
}

bool Scene2::ImportH5PartSph(size_t step, number& time)
{
    static std::string filename = getSphPartFileName();
    static bool multiple_files = false;

    if (multiple_files)
        if (step >= 1000) {
            filename = getSphPartFileName().substr(0, getSphPartFileName().find_first_of('.')) + "_" + std::to_string(step/1000) + "." + getSphPartFileName().substr(getSphPartFileName().find_first_of('.') + 1);
        } else {
            filename = getSphPartFileName();
        }

    H5PartFile* part_file = H5PartOpenFile(filename.c_str(), H5PART_READ);

    if (H5PartSetStep(part_file, step) != H5PART_SUCCESS) {
        H5PartCloseFile(part_file);
        if (step >= 1000) {
            filename = getSphPartFileName().substr(0, getSphPartFileName().find_first_of('.')) + "_" + std::to_string(step/1000) + "." + getSphPartFileName().substr(getSphPartFileName().find_first_of('.') + 1);
        } else {
            filename = getSphPartFileName();
        }
        part_file = H5PartOpenFile(filename.c_str(), H5PART_READ);
        if (H5PartSetStep(part_file, step) != H5PART_SUCCESS) {
            H5PartCloseFile(part_file);
            return false;
        } else {
            multiple_files = true;
        }
    }

    size_t npart = H5PartGetNumParticles(part_file);

    H5PartReadStepAttrib(part_file, "time", &time);

    std::unique_ptr<unsigned int[]> id(new unsigned int[npart]);
    std::unique_ptr<double[]> x(new double[npart]);
    std::unique_ptr<double[]> y(new double[npart]);
    std::unique_ptr<double[]> z(new double[npart]);
    std::unique_ptr<double[]> vx(new double[npart]);
    std::unique_ptr<double[]> vy(new double[npart]);
    std::unique_ptr<double[]> density(new double[npart]);
    std::unique_ptr<double[]> pressure(new double[npart]);
    std::unique_ptr<unsigned int[]> type(new unsigned int[npart]);

    H5PartReadDataInt32(part_file, "id", (int*)id.get());
    H5PartReadDataFloat64(part_file, "x", x.get());
    H5PartReadDataFloat64(part_file, "y", y.get());
    H5PartReadDataFloat64(part_file, "z", z.get());
    H5PartReadDataFloat64(part_file, "vx", vx.get());
    H5PartReadDataFloat64(part_file, "vy", vy.get());
    H5PartReadDataFloat64(part_file, "density", density.get());
    H5PartReadDataFloat64(part_file, "pressure", pressure.get());
    H5PartReadDataInt32(part_file, "type", (int*)type.get());

    H5PartCloseFile(part_file);

    while (getParticles()->size() < npart)
        addParticle(new Particle2(0, Vector2()));

    while (getParticles()->size() > npart)
        getParticles()->pop_back();

    for (size_t i = 0; i < getParticles()->size(); i++) {
        Particle2 *ps = getParticles()->at(i);

        ps->id = id[i];
        ps->pos.x = x[i];
        ps->pos.y = y[i];
        ps->rho = density[i];
        ps->p = pressure[i];
        ps->v.x = vx[i];
        ps->v.y = vy[i];
        ps->type = type[i];
    }

    return true;
}

void mssleep(long ms)
{
#ifdef Q_OS_WIN
    Sleep(uint(ms));
#else
    struct timespec ts = { ms / 1000, (ms % 1000) * 1000 * 1000 };
    nanosleep(&ts, NULL);
#endif
}

void Scene2::eraseParticles()
{
    particles->erase(std::remove_if(
                         particles->begin(),
                         particles->end(),
                         is_out), particles->end());

//    grid->update(particles);
    cache_neighbours();
}

bool Scene2::isInflowActive() const
{
    return active_inflow;
}

bool Scene2::hasInflow() const
{
    return inflow;
}

void Scene2::exportNumSteps()
{
    if (this->mode == eCalculationMode::eCalcSPH) {
        H5PartFile* part_file = H5PartOpenFile(getSphPartFileName().c_str(), H5PART_APPEND);
        if (part_file) {
            H5PartWriteFileAttribInt64(part_file, "num_steps", record_step);
            H5PartCloseFile(part_file);
        }
    }
}

bool Scene2::iterate(size_t step, number& time)
{
    // Import
    if (false == this->ImportStepFromH5Part(step, time))
        return false;
    // Calculate
    if (false == this->step(step, time))
        return false;
    // Export
    this->ExportStepToH5Part(time);
    return true;
}

bool Scene2::step(size_t step, number& time)
{
    this->time = time;

    number qd_factor = 1;
    number bsb5_factor = 1;
    number tkn_factor = 1;

    switch (this->mode) {
    case eCalculationMode::eCalcSPH:
        return stepSph(step, qd_factor, time);
        break;
    case eCalculationMode::eVisualize:
        return true;
        break;
    default:
        throw std::runtime_error("Mode not supported! Use <sph>, <asm> or <visualize>.");
        break;
    }
}

bool Scene2::stepSph(size_t step, number qd_factor, number& time) {
    // PRE
    handle_inflow(time, qd_factor);

    handle_open_outflow();
    handle_periodic_boundaries();

    eraseParticles();

    if (time < configuration.sph.t_damp)
        sph.damp = calc_damping_factor(time);
    else
        sph.damp = 1.0;

    // CALC
    time += sph.step(*this);

    if (0 == step % 10000) {
        std::cout << "SPH iteration " << step << "\tSimulation Time: " << time << " sec" <<std::endl;
    }

    return true;
}

bool Scene2::stepVisualize(size_t step, number& time) {

    // post processing
    if (measure) {
        handle_hcounter();
        handle_qcounter();
        //handle_kineticenergycounter();
        //handle_velocity_profile();
    }

    return ImportH5PartSph(step, time);
}

/**
  * checks if particle p is out of the scene.
  */
bool Scene2::outOfScene(const Vector2& pos) const {
    return (pos.x + sph.getParameters().originx < 0 || pos.x + sph.getParameters().originx > configuration.scene.width || pos.y + sph.getParameters().originy < 0 || pos.y + sph.getParameters().originy > configuration.scene.height);
}

std::string Scene2::getSphPartFileName() const {
    return sph_part_file_name;
}

std::string Scene2::getAsmPartFileName() const {
    return asm_part_file_name;
}

std::string Scene2::getCsvFileName() const {
    return csv_file_name;
}

Particles *Scene2::getParticles() const {
    return particles;
}

void Scene2::addSPHImportFile(std::string part5importFile)
{
    sph_part_file_name = part5importFile + ".h5part";
    H5PartFile* part_file = H5PartOpenFile(sph_part_file_name.c_str(), H5PART_READ);
    if (part_file)
    {
        cout << "Successfully loaded h5part file for SPH data: " << sph_part_file_name << std::endl;

        if (H5PartReadFileAttrib(part_file, "num_steps", &numImportedSteps) < 0)
            numImportedSteps = 0;
        else
            numImportedSteps -= 1;

        char id[36];
        if (H5PartReadFileAttrib(part_file, "sph_id", &id) < 0)
            cout << "Could not read SPH ID from SPH H5Part file!" << endl;

        sph_id = string(id);
        ImportConfigurationFromH5Part(part_file);
        has_sph_import = true;

        H5PartCloseFile(part_file);

    } else {
        throw std::runtime_error("SPH import file could not be opened!");
    }
}

void Scene2::addASMImportFile(std::string part5importFile)
{
    asm_part_file_name = part5importFile; // + ".h5part";
    H5PartFile* part_file = H5PartOpenFile(asm_part_file_name.c_str(), H5PART_READ);
    if (part_file)
    {
        char tmp_sph_id[36];
        if (H5PartReadFileAttrib(part_file, "sph_id", &tmp_sph_id) < 0)
            cout << "Could not read SPH ID from ASM H5Part file!" << endl;
        if (sph_id.compare(string(tmp_sph_id)) != 0)
            cout << "SPH IDs does not match!" << endl;

        ImportConfigurationFromH5Part(part_file);
        cout << "Successfully loaded h5part file for ASM data" << std::endl;
        has_asm_import = true;

        H5PartCloseFile(part_file);
    } else {
        throw std::runtime_error("ASM import file could not be opened!");
    }
}

void Scene2::addSPHExportFile(std::string part5exportFile)
{
    sph_part_file_name = part5exportFile + ".h5part";
    H5PartFile* part_file = H5PartOpenFile(sph_part_file_name.c_str(), H5PART_WRITE);
    if (part_file) {

        // add unique SPH ID
        boost::uuids::uuid uuid = boost::uuids::random_generator()();
        std::string id = boost::uuids::to_string(uuid);
        if (H5PartWriteFileAttribString(part_file, "sph_id", id.c_str()) < 0)
            throw std::runtime_error("Could not write SPH ID to SPH export file!");

        // write configuration
        ExportConfigurationToH5Part(part_file);

        cout << "Successfully opened h5part file for writing SPH data: " << sph_part_file_name << std::endl;

        H5PartCloseFile(part_file);

    } else {
        throw std::runtime_error("SPH export file could not be opened!");
    }
}

void Scene2::addInflowFile(std::string inflowFile) {
    std::ifstream infile;
    infile.open(inflowFile.c_str());

    std::string line;
    number time = 0.0;
    while(std::getline(infile, line)) {
        std::istringstream s(line);
        std::string field;
        std::vector<number> inflow_vector;
        getline(s, field,';');
        time = ::atof(field.c_str());
        for (int i=0; i<3; i++) {
            getline(s, field,';');
            inflow_vector.push_back(::atof(field.c_str()));
        }
        inflow_map[time] = inflow_vector;
    }
}

size_t Scene2::getNumImportedSteps()
{
    return numImportedSteps;
}

void Scene2::setFrameDelay(double delay)
{
    frameDelay = delay;
}

double Scene2::getFrameDelay()
{
    return frameDelay;
}

/**
  * precache neighbour particles.
  */
inline void Scene2::cache_neighbours_unigrid() {
    grid->update(particles);
    if (neighbours.size() < particles->size())
        neighbours.resize(particles->size());

    number ds = sph.dt*sph.max_v;
    number distLimit = h + ds*1.001 + 10E-15;

#pragma omp parallel for
    for (long long i = 0; i < (long long)particles->size(); i++) {
        Particle2 *p = particles->at(i);
        std::vector<NeighbourData> &nd = neighbours[i];
        nd.clear();

        if (configuration.periodic_in.empty() || configuration.periodic_out.empty()) {
             grid->neighbours(i, nd, Vector2(0, 0), distLimit);
        }
        else if (configuration.periodic_out[0].dist(p->pos) <= h) {
            grid->neighbours(i, nd, configuration.periodic_out[0].b - configuration.periodic_in[0].b, distLimit);
        }
        else if (configuration.periodic_in[0].dist(p->pos) <= h) {
             grid->neighbours(i, nd, configuration.periodic_in[0].b - configuration.periodic_out[0].b, distLimit);
        }
        else {
            grid->neighbours(i, nd, Vector2(0, 0), distLimit);
        }
    }
}

inline void Scene2::cache_neighbours_bruteforce() {
    if (neighbours.size() < particles->size()) {
        neighbours.resize(particles->size());
    }

#pragma omp parallel for
    for (long long i = 0; i < (long long)particles->size(); i++) {
        neighbours[i].clear();
        Vector2 pos_i = particles->at(i)->pos;
        for (int j = 0; j < particles->size(); ++j) {
            Vector2 pos_j = particles->at(j)->pos;
            Vector2 dir = pos_i - pos_j;
            number dist = dir.length();
            if (dist <= h) {
                neighbours[i].push_back(NeighbourData(dir, dist, j, kernel, Vector2(0, 0)));
            }
        }
    }
}

inline void Scene2::cache_neighbours() {
    cache_neighbours_unigrid();
}

/**
  * shifts particles that touch the out boundary to the in boundary (periodic channel)
  */
inline void Scene2::handle_periodic_boundaries() {
    for (size_t b = 0; b < configuration.periodic_in.size(); b++) {
        const LineSegment &out = configuration.periodic_out[b];
        const LineSegment &in = configuration.periodic_in[b];
        for (long long i = 0; i < (long long)particles->size(); i++) {
            Particle2 *p = particles->at(i);
            if (p->is_fluid() || p->is_second_phase()) {
                if (out.cross(p->pos)) {
                    p->pos = in.b + (out.b-p->pos).length()*in.m.normalized();
                    p->v = p->v.length()*in.normal();
                }
            } else
                continue;
        }
    }

    for (size_t b = 0; b < configuration.rigid_periodic_in.size(); b++) {
        const LineSegment &out = configuration.rigid_periodic_out[b];
        const LineSegment &in = configuration.rigid_periodic_in[b];
        for (long long i = 0; i < (long long)particles->size(); i++) {
            Particle2 *p = particles->at(i);
            if (p->is_moving_wall() || p->is_rigid()) {
                if (out.cross(p->pos)) {
                    p->pos = in.b + (out.b-p->pos).length()*in.m.normalized();
                    p->v = p->v.length()*in.normal();
                }
            } else
                continue;
        }
    }
}

/**
    * determines type of inflow particles and handles inflow with constant velocity
    */
inline void Scene2::handle_inflow(number time, number qd_factor) {
    if (!configuration.horizontal_in.empty())
        handle_horizontal_inflow(qd_factor);

    if (!configuration.bottom_in.empty())
        handle_bottom_inflow(time);
}

inline void Scene2::handle_horizontal_inflow(number qd_factor) {
    for (long long i = 0; i < (long long)particles->size(); i++) {
        for (size_t b = 0; b <  configuration.horizontal_in.size(); b++) {
            const LineSegment &in =  configuration.horizontal_in[b];
            Particle2 *p = particles->at(i);
            if (!p->is_inflow() && (!p->is_inflownew())) continue;
            if (bottom_in_state)
                p->v.clear();
            else
                p->v = in.wall_velocity * qd_factor;
            if ((p->is_inflownew()) && (in.dist(p->pos) >= (configuration.scene.sample_dist))) {
                Particle2 *new_p = new Particle2(*p);
                new_p->pos -= Vector2(in.normal_distance(new_p->pos));
                p->unset(Particlegenerator::InFlowNew);
                p->set(Particlegenerator::InFlow);
                addParticle(new_p);
            }
            else if ((p->is_inflow()) && (in.dist(p->pos) >= (h + configuration.sph.sd))) {
                p->unset(Particlegenerator::InFlow);
                p->set(Particlegenerator::Fluid);
                p->set(Particlegenerator::Reactor1);
            }
        }
    }
}

inline void Scene2::handle_bottom_inflow(number time) {

    const number rho0_1 = sph.getParameters().rho0_1;
    const number current_water_level = water_level;
    number rho0;
    number c0;
    number gamma;
    number rho_part;
    water_level = 0;

    for (size_t b = 0; b < configuration.bottom_in.size(); b++) {
        const LineSegment &in = configuration.bottom_in[b];
        const number projected_g = fabs(sph.g*in.normal());
        if (in.is_moving) {
            rho0 = sph.getParameters().rho0_2;
            c0 = sph.getParameters().c2;
            gamma = sph.getParameters().gamma2;
            rho_part = gamma / (rho0*c0*c0);
        }
        else {
            rho0 = sph.getParameters().rho0_1;
            c0 = sph.getParameters().c1;
            gamma = sph.getParameters().gamma1;
            rho_part = gamma / (rho0*c0*c0);
        }

#pragma omp parallel for
        for (long long i = 0; i < (long long)particles->size(); i++) {
            Particle2 *p = particles->at(i);
            if (p->is_wall()) continue;

            if (fabs((in.b+number(0.5)*in.m-p->pos)*in.m.normalized()) <= (in.m.length()+configuration.sph.sd)/number(2.0))  {
                if (!(p->is_second_phase()))
                    water_level = (std::max)(water_level, in.normal_distance(p->pos).length());
            }
            else continue;

            if (p->is_inflownew()) {
                p->p = sph.damp*sph.damp*rho0_1*projected_g*(current_water_level-in.dist(p->pos));
                p->rho = rho0*std::pow(p->p*rho_part + 1, 1 / gamma);
                if ((time > configuration.sph.t_damp) && (!bottom_in_state))
                    p->v = in.wall_velocity;
                else
                    p->v = Vector2 (0,0);
                if (in.dist(p->pos) >= (configuration.sph.sd)) {
                    Particle2 *new_p = new Particle2(*p);
                    new_p->pos -= Vector2(in.normal_distance(p->pos));
                    p->unset(Particlegenerator::InFlowNew);
                    p->set(Particlegenerator::InFlow);
                    addParticle(new_p);
                }
            }
            else if (p->is_inflow()) {
                p->p = sph.damp*sph.damp*rho0_1*projected_g*(current_water_level-in.dist(p->pos));
                p->rho = rho0*std::pow(p->p*rho_part + 1, 1 / gamma);
                if ((time > configuration.sph.t_damp) && (!bottom_in_state))
                    p->v = in.wall_velocity;
                else
                    p->v = Vector2 (0,0);
                if (in.dist(p->pos) >= (h + configuration.scene.sample_dist)) {
                    p->unset(Particlegenerator::InFlow);
                    p->set(Particlegenerator::Fluid);
                }
            }
        }
    }
    //handle_surface_degassing(water_level, neighbours);
}

inline void Scene2::handle_surface_degassing(Neighbours neighbours) {
    for (size_t b = 0; b < configuration.bottom_in.size(); b++) {
#pragma omp parallel for
        for (long long i = 0; i < (long long)particles->size(); i++) {
            Particle2 *p = particles->at(i);
            if ((p->is_second_phase()) && !((p->is_inflow()) || (p->is_inflownew())) && (neighbours[i].size() < 7))
                p->set(Particlegenerator::Out);
        }
    }
}

/**
    * removes particles that cross the (open) outflow threshold from the scene
    */
inline void Scene2::handle_open_outflow() {
    if (configuration.open_out.empty())
        return;

    for (int i = 0; i < particles->size(); i++) {
        Particle2 *p = particles->at(i);
        if (p->is_wall()) continue;
        for (size_t b = 0; b < configuration.open_out.size(); b++) {
            const LineSegment &out = configuration.open_out[b];
            if (out.dist(p->pos) < (h + configuration.scene.sample_dist)) {
                p->unset(Particlegenerator::Fluid);
                if (!out.cross(p->pos)) continue;
                p->type = p->type | Particlegenerator::Out;
            }
        }
    }
}

void Scene2::handle_velocity_profile() {
    for (size_t b = 0; b < configuration.vprofile.size(); b++) {
        const LineSegment &vp = configuration.vprofile[b];
        for (size_t i = 0; i < particles->size(); i++) {
            Particle2 *p = particles->at(i);
            if (p->is_wall()) continue;

            if (vp.dist(p->pos) < configuration.scene.sample_dist*0.48)
                printf("%f, %f;", abs((p->pos - vp.m - vp.b)*vp.m.normalized()), p->v.length());
        }
        printf("-----------------------------------------------------------------------------\n");
    }
}

void Scene2::handle_qcounter() {
    for (size_t b = 0; b < configuration.qcounter.size(); b++) {
        const LineSegment &qc = configuration.qcounter[b];
        int counted = 0;
        number v = 0;
        number projection_distance = 0;
        number upper_projection = 0;
        number lower_projection = qc.m.length();

        for (int i = 0; i < particles->size(); i++) {
            Particle2 *p = particles->at(i);
            if ((p->is_wall()) || (qc.dist(p->pos) > h)) continue;
            counted++;
            v += p->v.length();
            projection_distance = (qc.b + qc.m - p->pos + qc.normal_distance(p->pos)).length();
            upper_projection = (std::max)(upper_projection, projection_distance);
            lower_projection = (std::min)(lower_projection, projection_distance);
        }
        if (counted) {
            v /= counted;
            cout << "average q for counter " << (upper_projection - lower_projection)*v << "\n";
        }
    }
}

inline void Scene2::handle_hcounter() {
    for (size_t b = 0; b < configuration.hcounter.size(); b++) {
        const LineSegment &hc = configuration.hcounter[b];
        number projected_height = 0;

        for (int i = 0; i < particles->size(); i++) {
            Particle2 *p = particles->at(i);
            if ((p->is_wall()) || (hc.dist(p->pos) > h)) continue;
            projected_height = (std::max)(projected_height, (hc.b + hc.m - p->pos + hc.normal_distance(p->pos)).length());
        }
        cout << "weir head" << projected_height << "\n";
    }
}

inline void Scene2::handle_kineticenergycounter() {
    number dsum = 0;
    for (int i = 0; i < particles->size(); i++) {
        Particle2 *p = particles->at(i);
        if (p->is_wall()) continue;
        dsum += p->rho * p->volume * p->v * p->v;
    }
    cout << "total kinetic energy " << dsum / 2.0 << "\n";
}
