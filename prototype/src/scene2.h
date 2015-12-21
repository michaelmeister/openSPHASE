#ifndef SCENE2_H
#define SCENE2_H

#include "stable.h"

#include <memory>
#include <vector>
#include <set>
#include <boost/foreach.hpp>
#include <limits>
#include <iostream>
using namespace std;

#include "unigrid.h"
#include "sph.h"
#include "neighbourdata.h"
#include "kernels.h"

#include <string>
#include <sstream>
#include <assert.h>
#include <map>

#include <boost/bind.hpp>
#include <thrust/remove.h>
#include <stdio.h>

#include <configuration.h>
#include <jsonparser.h>
#include <particle2.h>
#include <configurationparser.h>

struct H5PartFile;

class Scene2 {
    friend class ParticleRenderer;
    friend class JSONParser;
    friend class MainWindow;
public:

    enum class eCalculationMode {
        eCalcSPH,
        eVisualize
    };

    Scene2(eCalculationMode mode, const Sceneparser::Configuration<number>& configuration, Particles* particles, bool single_file);
    ~Scene2();

    void check_overlapping_particles();

    inline number calc_damping_factor(number time);
    bool iterate(size_t step, number& time);// returns false on break up
    bool step(size_t step, number& time);	// returns false on break up
    Particles *getParticles() const;

    void addSPHImportFile(std::string part5importFile);
    void addASMImportFile(std::string part5importFile);
    void addSPHExportFile(std::string part5exportFile);
    void addInflowFile(std::string inflowFile);

    void dumpToH5PartFile(number time);
    size_t getNumImportedSteps();
    void setFrameDelay(double delay);
    double getFrameDelay();
    Sceneparser::Sph<number> getSPHParameters();
    void setParameters(Sceneparser::Configuration<number> conf);

    bool isDumpASM();
    void setDumpASM(bool value);

    bool ImportStepFromH5Part(size_t step, number& time);
    void ExportStepToH5Part(number time);

private:
    /**
      * calculates the mass for the particles according to rho0 and the
      * particle distribution. also sets the  influence radius according to
      * particles volumes.
      */
    void init();

    bool stepSph(size_t step, number qd_factor, number& time);
    bool stepVisualize(size_t step, number& time);

    void addParticle(Particle2*);

    void ExportConfigurationToH5Part(H5PartFile* part_file);
    void ImportConfigurationFromH5Part(H5PartFile* part_file);
    void ExportStepToCSV(number time);
    bool ImportH5PartSph(size_t step, number& time);
    void ExportSphStepToH5Part(number time);

    /**
      * precache neighbour particles.
      */
    inline void cache_neighbours_unigrid();
    inline void cache_neighbours_bruteforce();
    inline void cache_neighbours();

    /**
      * shifts particles that touch the out boundary to the in boundary (periodic channel)
      */
    inline void handle_periodic_boundaries();

    /**
        * determines type of inflow particles and handles inflow with constant velocity
        */
    inline void handle_inflow(number time, number qd_factor);
    inline void handle_horizontal_inflow(number qd_factor);
    inline void handle_bottom_inflow(number time);
    inline void handle_surface_degassing(Neighbours neighbours);

    /**
        * removes particles that cross the (open) outflow threshold from the scene
        */
    inline void handle_open_outflow();
    void handle_velocity_profile();
    void handle_qcounter();
    inline void handle_hcounter();
    inline void handle_kineticenergycounter();

    const eCalculationMode mode;
    Sceneparser::Configuration<number> configuration;

    number      flow1;
    number      flow2;
    Particles* const particles;
    Neighbours	neighbours;
    number		h;//cached values;
    number      water_level;
    bool		measure, bottom_in_state, active_inflow, inflow;
    SPH			sph;

    UniGrid*	grid;
    std::string sph_part_file_name;
    std::string asm_part_file_name;
    std::string csv_file_name;
    size_t		numImportedSteps;
    double		frameDelay;
    int         tmp_step;
    std::map<number, std::vector<number>> inflow_map;
    number      dump_time;

    std::string sph_id;

    bool        has_asm_import;
    bool        has_sph_import;

    int         record_step;

    bool        single_file;

public:
    number      time;
    const Kernel *kernel;

    number max_time() { return configuration.max_time; }
    number rec_step() { return configuration.rec_step; }

    bool outOfScene(const Vector2& pos) const;
    void eraseParticles();
    bool isInflowActive() const;
    bool hasInflow() const;
    bool hasImport() const;
    std::string getSphPartFileName() const;
    std::string getAsmPartFileName() const;
    std::string getCsvFileName() const;
    void exportNumSteps();
};


#endif // SCENE2_H
