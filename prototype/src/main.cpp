#ifndef NOGUI
#include <QApplication>
#include <QProcessEnvironment>
#include <QGLFormat>
#endif

#include <iostream>

#ifndef NOGUI
#include "mainwindow.h"
#endif

#include "scene2.h"
#include <QFile>

#include <jsonparser.h>
#include <particles.h>

using namespace std;

Scene2 *applyArguments(int argc, char *argv[]) {
    QString import_file = "";
    QString output_file = "";
    QString json_file = "";
    QString inflow_file = "";
    QString mode = "";
    bool single_file = false;
    double max_time = 0.0;
    double rec_step = 0.0;

    while (1) {
        static struct option long_options[] =
        {
            {"json",        required_argument, 0, 'j'},
            {"in",          required_argument, 0, 'i'},
            {"out",         required_argument, 0, 'o'},
            {"rec_step",    required_argument, 0, 'r'},
            {"max_time",    required_argument, 0, 'm'},
            {"inflow",      required_argument, 0, 'f'},
            {"mode",        required_argument, 0, 'x'},
            {"single_file", 0, 0, 's'},
            {0, 0, 0, 0}
        };
        /* getopt_long stores the option index here. */
        int option_index = 0;

        int c = getopt_long (argc, argv, ":x:i:o:r:m:j:fs", long_options, &option_index);

        /* Detect the end of the options. */
        if (c == -1)
            break;

        switch (c) {
        case 'j':
            json_file = optarg;
            break;

        case 'i':
            import_file = optarg;
            break;

        case 'o':
            output_file = optarg;
            break;

        case 'r':
            rec_step = atof(optarg);
            break;

        case 'm':
            max_time = atof(optarg);
            break;

        case 'f':
            inflow_file = optarg;
            break;

        case 'x':
            mode = optarg;
            break;

        case 's':
            single_file = true;
            break;

        case '?':
            /* getopt_long already printed an error message. */
            break;

        default:
            abort ();
        }
    }

    const std::string errmsg("insufficient arguments: --mode=[sph|asm|visualize] --in=<h5part-file> | --json=<json_file> [--out=<h5part-file> --rec_step=<record time step> --max_time=<max time> --inflow=<inflow-file>]\n");
    if (mode.isEmpty() || 0 == mode.compare("visualize")) {
#ifdef NOGUI
        throw std::runtime_error("Visual mode is only available in GUI binary (see define NOGUI)");
#endif
        throw std::runtime_error("Visual mode does not support command line parsing");
    } else if (0 == mode.compare("sph")) {
#ifndef NOGUI
        throw std::runtime_error("SPH mode is only available in NOGUI binary");
#endif
        if (json_file.isEmpty()) {
            throw std::runtime_error(errmsg + "json file needs to be provided");
        }
        if (output_file.isEmpty()) {
            throw std::runtime_error(errmsg + "output file needs to be provided");
        }
        if (rec_step <= 0.0) {
            throw std::runtime_error(errmsg + "rec step must be set to a value greater 0.0");
        }
        // set up scene
        QString json_file_path = json_file + ".json";
        Sceneparser::JSONParser<number> parser(json_file_path.toStdString(), "", output_file.toStdString(), rec_step, max_time, inflow_file.toStdString());
        Particlegenerator::ParticleHelper<number> hlp(parser.getConfiguration());
        auto tmp = hlp.getParticles();
        auto scene = new Scene2(Scene2::eCalculationMode::eCalcSPH, parser.getConfiguration(), tmp, single_file);
        return scene;
    }
    else {
        throw std::runtime_error(errmsg + "mode must be either set to sph, asm, or visualize");
    }
}

void nogui(int argc, char *argv[]) {
    Scene2 *scene = applyArguments(argc, argv);

    number current_time = 0.0;
    size_t step = 0;

    while (current_time < scene->max_time() || scene->max_time() == 0.0) {
        if (false == scene->iterate(step, current_time))
            break;

        step++;
    }

    scene->exportNumSteps();

    delete scene;
}

int main(int argc, char *argv[]) {

#ifdef NOGUI
    nogui(argc, argv);
    return 0;
#else
    QApplication app(argc, argv);

    QCoreApplication::setOrganizationDomain("uibk.ac.at");
    QCoreApplication::setOrganizationName("IUT");
    QCoreApplication::setApplicationName("Sphase");
    QCoreApplication::setApplicationVersion("0.1.0");

    QGLFormat glf = QGLFormat::defaultFormat();
    glf.setSampleBuffers(true);
    glf.setSamples(8);
    QGLFormat::setDefaultFormat(glf);

    MainWindow window;
    window.show();
    app.connect(&app, SIGNAL(lastWindowClosed()), &app, SLOT(quit()));
    app.exec();
    return 0;
#endif
}
