#include "sphthread.h"
#include <QDebug>
#include "scene2.h"
#include <chrono>

using namespace std::chrono;

SPHThread::SPHThread(Scene2 *scene)
    : scene(scene), running(false), time(0), time_steps(0), step(0) {
}

void SPHThread::run() {
    while (running && (scene->max_time() == 0.0 || time < scene->max_time()))
    {
        high_resolution_clock::time_point t1 = high_resolution_clock::now();
        double dt = time;
        bool canceled = !scene->iterate(step, time);
        drawing_mutex.lock();
        emit stepUpdated(step, canceled);
        if (canceled)
            break;

        step = step + 1;
        dt = time - dt;
        time_steps++;

        high_resolution_clock::time_point t2 = high_resolution_clock::now();

        auto duration = std::chrono::duration_cast<std::chrono::microseconds>( t2 - t1 ).count();
        auto delay = scene->getFrameDelay() * 1000000.0;

        if (delay > 0 && delay > duration) {
            QThread::msleep((delay - duration) / 1000.0);
        }
    }
}
