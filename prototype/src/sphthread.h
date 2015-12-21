#ifndef SPHTHREAD_H
#define SPHTHREAD_H

#include <QThread>
#include <mutex>

class Scene2;

class SPHThread : public QThread {
    Q_OBJECT
public:
    SPHThread(Scene2 *scene);

    void run();

    inline void finishedDrawing() {
        drawing_mutex.unlock();
    }

    void stop() {
        running = false;
    }

    void toggleRunning() {
        if (running) {
            running = false;
        }
        else {
            running = true;
            this->start();
        }
    }

    void setTime(double t) {
        time = t;
    }

    double getTime() const {
        return time;
    }

    long getTimeSteps() const {
        return time_steps;
    }

    void forceStep(int step){
        this->step = step;
    }

signals:
    void dump(double time, long time_step);
    void stepUpdated(long step, bool finished);

private:
    double time;
    double dump_time;
    bool running;
    Scene2 *scene;

    long time_steps;
    long step;

    std::mutex drawing_mutex;
};

#endif // SPHTHREAD_H
