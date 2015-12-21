#ifndef PARTICLERENDERER_H
#define PARTICLERENDERER_H

#include "stable.h"
#include <QGLViewer/qglviewer.h>
#include <QTimer>
#include <QMutex>
#include <set>

class Scene2;
class SPHThread;
struct H5PartFile;
class QGLShader;
class QGLShaderProgram;

class ParticleRenderer : public QGLViewer {
    Q_OBJECT
public:

    enum RenderValue {
        Density,
        Pressure,
        Velocity
    };

    ParticleRenderer(QWidget *parent = 0);
    virtual ~ParticleRenderer();
    void draw();
    void postDraw();
    void postSelection(const QPoint& point);
    void init();
    void setup(Scene2 *scene);
    void getRenderValueRange(double &min, double &max);
    void setRenderValue(RenderValue value);
    void setRenderValueRange(double min, double max);
    void setVectorScale(double s);
    void setDrawVectors(bool enabled);
    void setDrawDots(bool enabled);
    void setNormalizeVectors(bool enabled);
    void setVectorNumberFactor(unsigned int value); // [1-16]
    void setPerVelocityPerParticle(bool value);
    void setDrawGrid(bool value);
    void setDrawGridParticles(bool value);
    bool isRun();

    void finishedDrawing();

    /**
     * @brief should wall particles be colored or just rendered black
     * @param color
     */
    void colorWallParticles(bool color);
    void colorFluidsDifferent(bool color);
    void colorSelectedCell(bool color);

    void toggleStartSimulation();
    void forceStep(size_t step);

    void updateTime(double time);
    void setReactorAreas(bool value);
signals:
    void particleSelected(int id);
    void cellSelected(std::set<int> ids);
    void stepUpdated(long step, bool canceled);
    void showGridChanged();
    void selectCellsChanged();
    void showGridParticlesChanged();

public slots:
    void printFPS();
    void setSelected(int selected, bool update_gl = false);
    void drawGrid();
    void setSelectCells(bool value);
protected:
    void closeEvent(QCloseEvent *);

private:
    void drawParticleDots(bool drawGridParticles);
    bool drawParticleDot(Particle2 *p, int i, number l, bool isGridParticle, float range);
    void drawPerParticleVelocities();
    void drawPerCellVelocities();

    void keyPressEvent(QKeyEvent *);
    void mousePressEvent(QMouseEvent* e);

    void renderBoundaries() const;
    void renderSceneBoundaries() const;
    void renderTime();
    void prepareRainbowTexture();

    Scene2*		scene;
    QTimer		timer;
    QTimer		fps_timer;
    int			selected;
    QMutex		sph_data;
    bool		drag;
    QPoint		drag_start;
    GLuint		rainbow_texture;
    GLuint		black_texture;
    GLuint      blue_texture;
    GLuint      grey_texture;
    GLuint		red_texture;
    bool		color_wall_particles;
    bool        color_selected_cell;
    bool		fluidSeparationColors;
    SPHThread*	thread;
    QGLShader*	vertex;
    QGLShader*	fragment;
    QGLShaderProgram*	prog;
    RenderValue render_value;
    double		rv_min;
    double		rv_max;
    double		vectorScale;
    bool		drawVectors;
    bool		drawDots;
    bool		normalizeVectors;
    unsigned int vectorNumberFactor;
    bool		velocityPerParticle;
    bool        showGrid;
    bool        showGridParticles;
    bool        selectCells;
    bool        isRunning;
    std::set<int>    selectedCells;
    void drawNeighbourParticles(Particle2 *p, number l);
};

#endif // PARTICLERENDERER_H
