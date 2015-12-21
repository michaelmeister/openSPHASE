#include "particlerenderer.h"

#include <QDebug>
#include <QKeyEvent>

#include "scene2.h"
#include "sphthread.h"
#include <boost/foreach.hpp>
#include <H5Part.h>

#include <QTemporaryFile>
#include <QDir>
#include <QGLShader>
#include <QGLShaderProgram>
#include <QGLFunctions>
#ifdef __APPLE__
#include <OpenGL/glu.h>
#else
#include <GL/glu.h>
#endif

#ifdef _WIN32
// undefine what is already defined in qglshaderprogram.h
#undef GL_LINES_ADJACENCY_EXT
#undef GL_LINE_STRIP_ADJACENCY_EXT
#undef GL_TRIANGLES_ADJACENCY_EXT
#undef GL_TRIANGLE_STRIP_ADJACENCY_EXT
#include <glext\glext.h>

PFNGLMULTITEXCOORD1FPROC glMultiTexCoord1f = NULL;
PFNGLMULTITEXCOORD2FPROC glMultiTexCoord2f = NULL;
#endif // _WIN32

#define CGL \
{ \
    GLenum error = glGetError(); \
    if (error != GL_NO_ERROR) { \
    /*cout << __FILE__ << ":" << __LINE__ << gluErrorString(error) << std::endl;*/ \
    } \
}

ParticleRenderer::ParticleRenderer(QWidget *parent)
    : QGLViewer(parent), scene(0), selected(-1)
{
    connect(&fps_timer, SIGNAL(timeout()), SLOT(printFPS()));
    fps_timer.setInterval(1000);

    setMouseBinding(Qt::NoModifier, Qt::LeftButton, NO_CLICK_ACTION);

    setSnapshotFileName("out.png");
    setSnapshotFormat("PNG");
    setSnapshotCounter(0);

    render_value = Pressure;
    rv_min = 0;
    rv_max = 1;
    color_wall_particles = false;
    fluidSeparationColors = false;
    vectorNumberFactor = 1;
    velocityPerParticle = true;
    showGrid = false;
    showGridParticles = false;
    selectCells = false;
    color_selected_cell = false;
    selectedCells.clear();
    isRunning = false;
    normalizeVectors = true;
    vectorScale = 1.0;
}

void ParticleRenderer::mousePressEvent(QMouseEvent* e) {
    QGLViewer::mousePressEvent(e);
}

ParticleRenderer::~ParticleRenderer() {
    //delete scene; //WHAT THE FUCK destructor is called before closeEvent ... NO PLAN WTF happens here
}

void ParticleRenderer::setSelected(int selected, bool update_gl) {
    // return if selection has not changed
    if (this->selected == selected) {
        return;
    }

    // unset old selection flag for particle
    if (this->selected >= 0) {
        Particle2 *p = scene->getParticles()->at(this->selected);
        p->unset_selected();
    }

    // assign new selection index to particle renderer
    this->selected = selected;

    // set new selection flag for particle
    if (selected >= 0) {
        Particle2 *p = scene->getParticles()->at(this->selected);
        p->set_selected();
    }

    if (update_gl)
        updateGL();
}

#define SAMPLES 512
void ParticleRenderer::prepareRainbowTexture() {
    glEnable(GL_TEXTURE_1D);

    GLubyte black[1][4];


    black[0][0] = 0.0;
    black[0][1] = 0.0;
    black[0][2] = 0.0;
    black[0][3] = 255;

    glGenTextures(1, &black_texture);
    CGL
            glBindTexture(GL_TEXTURE_1D, black_texture);

    CGL
            glTexImage1D(GL_TEXTURE_1D, 0, GL_RGBA, 1, 0, GL_RGBA, GL_UNSIGNED_BYTE, black);
    CGL
            Q_ASSERT(glGetError() == GL_NO_ERROR);
    CGL
            glTexParameteri(GL_TEXTURE_1D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
    glTexParameteri(GL_TEXTURE_1D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
    glTexParameteri(GL_TEXTURE_1D, GL_TEXTURE_WRAP_S, GL_CLAMP_TO_EDGE);
    CGL
            Q_ASSERT(glIsTexture(black_texture));

    GLubyte blue[1][4];


    blue[0][0] = 0;
    blue[0][1] = 0.0;
    blue[0][2] = 255.0;
    blue[0][3] = 180.0;

    glGenTextures(1, &blue_texture);
    CGL
            glBindTexture(GL_TEXTURE_1D, blue_texture);

    CGL
            glTexImage1D(GL_TEXTURE_1D, 0, GL_RGBA, 1, 0, GL_RGBA, GL_UNSIGNED_BYTE, blue);
    CGL
            Q_ASSERT(glGetError() == GL_NO_ERROR);
    CGL
            glTexParameteri(GL_TEXTURE_1D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
    glTexParameteri(GL_TEXTURE_1D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
    glTexParameteri(GL_TEXTURE_1D, GL_TEXTURE_WRAP_S, GL_CLAMP_TO_EDGE);
    CGL
            Q_ASSERT(glIsTexture(blue_texture));

    GLubyte grey[1][4];


    grey[0][0] = 200.0;
    grey[0][1] = 200.0;
    grey[0][2] = 200.0;
    grey[0][3] = 200.0;

    glGenTextures(1, &grey_texture);
    CGL
            glBindTexture(GL_TEXTURE_1D, grey_texture);

    CGL
            glTexImage1D(GL_TEXTURE_1D, 0, GL_RGBA, 1, 0, GL_RGBA, GL_UNSIGNED_BYTE, grey);
    CGL
            Q_ASSERT(glGetError() == GL_NO_ERROR);
    CGL
            glTexParameteri(GL_TEXTURE_1D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
    glTexParameteri(GL_TEXTURE_1D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
    glTexParameteri(GL_TEXTURE_1D, GL_TEXTURE_WRAP_S, GL_CLAMP_TO_EDGE);
    CGL
            Q_ASSERT(glIsTexture(grey_texture));

    GLubyte red[1][4];


    red[0][0] = 255.0;
    red[0][1] = 0.0;
    red[0][2] = 0.0;
    red[0][3] = 255.0;

    glGenTextures(1, &red_texture);
    CGL
            glBindTexture(GL_TEXTURE_1D, red_texture);

    CGL
            glTexImage1D(GL_TEXTURE_1D, 0, GL_RGBA, 1, 0, GL_RGBA, GL_UNSIGNED_BYTE, red);
    CGL
            Q_ASSERT(glGetError() == GL_NO_ERROR);
    CGL
            glTexParameteri(GL_TEXTURE_1D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
    glTexParameteri(GL_TEXTURE_1D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
    glTexParameteri(GL_TEXTURE_1D, GL_TEXTURE_WRAP_S, GL_CLAMP_TO_EDGE);
    CGL
            Q_ASSERT(glIsTexture(red_texture));

    GLubyte data[SAMPLES][4];

    double h_inc = 0.7 / float(SAMPLES);

    for (int i = 0; i < SAMPLES; i++) {
        QColor c = QColor::fromHsvF(h_inc*i+0.3, 1.0, 1.0);
        data[i][0] = c.red();
        data[i][1] = c.green();
        data[i][2] = c.blue();
        data[i][3] = c.alpha();
    }

    glGenTextures(1, &rainbow_texture);
    CGL
            glBindTexture(GL_TEXTURE_1D, rainbow_texture);
    CGL
            glTexImage1D(GL_TEXTURE_1D, 0, GL_RGBA, SAMPLES, 0, GL_RGBA, GL_UNSIGNED_BYTE, data);
    CGL
            Q_ASSERT(glGetError() == GL_NO_ERROR);
    CGL
            glTexParameteri(GL_TEXTURE_1D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
    glTexParameteri(GL_TEXTURE_1D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
    glTexParameteri(GL_TEXTURE_1D, GL_TEXTURE_WRAP_S, GL_CLAMP_TO_EDGE);
    CGL
            Q_ASSERT(glIsTexture(rainbow_texture));

    glDisable(GL_TEXTURE_1D);
}

bool ParticleRenderer::isRun()
{
    return isRunning;
}

void ParticleRenderer::finishedDrawing()
{
    thread->finishedDrawing();
}

void ParticleRenderer::toggleStartSimulation()
{
    if (isRunning)
        isRunning = false;
    else
        isRunning = true;

    timer.isActive() ? timer.stop() : timer.start();
    fps_timer.isActive() ? fps_timer.stop() : fps_timer.start();
    thread->toggleRunning();
}

void ParticleRenderer::keyPressEvent(QKeyEvent *e) {
    /*if (e->key() == Qt::Key_Space)
    {
        toggleStartSimulation();
        return;
    }*/
    if (e->key() == Qt::Key_M) {
        if (scene->measure == true)
            scene->measure = false;
        else
            scene->measure = true;
        return;
    }
    if (e->key() == Qt::Key_V) {
        scene->handle_velocity_profile();
        return;
    }
    if (e->key() == Qt::Key_Q) {
        scene->handle_qcounter();
        return;
    }
    if (e->key() == Qt::Key_P) {
        if (selected >= 0)
            cout << scene->particles->at(selected)->pos.x << ", " << scene->particles->at(selected)->pos.y << endl;
        return;
    }
    if (e->key() == Qt::Key_G) {
        emit showGridChanged();
    }
    if (e->key() == Qt::Key_C) {
        emit selectCellsChanged();
    }
    if (e->key() == Qt::Key_T) {
        toggleTextIsEnabled();
        updateGL();
    }
    if (e->key() == Qt::Key_F) {
        emit showGridParticlesChanged();
    }

    QGLViewer::keyPressEvent(e);
}

static const char *vertex_source =
        "void main(void) {"
        "   gl_FrontColor = gl_Color;"
        "   gl_Position = gl_ModelViewProjectionMatrix * gl_Vertex;"
        "   gl_TexCoord[0]  = gl_MultiTexCoord0;"
        "   gl_TexCoord[1]  = gl_MultiTexCoord1;"
        "}"
        ;

static const char *fragment_source =
        "uniform sampler1D myTexture;"
        "void main(void) {"
        "   vec4 tc = gl_TexCoord[1];"
        "   vec4 c = texture1D(myTexture, gl_TexCoord[0].s);"
        "   tc.s -= 0.5;"
        "   tc.t -= 0.5;"
        "   if ((tc.s * tc.s + tc.t * tc.t) < 0.25)"
        "       gl_FragColor = c;"
        "   else"
        "       gl_FragColor = vec4(0,0,0,0);"
        "}"
        ;

void ParticleRenderer::init() {
    glClearColor(1.0, 1.0, 1.0, 0.0); CGL
            glPointSize(4); CGL
            glInitNames(); CGL
            glDisable(GL_LIGHTING); CGL
            glEnable(GL_POINT_SMOOTH); CGL
            glEnable(GL_BLEND); CGL
            glEnable(GL_DEPTH_TEST); CGL
            glDepthFunc(GL_LEQUAL); CGL

            glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA); CGL
            glHint(GL_POINT_SMOOTH_HINT, GL_NICEST); CGL

            prepareRainbowTexture();

    bool r;

    vertex = new QGLShader(QGLShader::Vertex);
    r = vertex->compileSourceCode(vertex_source);
    Q_ASSERT(r);

    fragment = new QGLShader(QGLShader::Fragment);
    r = fragment->compileSourceCode(fragment_source);
    Q_ASSERT(r);

    prog = new QGLShaderProgram();
    r = prog->addShader(vertex);
    Q_ASSERT(r);

    r = prog->addShader(fragment);
    Q_ASSERT(r);

    r = prog->link();
    Q_ASSERT(r);

    glEnable(GL_MULTISAMPLE);

    GLint bufs;
    GLint samples;
    glGetIntegerv(GL_SAMPLE_BUFFERS, &bufs);
    glGetIntegerv(GL_SAMPLES, &samples);
    qDebug("Have %d buffers and %d samples", bufs, samples);

#ifdef _WIN32
    glMultiTexCoord1f = (PFNGLMULTITEXCOORD1FPROC)wglGetProcAddress("glMultiTexCoord1f");
    glMultiTexCoord2f = (PFNGLMULTITEXCOORD2FPROC)wglGetProcAddress("glMultiTexCoord2f");
#endif // _WIN32
}

void ParticleRenderer::drawGrid() {
    Q_ASSERT(scene);
    glColor3f(0.7, 0.7, 0.7);
    glLineWidth(0.5);
    glBegin(GL_LINES);
    double width = scene->configuration.scene.width;
    double height = scene->configuration.scene.height;
    double gridOriginX = scene->configuration.sm.grid.origin_x;
    double gridOriginY = scene->configuration.sm.grid.origin_y;
    double dx = scene->configuration.scene.sample_dist*scene->configuration.sm.grid.size;
    double startX = fmod(gridOriginX, dx);
    double startY = fmod(gridOriginY, dx);

    for (double x = startX; x < width; x += dx) {
        glVertex2f(x, 0);
        glVertex2f(x, height);
    }

    for (double y = startY; y < height; y += dx) {
        glVertex2f(0, y);
        glVertex2f(width, y);
    }

    glEnd();

}


void ParticleRenderer::setup(Scene2 *scene) {
    setSceneBoundingBox(qglviewer::Vec(0, 0, 0), qglviewer::Vec(scene->configuration.scene.width, scene->configuration.scene.height, 0.0));
    showEntireScene();
    this->scene = scene;
    thread = new SPHThread(scene);
    QObject::connect(thread, SIGNAL(stepUpdated(long, bool)), this, SIGNAL(stepUpdated(long, bool)), Qt::QueuedConnection);
}

void ParticleRenderer::getRenderValueRange(double &min, double &max) {
    if (render_value == Velocity) {
        scene->sph.getVelocityRange(min, max);
    }

    if (render_value == Density) {
        scene->sph.getDensityRange(min, max);
    }

    if (render_value == Pressure) {
        scene->sph.getPressureRange(min, max);
    }
}

void ParticleRenderer::renderBoundaries() const {
    glColor3f(.0, .0, .0);

    auto sceneconfiguration = scene->configuration;
    std::vector<LineSegment> boundaries(sceneconfiguration.solid_boundaries);
    boundaries.insert(boundaries.end(), sceneconfiguration.periodic_in.begin(), sceneconfiguration.periodic_in.end());
    boundaries.insert(boundaries.end(), sceneconfiguration.periodic_out.begin(), sceneconfiguration.periodic_out.end());
    boundaries.insert(boundaries.end(), sceneconfiguration.rigid_periodic_in.begin(), sceneconfiguration.rigid_periodic_in.end());
    boundaries.insert(boundaries.end(), sceneconfiguration.rigid_periodic_out.begin(), sceneconfiguration.rigid_periodic_out.end());
    boundaries.insert(boundaries.end(), sceneconfiguration.horizontal_in.begin(), sceneconfiguration.horizontal_in.end());
    boundaries.insert(boundaries.end(), sceneconfiguration.open_out.begin(), sceneconfiguration.open_out.end());
    boundaries.insert(boundaries.end(), sceneconfiguration.vprofile.begin(), sceneconfiguration.vprofile.end());
    boundaries.insert(boundaries.end(), sceneconfiguration.qcounter.begin(), sceneconfiguration.qcounter.end());
    boundaries.insert(boundaries.end(), sceneconfiguration.hcounter.begin(), sceneconfiguration.hcounter.end());

    for (std::vector<LineSegment>::iterator it = boundaries.begin();
         it != boundaries.end();
         ++it) {
        glBegin(GL_LINES);
        LineSegment&ls = *it;
        Vector2 end = ls.b + ls.m;
        glVertex2f(ls.b.x, ls.b.y);
        glVertex2f(end.x, end.y);
        glEnd();
    }
}

void ParticleRenderer::renderSceneBoundaries() const {
    glColor3f(.9, .9, .9);

    glBegin(GL_LINE_STRIP);
    glVertex2f(0 - this->scene->sph.getParameters().originx, 0 - this->scene->sph.getParameters().originy);
    glVertex2f(0 - this->scene->sph.getParameters().originx, this->scene->configuration.scene.height - this->scene->sph.getParameters().originy);
    glVertex2f(this->scene->configuration.scene.width - this->scene->sph.getParameters().originx, this->scene->configuration.scene.height - this->scene->sph.getParameters().originy);
    glVertex2f(this->scene->configuration.scene.width - this->scene->sph.getParameters().originx, 0 - this->scene->sph.getParameters().originy);
    glVertex2f(0 - this->scene->sph.getParameters().originx, 0 - this->scene->sph.getParameters().originy);
    glEnd();
}

void ParticleRenderer::updateTime(double time)
{
    thread->setTime(time);
    updateGL();
}

void ParticleRenderer::renderTime() {
    glColor3f(.0f, .0f, .0f);
    QFont f;
    QFontMetrics fm(f);
    renderText(10, height() - fm.height(), QString("time %1 ").arg(thread->getTime()));
}

void ParticleRenderer::setDrawGrid(bool value) {
    showGrid = value;
    updateGL();
}

void ParticleRenderer::setReactorAreas(bool value) {
    //showReactorAreas = value;
    updateGL();
}

void ParticleRenderer::setDrawGridParticles(bool value) {
    showGridParticles = value;
    updateGL();
}

void ParticleRenderer::draw() {
    if (!scene) return;
    QMutexLocker locker(&sph_data); //protect sph data against multiple access
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
    glMatrixMode(GL_MODELVIEW);

    CGL

            glColor3f(0.0, 0.0, 0.0);

    if (showGrid) {
        drawGrid();
    }

    if (this->drawDots) {
        drawParticleDots(showGridParticles);
    }

    if (drawVectors) {
        if (velocityPerParticle) {
            drawPerParticleVelocities();
        } else {
            drawPerCellVelocities();
        }
    }

    if (selected < 0 || selected >= scene->particles->size()) {
        return;
    }

//    BOOST_FOREACH(NeighbourData &nd, scene->neighbours[selected]) {
//        Particle2 *p = scene->particles->at(nd.j);
//        if (p == scene->particles->at(selected)) {
//            glColor3f(1.0, 0.0, 0.0);
//        } else {
//            glColor3f(0.0, 0.0, 1.0);
//        }
//        glBegin(GL_POINTS);
//        glVertex2f(p->pos.x, p->pos.y);
//        glEnd();
//    }
}

bool ParticleRenderer::drawParticleDot(Particle2 *p, int i, number l, bool isGridParticle = false, float range = 0) {
    bool is_selected = false;
    if (p->selected) {
        is_selected = true;
        if (i != selected)
            setSelected(i, false);
    }

    glPushName(i);

    if (!color_wall_particles && p->is_wall()) {
        glBindTexture(GL_TEXTURE_1D, black_texture);
    }
    else if (fluidSeparationColors && p->is_second_phase()) {
        glBindTexture(GL_TEXTURE_1D, red_texture);
    }
    else if (isGridParticle) {
        glBindTexture(GL_TEXTURE_1D, black_texture);
    }
    else {
        glBindTexture(GL_TEXTURE_1D, rainbow_texture);
    }

    if (showGrid && color_selected_cell && p->is_fluid())
        glBindTexture(GL_TEXTURE_1D, red_texture);

    //glColor3f(0.0, 1.0, 0.0);
    double value;
    switch (render_value) {
    case Velocity:
        value = p->v.length();
        break;
    case Pressure:
        value = p->p;
        break;
    default:
    case Density:
        value = p->rho;
        break;
    }

    value = (value - rv_min) * range;

    Vector2 pos = p->pos;
    glBegin(GL_QUADS);
    glMultiTexCoord2f(GL_TEXTURE1, .0, .0); glMultiTexCoord1f(GL_TEXTURE0, value); glVertex3f(pos.x - l, pos.y - l, 0.0);
    glMultiTexCoord2f(GL_TEXTURE1, 1., .0); glMultiTexCoord1f(GL_TEXTURE0, value); glVertex3f(pos.x + l, pos.y - l, 0.0);
    glMultiTexCoord2f(GL_TEXTURE1, 1., 1.); glMultiTexCoord1f(GL_TEXTURE0, value); glVertex3f(pos.x + l, pos.y + l, 0.0);
    glMultiTexCoord2f(GL_TEXTURE1, 0., 1.); glMultiTexCoord1f(GL_TEXTURE0, value); glVertex3f(pos.x - l, pos.y + l, 0.0);
    glEnd();
    glPopName();
    CGL

    return is_selected;
}

void ParticleRenderer::drawNeighbourParticles(Particle2 *p, number l) {
    glBindTexture(GL_TEXTURE_1D, red_texture);

    double value = 0.0;

    Vector2 pos = p->pos;
    glBegin(GL_QUADS);
    glMultiTexCoord2f(GL_TEXTURE1, .0, .0); glMultiTexCoord1f(GL_TEXTURE0, value); glVertex3f(pos.x - l, pos.y - l, 0.0);
    glMultiTexCoord2f(GL_TEXTURE1, 1., .0); glMultiTexCoord1f(GL_TEXTURE0, value); glVertex3f(pos.x + l, pos.y - l, 0.0);
    glMultiTexCoord2f(GL_TEXTURE1, 1., 1.); glMultiTexCoord1f(GL_TEXTURE0, value); glVertex3f(pos.x + l, pos.y + l, 0.0);
    glMultiTexCoord2f(GL_TEXTURE1, 0., 1.); glMultiTexCoord1f(GL_TEXTURE0, value); glVertex3f(pos.x - l, pos.y + l, 0.0);
    glEnd();
    glPopName();
    CGL
}

void ParticleRenderer::drawParticleDots(bool drawGridParticles)
{
    float range = 1.0 / (rv_max - rv_min);
    number l = scene->configuration.scene.sample_dist / 2;
    bool r = prog->bind();
    bool is_selected = false;
    Q_ASSERT(r);

    glEnable(GL_TEXTURE_1D);

    for (long long i = 0; i < (long long)scene->particles->size(); i++) {
        Particle2 *p = scene->particles->at(i);
        if (drawParticleDot(p, i, l, false, range))
            is_selected = true;
    }


    if (!is_selected)
        setSelected(-1, false);

//    for (int i = 0; i < selectedNeighbours.size(); i++) {
//        std::vector<int> neighbours = scene->a.getGridParticleNeighbours(selectedNeighbours[i]);
//        if (!neighbours.empty()) {
//            for (int j = 0; j < neighbours.size(); j++) {
//                Particle2 *n = scene->particles->at(neighbours.at(j));
//                drawNeighbourParticles(n, l);
//            }
//        }
//    }

    glDisable(GL_TEXTURE_1D);

    QGLFunctions f(this->context());
    f.glUseProgram(0);
}

void drawVector(const Vector2& pos, const Vector2& shaft)
{
    const double arrow_length = 0.6;
    const double arrow_width = 0.2;

    const double leftx = shaft.x * arrow_length + shaft.y * arrow_width;
    const double lefty = shaft.y * arrow_length - shaft.x * arrow_width;

    const double rightx = shaft.x * arrow_length - shaft.y * arrow_width;
    const double righty = shaft.y * arrow_length + shaft.x * arrow_width;

    glLineWidth(2.5f);
    glBegin(GL_LINE_STRIP);
    glVertex2f(pos.x, pos.y);
    glVertex2f(pos.x + shaft.x, pos.y + shaft.y);
    glVertex2f(pos.x + leftx, pos.y + lefty);
    glEnd();
    glBegin(GL_LINES);
    glVertex2f(pos.x + rightx, pos.y + righty);
    glVertex2f(pos.x + shaft.x, pos.y + shaft.y);
    glEnd();
    glLineWidth(1.0f);
}

void ParticleRenderer::drawPerParticleVelocities()
{
    for (size_t i = 0; i < scene->particles->size(); i++) {
        if (i % vectorNumberFactor != 0)
            continue;

        Particle2 *p = scene->particles->at(i);

        if (p->is_wall())
            continue;

        glPushName((int)i);

        Vector2 head = p->v;
        if (normalizeVectors)
            head = head.normalized(scene->configuration.sph.c1);

        drawVector(p->pos, head * vectorScale);
        glPopName();
    }
}

// not beautiful, but this way it gets initialized only once
std::vector<Particles> g_cells;

void ParticleRenderer::drawPerCellVelocities()
{
    // init cells
    Sceneparser::Vector2<size_t> grid_dim = {(size_t)std::ceil(scene->configuration.scene.width / (scene->configuration.scene.neighbours*scene->configuration.scene.sample_dist)), (size_t)std::ceil(scene->configuration.scene.height / (scene->configuration.scene.neighbours*scene->configuration.scene.sample_dist))};

    // we need a cell density above the unigrid, lets create one
    grid_dim = grid_dim / vectorNumberFactor;
    g_cells.resize(grid_dim.x*grid_dim.y);
    number gridOriginX = 0.0; // scene->getGridOriginX();
    number gridOriginY = 0.0; // scene->getGridOriginY();

    for (int p_i = 0; p_i < scene->particles->size(); p_i++)
    {
        Particle2 *p = scene->particles->at(p_i);
        number rel_pos_x = p->pos.x - gridOriginX;
        number rel_pos_y = p->pos.y - gridOriginY;
        rel_pos_x /= scene->configuration.scene.width;
        rel_pos_y /= scene->configuration.scene.height;
        const size_t idxx = (size_t)(rel_pos_x*grid_dim.x);
        const size_t idxy = (size_t)(rel_pos_y*grid_dim.y);
        g_cells[idxx + idxy * grid_dim.x].push_back(p);
    }

    Kernel* k = Kernel::get("wendland2", (std::min)(grid_dim.x, grid_dim.y) / 4.0);

    for (int c_i = 0; c_i < g_cells.size(); c_i++)
    {
        Particles& ps = g_cells.at(c_i);
        const size_t x = c_i % grid_dim.x;
        const size_t y = c_i / grid_dim.x;
        Vector2 cellCenter(	(x + 0.5) / grid_dim.x * scene->configuration.scene.width,
                                            (y + 0.5) / grid_dim.y * scene->configuration.scene.height);

        Vector2 avg_Velocity(0, 0);
        number sum_weighting = 0;
        number weighting = 0;
        for (int p_i = 0; p_i < ps.size(); p_i++)
        {
            Particle2* p = ps.at(p_i);
            if (p->is_wall())
                continue;
            weighting = k->f((cellCenter - p->pos).length());
            avg_Velocity += weighting * p->v;
            sum_weighting += weighting;
        }

        avg_Velocity = avg_Velocity / sum_weighting;

        if (normalizeVectors)
            avg_Velocity = avg_Velocity.normalized(scene->configuration.sph.c1);

        cellCenter.x += gridOriginX;
        cellCenter.y += gridOriginY;

        drawVector(cellCenter, avg_Velocity * vectorScale);
        // reset cell contend
        ps.clear();
    }

    delete k;
}

void ParticleRenderer::postDraw() {
    if (!scene) return;

    renderBoundaries();
    renderSceneBoundaries();
    startScreenCoordinatesSystem(true);

    glEnable(GL_TEXTURE_1D);
    glBindTexture(GL_TEXTURE_1D, rainbow_texture);

    assert(glGetError() == GL_NO_ERROR);
    glBegin(GL_QUADS);
    glColor3f(1.0, 1.0, 1.0);
    glTexCoord1f(0);        glVertex2f(width()*0.05, height()*0.1);
    glTexCoord1f(0);        glVertex2f(width()*0.05 + 40, height()*0.1);

    glTexCoord1f(1.0);      glVertex2f(width()*0.05 + 40, height()*0.4);
    glTexCoord1f(1.0);      glVertex2f(width()*0.05, height()*0.4);
    glEnd();

    CGL

            glDisable(GL_TEXTURE_1D);
    stopScreenCoordinatesSystem();
    renderTime();
}

void ParticleRenderer::postSelection(const QPoint& point) {
    qDebug() << selectedName();

    if (selectedName() < 0 || selectedName() > scene->particles->size()) {
        setSelected(-1);
    } else {
        setSelected(selectedName());
    }

    emit particleSelected(selected);
}

void ParticleRenderer::printFPS() {
    static long ts = 0;
    long new_ts = thread->getTimeSteps();
    cout << "FPS: " << new_ts - ts << endl;
    ts = new_ts;
}

void ParticleRenderer::setSelectCells(bool value) {
    selectCells = value;
}

void ParticleRenderer::closeEvent(QCloseEvent *e)
{
    thread->stop();

    while (!thread->isFinished()) {}
    e->accept();
}

void ParticleRenderer::setRenderValue(RenderValue value) {
    this->render_value = value;
}

void ParticleRenderer::setRenderValueRange(double min, double max) {
    rv_min = min;
    rv_max = max;
}

void ParticleRenderer::setVectorScale(double s)
{
    vectorScale = s;
    repaint();
}

void ParticleRenderer::setDrawVectors(bool enabled)
{
    drawVectors = enabled;
    repaint();
}

void ParticleRenderer::setDrawDots(bool enabled)
{
    drawDots = enabled;
    repaint();
}

void ParticleRenderer::colorWallParticles(bool color) {
    this->color_wall_particles = color;
}

void ParticleRenderer::colorSelectedCell(bool color) {
    this->color_selected_cell = color;
}

void ParticleRenderer::colorFluidsDifferent(bool color) {
    this->fluidSeparationColors = color;
}

void ParticleRenderer::setNormalizeVectors(bool enabled) {
    normalizeVectors = enabled;
    repaint();
}

void ParticleRenderer::setVectorNumberFactor(unsigned int mod) {
    vectorNumberFactor = mod;
    repaint();
}

void ParticleRenderer::setPerVelocityPerParticle(bool value) {
    velocityPerParticle = value;
    repaint();
}

void ParticleRenderer::forceStep(size_t step) {
    thread->forceStep((int)step);
}
