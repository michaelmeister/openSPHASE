#include "mainwindow.h"
#include "ui_mainwindow.h"
#include "scene2.h"
#include <QSettings>
#include <QDoubleSpinBox>
#include <QFileDialog>
#include <jsonparser.h>
#include <particles.h>
#include <H5Part.h>
#include <fstream>

MainWindow::MainWindow(QWidget *parent) :
    QMainWindow(parent),
    ui(new Ui::MainWindow), scene(nullptr) {
    ui->setupUi(this);

    this->addAction(ui->actionStartSimulation);
    ui->groupBoxASM->setEnabled(false);
    ui->groupBoxSPH->setEnabled(false);

    ui->sphfilebutton->setEnabled(true);
    ui->asmfilebutton->setEnabled(false);

    ui->textBrowserParams->setEnabled(true);

    QSettings settings;

    restoreGeometry(settings.value("mainWindowGeometry").toByteArray());
    restoreState(settings.value("mainWindowState").toByteArray());

    this->connect(ui->particle_renderer, SIGNAL(particleSelected(int)), SLOT(particle_selection_from_renderer(int)));
    this->connect(ui->particle_renderer, SIGNAL(cellSelected(std::set<int>)), SLOT(cell_selection_from_renderer(std::set<int>)));
    this->connect(ui->particle_renderer, SIGNAL(showGridChanged()), SLOT(show_grid_from_renderer()));
    this->connect(ui->particle_renderer, SIGNAL(showGridParticlesChanged()), SLOT(show_grid_particles_from_renderer()));
    this->connect(ui->particle_renderer, SIGNAL(selectCellsChanged()), SLOT(select_cells_from_renderer()));
    this->connect(ui->particle_renderer, SIGNAL(stepUpdated(long, bool)), SLOT(stepUpdated(long, bool)));
}

void MainWindow::setScene(Scene2* sc) {
    if (nullptr != this->scene)
        delete(this->scene);

    this->scene = sc;

    ui->particle_renderer->setup(scene);

    ui->selected->setMinimum(0);
    ui->selected->setMaximum((int)scene->particles->size());

    on_drawDots_clicked(ui->drawDots->isChecked());
    on_drawVectors_clicked(ui->drawVectors->isChecked());
    on_vector_scale_valueChanged(ui->vector_scale->value());
    on_normalizeVectors_clicked(ui->normalizeVectors->isChecked());

    int m = scene->getNumImportedSteps();
    ui->progressSlider->setMaximum(m);
    ui->progressSpin->setMaximum(m);
}

MainWindow::~MainWindow() {
    delete ui;
    delete scene;
}

void MainWindow::closeEvent(QCloseEvent *event) {
    QSettings settings;
    settings.setValue("mainWindowGeometry", saveGeometry());
    settings.setValue("mainWindowState", saveState());
}

void MainWindow::displayParticleValues(int id) {
    if (id >= 0 && scene->particles->size() > id) {
        Particle2 *p_i = scene->particles->at(id);
        ui->density->setText(QString::number(p_i->rho, 'f', 6));
        ui->pressure->setText(QString::number(p_i->p, 'f', 6));
        Vector2 v = p_i->v;
        ui->velocity->setText(QString("|%1| %2, %3").arg(v.length()).arg(v.x).arg(v.y));
        Vector2 f = p_i->f;
        ui->acceleration->setText(QString("|%1| %2, %3").arg(f.length()).arg(f.x).arg(f.y));
    } else {
        ui->density->setText("");
        ui->pressure->setText("");
        ui->velocity->setText("");
        ui->acceleration->setText("");
    }
}

void MainWindow::on_selected_valueChanged(int id) {
    displayParticleValues(id);
}

void MainWindow::particle_selection_from_renderer(int id) {
    ui->selected->setValue(id);
    displayParticleValues(id);
}

void MainWindow::cell_selection_from_renderer(std::set<int> ids) {
    selectedCells = ids;
}

void MainWindow::show_grid_from_renderer() {
    if (ui->pushButtonShowGrid->isChecked()) {
        ui->pushButtonShowGrid->setChecked(false);
        on_pushButtonShowGrid_clicked(false);
    } else {
        ui->pushButtonShowGrid->setChecked(true);
        on_pushButtonShowGrid_clicked(true);
    }
}

void MainWindow::show_grid_particles_from_renderer() {
    if (ui->pushButtonGridParticles->isChecked()) {
        ui->pushButtonGridParticles->setChecked(false);
        on_pushButtonGridParticles_clicked(false);
    } else {
        ui->pushButtonGridParticles->setChecked(true);
        on_pushButtonGridParticles_clicked(true);
    }
}

void MainWindow::select_cells_from_renderer() {
    if (ui->pushButtonSelectCells->isChecked()) {
        ui->pushButtonSelectCells->setChecked(false);
        on_pushButtonSelectCells_clicked(false);
    } else {
        ui->pushButtonSelectCells->setChecked(true);
        on_pushButtonSelectCells_clicked(true);
    }
}

void MainWindow::on_min_value_valueChanged(double min) {
    ui->particle_renderer->setRenderValueRange(min, ui->max_value->value());
    ui->particle_renderer->updateGL();
}

void MainWindow::on_max_value_valueChanged(double max) {
    ui->particle_renderer->setRenderValueRange(ui->min_value->value(), max);
    ui->particle_renderer->updateGL();
}

void MainWindow::on_value_comboBox_activated(int value) {
    ui->particle_renderer->setRenderValue((ParticleRenderer::RenderValue) value);
    ui->particle_renderer->updateGL();
    this->on_rangeButton_clicked();     // update range
}

void MainWindow::on_colorWallParticles_clicked(bool checked) {
    ui->particle_renderer->colorWallParticles(checked);
    ui->particle_renderer->updateGL();
}

void MainWindow::on_colorSelectedCell_clicked(bool checked) {
    ui->particle_renderer->colorSelectedCell(checked);
    ui->particle_renderer->updateGL();
}

void MainWindow::on_fluidSeparationColors_clicked(bool checked) {
    ui->particle_renderer->colorFluidsDifferent(checked);
    ui->particle_renderer->updateGL();
}

void MainWindow::on_drawVectors_clicked(bool checked)
{
    ui->vector_scale->setEnabled(checked);
    ui->scale_label->setEnabled(checked);
    ui->particle_renderer->setDrawVectors(checked);
}

void MainWindow::on_pushButtonShowGrid_clicked(bool checked)
{
    ui->particle_renderer->setDrawGrid(checked);
}

void MainWindow::on_pushButtonGridParticles_clicked(bool checked)
{
    ui->particle_renderer->setDrawGridParticles(checked);
}

void MainWindow::on_pushButtonExportASM_clicked(bool checked)
{
    //ui->particle_renderer->exportASMValues();
}

void MainWindow::on_pushButtonReactorAreas_clicked(bool checked)
{
    ui->particle_renderer->setReactorAreas(checked);
}

void MainWindow::on_drawDots_clicked(bool checked)
{
    ui->particle_renderer->setDrawDots(checked);
}

void MainWindow::on_vector_scale_valueChanged(double arg)
{
    ui->particle_renderer->setVectorScale(arg);
}

void MainWindow::on_normalizeVectors_clicked(bool checked)
{
    ui->particle_renderer->setNormalizeVectors(checked);
}

void MainWindow::on_velVectorSlider_valueChanged(int value)
{
    ui->velVectorNumberLabel->setText("x " + QString::number(value));
    ui->particle_renderer->setVectorNumberFactor(value);
}

void MainWindow::on_perParticle_clicked(bool checked)
{
    ui->velVectorSlider->setEnabled(false);
    ui->label_9->setEnabled(false);
    ui->particle_renderer->setPerVelocityPerParticle(checked);
}

void MainWindow::on_perCell_clicked(bool checked)
{
    ui->velVectorSlider->setEnabled(true);
    ui->label_9->setEnabled(true);
    ui->velVectorSlider->setMaximum(8);

    ui->particle_renderer->setPerVelocityPerParticle(!checked);
    on_velVectorSlider_valueChanged(ui->velVectorSlider->value());
}

void MainWindow::on_pushButtonSS_clicked(bool checked)
{
    //plotGraphs();
}

void MainWindow::on_pushButtonSNH_clicked(bool checked)
{
    //plotGraphs();
}

void MainWindow::on_pushButtonSNO_clicked(bool checked)
{
    //plotGraphs();
}

void MainWindow::on_pushButtonXS_clicked(bool checked)
{
    //plotGraphs();
}

void MainWindow::on_pushButtonXBH_clicked(bool checked)
{
    //plotGraphs();
}

void MainWindow::on_pushButtonSND_clicked(bool checked)
{
    //plotGraphs();
}

void MainWindow::on_pushButtonXBA_clicked(bool checked)
{
    //plotGraphs();
}

void MainWindow::on_pushButtonXND_clicked(bool checked)
{
    //plotGraphs();
}

void MainWindow::on_pushButtonSO_clicked(bool checked)
{
    //plotGraphs();
}

void MainWindow::on_pushButtonXP_clicked(bool checked)
{
    //plotGraphs();
}

void MainWindow::on_pushButtonReactor1_clicked(bool checked)
{
    //plotGraphs();
}

void MainWindow::on_pushButtonReactor2_clicked(bool checked)
{
    //plotGraphs();
}

void MainWindow::on_pushButtonSelectCells_clicked(bool checked) {
   ui->particle_renderer->setSelectCells(checked);
}


void MainWindow::on_buttonPlay_clicked()
{
    static QString start = ui->buttonPlay->toolTip();
    static QString stop = "Pause";
    QString text = ui->buttonPlay->toolTip();
    ui->buttonPlay->setIcon(QIcon(text == start ? ":/icons/pause" : ":/icons/play"));
    ui->buttonPlay->setToolTip(text == start ? stop : start);
    ui->particle_renderer->toggleStartSimulation();
    showModified();
}

void MainWindow::on_progressSlider_sliderMoved()
{
    ui->progressSpin->setValue(ui->progressSlider->value());
}

void MainWindow::on_progressSlider_valueChanged(int value)
{
    ui->progressSpin->setValue(value);
}

void MainWindow::on_progressSpin_valueChanged(int value)
{
    if (ui->progressSlider->value() != value)
        ui->progressSlider->setValue(value);
    if (!ui->particle_renderer->isRun())
        updateTimestep(value);
}

void MainWindow::updateTimestep(int value)
{
    ui->particle_renderer->forceStep(value);
    number time = 0.0;
    scene->ImportStepFromH5Part(value, time);
    ui->particle_renderer->updateTime(time);
}

void MainWindow::on_frameDelay_valueChanged(double value)
{
    scene->setFrameDelay(value);
}

void MainWindow::stepUpdated(long step, bool canceled)
{
    ui->particle_renderer->updateGL();
    if (ui->progressSlider->maximum() < step) {
        ui->progressSlider->setMaximum(step);
        ui->progressSpin->setMaximum(step);
    }
    ui->progressSlider->setValue((int)step);
    ui->progressSpin->setValue((int)step);

    if (canceled)
        on_buttonPlay_clicked();

    ui->particle_renderer->finishedDrawing();
}

void MainWindow::on_rangeButton_clicked() {
    double min, max;
    ui->particle_renderer->getRenderValueRange(min, max);
    ui->min_value->setValue(min);
    ui->max_value->setValue(max);
    double step = (max - min) / 100.0;
    ui->min_value->setSingleStep(step);
    ui->max_value->setSingleStep(step);
}

void MainWindow::on_actionStartSimulation_triggered() {
    on_buttonPlay_clicked();
}

void MainWindow::on_sphfilebutton_clicked()
{
    QString path = QFileDialog::getOpenFileName(this, tr("SPH file"));
    if ( path.isNull() == false && path.endsWith("h5part"))
    {
        try {
            auto s = createScene(path.toStdString());

            watcher = new QFileSystemWatcher(this);
            watcher->addPath(path);
            QObject::connect(watcher, SIGNAL(fileChanged(QString)), this, SLOT(showModified()));

            path.chop(7);
            s->addSPHImportFile(path.toStdString());
            setScene(s);

            ui->sphInputLineEdit->setText(path);
            ui->groupBoxSPH->setEnabled(true);
            ui->asmfilebutton->setEnabled(true);
            updatePlaybackControls(); // enable / disable playback functionality
            number current_time = 0.0;
            scene->ImportStepFromH5Part(0, current_time);
            ui->labelNumSteps->setText("Steps: " + QString::number(scene->getNumImportedSteps()));
            Sceneparser::ConfigurationParser<number> parser(scene->configuration);
            std::string json = parser.getJson();
            ui->textBrowserParams->setText(QString(json.c_str()));
            //plotGraphs();
            ui->groupBoxASM->setEnabled(false);
            ui->asmInputLineEdit->setText("");
        } catch (const runtime_error& error) {
            QMessageBox::critical(this, "Error", error.what());
        }
    } else {
        QMessageBox::critical(this, "Error", "File could not be opened!");
    }
}

void MainWindow::showModified()
{
    if (!ui->particle_renderer->isRun()) {
        ui->progressSlider->setEnabled(true);
        ui->progressSpin->setEnabled(true);
    }
}

Scene2* MainWindow::createScene(std::string sph_file_path) {
    H5PartFile* part_file = H5PartOpenFile(sph_file_path.c_str(), H5PART_READ);
    if (part_file)
    {
        h5part_int64_t json_size(0);
        h5part_int64_t status = H5PartReadFileAttrib(part_file, "json_blob_size", &json_size);
        if (status != H5PART_SUCCESS) {
            return nullptr;
        }

        char* json = new char[json_size+1];
        status = H5PartReadFileAttrib(part_file, "json_blob", json);
        if (status != H5PART_SUCCESS) {
            return nullptr;
        }
        json[json_size] = '\0';

        H5PartCloseFile(part_file);

        std::string fileName = "tmpConfigurationFile.json";

        std::ofstream outfile;
        outfile.open(fileName, ios::out);
        outfile.write(json, sizeof(char)*json_size);
        outfile.close();

        delete[] json;

        // load scene
        Sceneparser::JSONParser<number> parser(fileName, "", "", 0, 0, "");
        Particlegenerator::ParticleHelper<number> hlp(parser.getConfiguration());
        auto s = new Scene2(Scene2::eCalculationMode::eVisualize, parser.getConfiguration(), hlp.getParticles(), false);
        return s;
    } else {
        throw std::runtime_error("SPH file could not be opened!");
    }
    return nullptr;
}

void MainWindow::on_asmfilebutton_clicked()
{
    QString path = QFileDialog::getOpenFileName(this, tr("ASM file"));
    if ( path.isNull() == false  && path.endsWith("h5part"))
    {
        try {
            scene->addASMImportFile(path.toStdString());
            ui->asmInputLineEdit->setText(path);
            ui->groupBoxASM->setEnabled(true);
            Sceneparser::ConfigurationParser<number> parser(scene->configuration);
            std::string json = parser.getJson();
            ui->textBrowserParams->setText(QString(json.c_str()));
            updatePlaybackControls(); // enable / disable playback functionality
            //plotGraphs();
        } catch (const runtime_error& error) {
            QMessageBox::critical(this, "Error", error.what());
        }
    } else {
        QMessageBox::critical(this, "Error", "ASM file could not be opened!");
    }
}


void MainWindow::updatePlaybackControls()
{
    bool enabled = (nullptr != this->scene);

    this->ui->buttonPlay->setEnabled(enabled);
    this->ui->delayLabel->setEnabled(enabled);
    this->ui->frameDelay->setEnabled(enabled);
    this->ui->progressSlider->setEnabled(enabled);
    this->ui->progressSpin->setEnabled(enabled);
}
