#ifndef MAINWINDOW_H
#define MAINWINDOW_H

#include <QMainWindow>
#include <set>
#include <map>
#include <QFileSystemWatcher>

class Scene2;

namespace Ui {
class MainWindow;
}

class MainWindow : public QMainWindow
{
    Q_OBJECT

public:
    explicit MainWindow(QWidget *parent = 0);
    ~MainWindow();

    void closeEvent(QCloseEvent *event);

private slots:
    void on_selected_valueChanged(int id);

    void particle_selection_from_renderer(int id);

    void cell_selection_from_renderer(std::set<int> ids);

    void select_cells_from_renderer();

    void show_grid_from_renderer();

    void on_min_value_valueChanged(double arg1);

    void on_max_value_valueChanged(double arg1);

    void on_value_comboBox_activated(int index);

    void on_colorWallParticles_clicked(bool checked);

    void on_colorSelectedCell_clicked(bool checked);

    void on_fluidSeparationColors_clicked(bool checked);

    void on_drawVectors_clicked(bool checked);

    void on_drawDots_clicked(bool checked);

    void on_vector_scale_valueChanged(double arg);

    void on_normalizeVectors_clicked(bool checked);

    void on_velVectorSlider_valueChanged(int value);

    void on_perParticle_clicked(bool checked);

    void on_perCell_clicked(bool checked);

    void on_buttonPlay_clicked();

    void on_progressSlider_sliderMoved();

    void on_progressSpin_valueChanged(int value);

    void on_frameDelay_valueChanged(double value);

    void stepUpdated(long step, bool canceled);

    void on_rangeButton_clicked();

    void on_pushButtonSS_clicked(bool checked);

    void on_pushButtonSNH_clicked(bool checked);

    void on_pushButtonSNO_clicked(bool checked);

    void on_pushButtonXS_clicked(bool checked);

    void on_pushButtonXBH_clicked(bool checked);

    void on_pushButtonSND_clicked(bool checked);

    void on_pushButtonXBA_clicked(bool checked);

    void on_pushButtonXND_clicked(bool checked);

    void on_pushButtonSO_clicked(bool checked);

    void on_pushButtonXP_clicked(bool checked);

    void on_pushButtonReactor1_clicked(bool checked);

    void on_pushButtonReactor2_clicked(bool checked);

    void on_pushButtonSelectCells_clicked(bool checked);

    void on_pushButtonShowGrid_clicked(bool checked);

    void on_pushButtonGridParticles_clicked(bool checked);

    void on_pushButtonExportASM_clicked(bool checked);

    void on_sphfilebutton_clicked();

    void on_asmfilebutton_clicked();

    void show_grid_particles_from_renderer();

    void showModified();

    void on_actionStartSimulation_triggered();

    void on_progressSlider_valueChanged(int value);

    void on_pushButtonReactorAreas_clicked(bool checked);
private:
    Ui::MainWindow *ui;
    Scene2 *scene;
    std::set<int> selectedCells;
    std::map<int, QColor> colorMap;
    QFileSystemWatcher* watcher;

    void displayParticleValues(int id);
    void setScene(Scene2* sc);
    void updatePlaybackControls();
    Scene2* createScene(std::string sph_file_path);
    void updateTimestep(int value);
};

#endif // MAINWINDOW_H
