#include "configurationparser.h"
#include <QJsonDocument>
#include <QJsonObject>
#include <QJsonValue>
#include <QJsonArray>

namespace Sceneparser {


namespace
{
template <typename number>
QJsonArray parseLineSegmentVector(std::vector<LineSegment<number>> tmpLineSegmentVector) {
    QJsonArray list;
    LineSegment<number> lineStruct;
    for (int i = 0; i < tmpLineSegmentVector.size(); i++) {
        QJsonObject map;

        lineStruct = tmpLineSegmentVector[i];

        QJsonArray velocity;
        velocity.append(lineStruct.wall_velocity.x);
        velocity.append(lineStruct.wall_velocity.y);
        map["velocity"] = velocity;

        map["moving"] = lineStruct.is_moving;

        QJsonObject from;
        from["x"] = lineStruct.b.x;
        from["y"] = lineStruct.b.y;
        map["from"] = from;

        QJsonObject to;
        to["x"] = lineStruct.m.x;
        to["y"] = lineStruct.m.y;
        map["to"] = to;

        list.append(map);
    }

    return list;
}

template <typename number>
QJsonArray parseAreaStructVector(const Configuration<number>& c, std::vector<AreaStruct<number>> tmpAreaStructVector) {
    QJsonArray list;
    AreaStruct<number> areaStruct;

    for (int i = 0; i < tmpAreaStructVector.size(); i++) {
        QJsonObject area;
        areaStruct = tmpAreaStructVector[i];
        area["x"] = areaStruct.x + c.scene.originx;
        area["y"] = areaStruct.y + c.scene.originy;
        area["width"] = areaStruct.width;
        area["height"] = areaStruct.height;
        list.push_back(area);
    }

    return list;
}
}

template <typename number>
std::string ConfigurationParser<number>::getJson() {
    return json;
}

template <typename number>
ConfigurationParser<number>::ConfigurationParser(const Configuration<number> configuration)
    : c(configuration) {
    createJson();
}

template <typename number>
void ConfigurationParser<number>::createJson() {

    // TODO check if fields are available

    QJsonObject root;

    root["export_file"] = QString::fromStdString(c.export_file);
    root["import_file"] = QString::fromStdString(c.import_file);
    root["inflow_file"] = QString::fromStdString(c.inflow_file);
    root["max_time"] = c.max_time;
    root["json_file"] = QString::fromStdString(c.json_file);
    root["rec_step"] = c.rec_step;

    QJsonObject scene;
    scene["height"] = c.scene.height;
    scene["hydrostatic_pressure"] = c.scene.hydrostatic_pressure;
    scene["neighbours"] = c.scene.neighbours;
    scene["origin_x"] = c.scene.originx;
    scene["origin_y"] = c.scene.originy;
    scene["sampling_dist"] = c.scene.sample_dist;
    scene["width"] = c.scene.width;
    root["scene"] = scene;

    QJsonObject sph;
    sph["alpha"] = c.sph.alpha;
    sph["c1"] = c.sph.c1;
    sph["c2"] = c.sph.c2;
    sph["d_air"] = c.sph.D_air;
    sph["d_water"] = c.sph.D_water;
    sph["epsilon_repulsion"] = c.sph.epsilon_repulsion;
    sph["epsilon_xsph"] = c.sph.epsilon_xsph;
    QJsonArray g;
    g.append(c.sph.g.x);
    g.append(c.sph.g.y);
    sph["g"] = g;
    sph["gamma1"] = c.sph.gamma1;
    sph["gamma2"] = c.sph.gamma2;
    sph["moving_least_squares"] = c.sph.moving_least_squares;
    sph["no_slip"] = c.sph.no_slip;
    sph["nu1"] = c.sph.nu1;
    sph["nu2"] = c.sph.nu2;
    sph["origin_x"] = c.sph.originx;
    sph["origin_y"] = c.sph.originy;
    sph["rho0_1"] = c.sph.rho0_1;
    sph["rho0_2"] = c.sph.rho0_2;
    sph["sd"] = c.sph.sd;
    sph["shepard"] = c.sph.shepard;
    sph["t_damp"] = c.sph.t_damp;
    sph["volume"] = c.sph.volume;
    QJsonObject stirrer;
    stirrer["is_active"] = c.sph.stirrer.is_active;
    stirrer["length_x"] = c.sph.stirrer.length_x;
    stirrer["length_y"] = c.sph.stirrer.length_y;
    stirrer["origin_x"] = c.sph.stirrer.origin_x;
    stirrer["origin_y"] = c.sph.stirrer.origin_y;
    stirrer["type"] = c.sph.stirrer.type;
    stirrer["velocity"] = c.sph.stirrer.velocity;
    QJsonObject stirrer2;
    stirrer2["is_active"] = c.sph.stirrer.is_active;
    stirrer2["length_x"] = c.sph.stirrer.length_x;
    stirrer2["length_y"] = c.sph.stirrer.length_y;
    stirrer2["origin_x"] = c.sph.stirrer.origin_x;
    stirrer2["origin_y"] = c.sph.stirrer.origin_y;
    stirrer2["type"] = c.sph.stirrer.type;
    stirrer2["velocity"] = c.sph.stirrer.velocity;
    sph["stirrer"] = stirrer;
    sph["stirrer2"] = stirrer2;
    root["sph"] = sph;

    QJsonObject sm;
    sm["ba"] = c.sm.ba;
    sm["bh"] = c.sm.bh;
    sm["eg"] = c.sm.eg;
    sm["eh"] = c.sm.eh;
    sm["fp"] = c.sm.fp;
    sm["height"] = c.sm.height;
    sm["ixb"] = c.sm.ixb;
    sm["ixe"] = c.sm.ixe;
    sm["ixp"] = c.sm.ixp;
    sm["ka"] = c.sm.ka;
    sm["kh"] = c.sm.kh;
    sm["Knh"] = c.sm.Knh;
    sm["Kno"] = c.sm.Kno;
    sm["Koa"] = c.sm.Koa;
    sm["Koh"] = c.sm.Koh;
    sm["Ks"] = c.sm.Ks;
    sm["Kx"] = c.sm.Kx;
    sm["ma"] = c.sm.ma;
    sm["mh"] = c.sm.mh;
    sm["rho7"] = c.sm.rho7;
    sm["sd"] = c.sm.sd;
    sm["SRT"] = c.sm.SRT;
    sm["time_scaling"] = c.sm.time_scaling;
    sm["width"] = c.sm.width;
    sm["ya"] = c.sm.ya;
    sm["yh"] = c.sm.yh;

    QJsonObject grid;
    grid["origin_x"] = c.sm.grid.origin_x;
    grid["origin_y"] = c.sm.grid.origin_y;
    grid["size"] = c.sm.grid.size;
    grid["neighbours"] = c.sm.grid.neighbours;
    sm["grid"] = grid;

    QJsonObject inflow_area;
    inflow_area["x"] = c.sm.inflow_area.x;
    inflow_area["y"] = c.sm.inflow_area.y;
    inflow_area["width"] = c.sm.inflow_area.width;
    inflow_area["height"] = c.sm.inflow_area.height;
    sm["inflow_area"] = inflow_area;

    QJsonObject outflow_area;
    outflow_area["x"] = c.sm.outflow_area.x;
    outflow_area["y"] = c.sm.outflow_area.y;
    outflow_area["width"] = c.sm.outflow_area.width;
    outflow_area["height"] = c.sm.outflow_area.height;
    sm["outflow_area"] = outflow_area;

    QJsonObject reactor1;
    reactor1["oxygen"] = c.sm.reactor1.oxygen;
    reactor1["cstr"] = c.sm.reactor1.cstr;
    reactor1["cstr_volume"] = c.sm.reactor1.cstr_volume;
    reactor1["inflow"] = c.sm.reactor1.inflow;
    QJsonObject area1;
    area1["x"] = c.sm.reactor1.area.x;
    area1["y"] = c.sm.reactor1.area.y;
    area1["width"] = c.sm.reactor1.area.width;
    area1["height"] = c.sm.reactor1.area.height;
    reactor1["area"] = area1;
    sm["reactor1"] = reactor1;

    QJsonObject reactor2;
    reactor2["oxygen"] = c.sm.reactor2.oxygen;
    reactor2["cstr"] = c.sm.reactor2.cstr;
    reactor2["cstr_volume"] = c.sm.reactor2.cstr_volume;
    reactor2["inflow"] = c.sm.reactor2.inflow;
    QJsonObject area2;
    area2["x"] = c.sm.reactor2.area.x;
    area2["y"] = c.sm.reactor2.area.y;
    area2["width"] = c.sm.reactor2.area.width;
    area2["height"] = c.sm.reactor2.area.height;
    reactor2["area"] = area2;
    sm["reactor2"] = reactor2;

    root["asm"] = sm;

    root["fluid1"] = parseAreaStructVector(c, c.fluid1);
    root["fluid2"] = parseAreaStructVector(c, c.fluid2);
    root["fluid_reactor1"] = parseAreaStructVector(c, c.fluid_reactor1);
    root["fluid_reactor2"] = parseAreaStructVector(c, c.fluid_reactor2);

    root["solid_boundaries"] = parseLineSegmentVector(c.solid_boundaries);

    root["periodic_in"] = parseLineSegmentVector(c.periodic_in);
    root["periodic_out"] = parseLineSegmentVector(c.periodic_out);
    root["rigid_periodic_in"] = parseLineSegmentVector(c.rigid_periodic_in);
    root["rigid_periodic_out"] = parseLineSegmentVector(c.rigid_periodic_out);

    root["horizontal_in"] = parseLineSegmentVector(c.horizontal_in);
    root["bottom_in"] = parseLineSegmentVector(c.bottom_in);
    root["open_out"] = parseLineSegmentVector(c.open_out);
    root["vprofile"] = parseLineSegmentVector(c.vprofile);
    root["q_counter"] = parseLineSegmentVector(c.qcounter);
    root["h_counter"] = parseLineSegmentVector(c.hcounter);

    QJsonDocument doc(root);
    json = QString::fromUtf8(doc.toJson().data()).toStdString();
}

template class ConfigurationParser<float>;
template class ConfigurationParser<double>;

}
