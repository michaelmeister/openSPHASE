#ifndef JSONPARSER_H
#define JSONPARSER_H

#include "configuration.h"
#include <string>

namespace Sceneparser {

template <typename number>
class JSONParser {
public:
    JSONParser(const std::string& json_file, const std::string& import_file, const std::string& export_file, const number& rec_step, const number& max_time, const std::string& inflow_file);
    Configuration<number> getConfiguration();

private:
    const std::string json_file;
    const std::string import_file;
    const std::string export_file;
    const number rec_step;
    const number max_time;
    const std::string inflow_file;

    Configuration<number> c;

    void createConfiguration();
};

}

#endif // JSONPARSER_H
