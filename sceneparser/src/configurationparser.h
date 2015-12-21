#ifndef CONFIGURATIONPARSER_H
#define CONFIGURATIONPARSER_H

#include "configuration.h"
#include <string>

namespace Sceneparser {

template <typename number>
class ConfigurationParser {
public:
    ConfigurationParser(const Configuration<number> configuration);
    std::string getJson();

private:
    const Configuration<number> c;

    std::string json;

    void createJson();
};

}

#endif // CONFIGURATIONPARSER_H
