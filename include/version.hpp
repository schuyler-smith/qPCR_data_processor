#ifndef VERSION
#define VERSION

#include <string>

std::string versioning() {
    std::string version = "\n";
    version += "smartchip_analyzer";
    version += " ";
    version += "version:";
    version += " ";
    version += "1.0.2";
    version += "  ";
    version += __DATE__;
    version += "\n\n";
    return version;
}

#define PROGRAM_VERSION versioning()

#endif // VERSION