/*
 *
 * Author:      Schuyler D. Smith
 * Function:    versioning
 * Purpose:     creates verson information
 *
 */

#ifndef VERSION
#define VERSION

#include <string>

std::string versioning() {
    std::string version = "\n";
    version += "smartchip_analyzer version: ";
    version += "1.0.1";
    version += "  ";
    version += "March 31, 2023";
    version += "\n\n";
    return version;
}

#define PROGRAM_VERSION versioning()

#endif // VERSION