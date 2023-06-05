#ifndef MACHINE_REASONABLE_ERROR

#define __REAL__  double
#define MACHINE_REASONABLE_ERROR 1e-12

//#define     _TTYPROMPTCOLOR_
//#define     _FUNCTION_EXPORT_


const double PI = 3.1415926535897932384626433832795;
const double EE = 2.7182818284590452353602874713527;
const double COULCOOEF = 138.9354846;

const int MAX_SOL = 100;                // Max atom site number
const int MAX_CMD_PARAMS = MAX_SOL;     // Max parameter number for a command
const int MAX_THREADS = 500;            // Max number of forks or threads
const int MAX_DIIS = 100;               // Max parameter number for a command
const int MAX_INCLUDE_RECURSIVE = 100;  // maximum include recursive levels

const int MAX_PATH = 1024;              // Maximum length of filename/path strings
const int MAX_NAME = 64;                // Maximum length of molecule/atom names
const int MAX_WORD = 1000;              // Maximum number of words in line analysis

const int MAX_RENAME_COUNT = 10000;     // Maximum number of renames
const int MAX_GVV_FILES = 3;

#ifdef PACKAGE_VERSION
  #define DISTRIBUTE_VERSION    PACKAGE_VERSION
#else
  #define DISTRIBUTE_VERSION    "1.2.5"
#endif
#ifndef DISTRIBUTE_VERSION
    const char * szLicence = "";
#else
    const char * szLicence = "\
 # EPRISM3D is free software. You can use, modify or redistribute under the\n\
 # terms of the GNU Lesser General Public License v3:\n\
 # https://www.gnu.org/licenses/lgpl-3.0.en.html\n\
";
#endif

#endif
