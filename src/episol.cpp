const char * software_name = "episol";
const char * software_version = "1.2.4.331";
const char * copyright_string = "(c) 2023 Cao Siqin";

#include    "header.h"
#include    "main-header.h"

#if defined(_GROMACS4_) || defined(_GROMACS5_) || defined(_GROMACS2016_) || defined(_GROMACS2018_)
  #define _GROMACS_
#endif
#include    <errno.h>
#include    <stdio.h>
#include    <stdlib.h>
#include    <stdint.h>
#include    <string.h>
#include    <math.h>
#include    <signal.h>
#include    <fcntl.h>
#include    <ctype.h>
#include    <unistd.h>
#include    <time.h>
#include    <libgen.h>
#include    <dirent.h>
#include    <sys/time.h>
#include    <sys/types.h>
#include    <sys/wait.h>
#include    <sys/stat.h>
#include    <sys/mman.h>
#include    <sys/resource.h>


#include    "String2.cpp"
#include    "Vector.cpp"
#include    "Element.h"
#include    "PDBAtom.cpp"
#include    "read_frame_abr.cpp"
#ifdef _LIBZ_
  #include  <zlib.h>
#endif
#include    "fftw3.h"
#ifndef MAX_MEMORYS
    const int MAX_MEMORYS = 65536;
    int _memory_blk_total = 0; size_t _memory_total = 0; size_t _memory_last_allocated = 0; bool _ignore_memory_capacity = false;
#endif

namespace GENSOLVENT {
    #include "gensolvent.cpp"
}

namespace GMXTOP2SOLUTE {
    #include "gmxtop2solute.cpp"
}

namespace EPRISM3D {
    #include "eprism3d.cpp"
}

namespace TS4SDUMP {
    #include "ts4sdump.cpp"
}

const char * szHelp = "\
  Sub commands:\n\
    eprism3d              perform the 3D-RISM calculations\n\
    gmxtop2solute         translate GROMACS TOP file to solute file\n\
    gensolvent            generate the intital solvent file\n\
    ts4sdump              dump the data of a TS4S file (.ts4s) to text\
";

int main(int argc, char * argv[]){
    bool success = true; //FILE * flog = stdout;
    #ifdef DISTRIBUTE_VERSION
        software_version = DISTRIBUTE_VERSION;
    #endif
    GMXTOP2SOLUTE::software_version = software_version;
    GMXTOP2SOLUTE::copyright_string = copyright_string;
    GENSOLVENT::software_version = software_version;
    GENSOLVENT::copyright_string = copyright_string;
    EPRISM3D::software_version = software_version;
    EPRISM3D::copyright_string = copyright_string;
    TS4SDUMP::software_version = software_version;
    TS4SDUMP::copyright_string = copyright_string;

    StringNS::string key = ""; if (argc>1) key = argv[1];

    if (argc<=1 || key=="-h" || key=="-help" || key=="--help"){
        printf("%s %s %s\n", software_name, software_version, copyright_string);
        printf("%s\n", szHelp);
    } else if (key=="gmxtop2solute" || key=="top2solute" || key=="top"){
        return GMXTOP2SOLUTE::main(argc-1, &argv[1]);
    } else if (key=="gensolvent" || key=="gv"){
        return GENSOLVENT::main(argc-1, &argv[1]);
    } else if (key=="eprism3d" || key=="rism3d" || key=="rismhi3d" || key=="3d"){
        return EPRISM3D::main(argc-1, &argv[1]);
    } else if (key=="ts4sdump" || key=="ts4s" || key=="ts"){
        return TS4SDUMP::main(argc-1, &argv[1]);
    } else {
        fprintf(stderr, "%s : error : unrecognizable executive: \"%s\"\n", software_name, argv[1]); success = false;
    }

    if (!success) return 1; return 0;
}

// cat `ls *.h *.cpp | grep -v heatmap` | grep #define | grep \# | grep -v // | grep define | tr '(' ' ' | awk '{print $2}' | sort | uniq | grep -v nullptr | awk '{printf("#undef %s\n",$1)}' > undef.h


