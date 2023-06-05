const char * software_name = "gmxtop2solute";
const char * software_version = "1.2.4.331";
const char * copyright_string = "(c) 2023 Cao Siqin";

#include    <errno.h>
#include    <stdio.h>
#include    <stdlib.h>
#include    <stdint.h>
#include    <string.h>
#include    <math.h>
#include    <signal.h>
#include    <fcntl.h>
#include    <ctype.h>
#include    <time.h>
#include    <sys/time.h>
#include    <sys/types.h>
#include    <sys/wait.h>
#include    <sys/stat.h>
#include    <sys/mman.h>
#include    <sys/resource.h>

#include    "header.h"
#include    "main-header.h"
#include    "String2.cpp"

#include    "Element.h"
#include    "main-atom-lists.h"
#include    "read_top.h"


//>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
//>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
//>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
const char * szHelp = "\
  The input/output files:\n\
    -p, -top              topology file, TOP\n\
    -ffpath, -include     forcefield folder, multiple separated with \":\"\n\
    -o                    output file, default: screen\n\
    -debug                show debug info\n\
    -excl                 exclude group, default: SOL\n\
    -use-atom-name        (-an) use atom name, not atom type (default)\n\
    -use-atom-type        (-at) use atom type, not atom name\n\
    -solvent-format       (-for-gensolvent) output as solvent format\n\
    -bond, -no-bond[s]    show/hide bond information, default on\n\
    -original-ri          (-ri) use the residue numbers in the top file\n\
    -reindex-residue[s]   (-rr) reindex residue numbers (default)\n\
    -default              reset all options\n\
  The output format:\n\
    mole_name atom_name mass charge sigma epsilon\n\
";


//>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
//>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
//>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

int debug_level = 0; bool show_atom_spc_name = true; bool reindex_residues = true; bool solvent_format = false;
    bool allow_bond = true; bool allow_index = true;
char * info_file_name = (char*)"";
char szfn_ffpath[MAX_PATH];
char szfn_top[MAX_PATH]; int i_param_szfn_top = -1;
char szfn_out[MAX_PATH];
char excl_grps[10][MAX_PATH]; int n_excl_grp = 0;
int analysis_parameter_line(char * argv[], int * argi, int argc, char * script_name, int script_line){
    int ret = 0; int i = *argi; bool analysis_script = !script_name? false : (!script_name[0]? false : true);
    StringNS::string key = argv[i];
    if (!analysis_script && (key == "-h" || key == "-help" || key == "--h" || key == "--help")){ ret = 2;
    } else if (!analysis_script && (key == "-version" || key == "--version")){ ret = 3;
    } else if (key == "-p" || key == "--p" || key == "-top" || key == "--top"){ if (i+1<argc){ i++; strcpy(szfn_top, argv[i]); i_param_szfn_top = i; }
    } else if (key == "-o" || key == "--o"){ if (i+1<argc){ i++; strcpy(szfn_out, argv[i]); }
    } else if (key == "-excl"){
        if (i+1<argc){ i++;
            if (n_excl_grp<10){
                strcpy(excl_grps[n_excl_grp], argv[i]);
                n_excl_grp ++;
            }
        }
    } else if (key == "-debug"){ debug_level = 1;
    } else if (key=="-ffpath" || key=="--ffpath" || key=="-ffpath" || key=="-include" || key=="--include" || key=="include"){ if (i+1<argc){
        i++; strcpy(szfn_ffpath,argv[i]); for (int j=0; j<MAX_PATH && szfn_ffpath[j]; j++) if (szfn_ffpath[j]==':') szfn_ffpath[j] = 0;
      }
    } else if (key=="-an"||key=="--an"||key=="-use-atom-name"||key=="--use-atom-name"||key=="-use_atom_name"||key=="--use_atom_name"){
        show_atom_spc_name = true;
    } else if (key=="-at"||key=="--at"||key=="-use-atom-type"||key=="--use-atom-type"||key=="-use_atom_type"||key=="--use_atom_type"){
        show_atom_spc_name = false;
    } else if (key=="-reindex-residue"||key=="--reindex-residue"||key=="-reindex_residue"||key=="--reindex_residue"||key=="-reindex-residues"||key=="--reindex-residues"||key=="-reindex_residues"||key=="--reindex_residues"||key=="-rr"||key=="--rr"){
        reindex_residues = true;
    } else if (key=="-original-residue-index"||key=="--original-residue-index"||key=="-original_residue_index"||key=="--original_residue_index"||key=="-original-ri"||key=="--original-ri"||key=="-original_ri"||key=="--original_ri"||key=="-ri"||key=="--ri"){
        reindex_residues = false;
    } else if (key=="-solvent-format"||key=="--solvent-format"||key=="-solvent_format"||key=="--solvent_format"||key=="-for-gensolvent"||key=="--solvent-format"||key=="-for_gensolvent"||key=="--solvent_format"){
        solvent_format = true;
    } else if (key=="-nb"||key=="--nb"||key=="-no-bond"||key=="--no-bond"||key=="-no-bonds"||key=="--no-bonds"){
        allow_bond = false;
    } else if (key=="-bond"||key=="--bond"||key=="-bond"){
        allow_bond = true;
    } else if (key=="-ni"||key=="--ni"||key=="-no-index"||key=="--no-index"){
        allow_index = false;
    } else if (key=="-default"||key=="--default"||key=="-default-format"||key=="--default-format"){
        show_atom_spc_name = true; solvent_format = false;
        allow_bond = false; allow_index = true;
    } else {
        i_param_szfn_top = i;
        strcpy(szfn_top, argv[i]);
    }
    *argi = i;
    return ret;
}
int analysis_params(int argc, char * argv[]){
    bool success = true; int error = 0;
  // analysis command line params
    if (argc<2){ printf("%s %s\n", software_name, software_version); return 0; }
    for (int i=1; i<argc; i++){
        if (argv[i][0] == '#' || argv[i][0] == ';'){
        } else {
            int _error = analysis_parameter_line(argv, &i, argc, (char*)"", i);
            if (_error){ success = false; error |= _error; }
        }
    }
    if (error == 2){
        printf("%s %s %s\n", software_name, software_version, copyright_string);
        printf("%s", szHelp);
        printf("%s", szLicence);
    } else if (error == 3){
        printf("%s\n", software_version);
    }
    if (!success) return error;
  // prepare other params

    return 0;
}



int main(int argc, char * argv[]){
    bool success = true; int error = 0; szfn_out[0] = 0;
    #ifdef DISTRIBUTE_VERSION
        software_version = DISTRIBUTE_VERSION;
    #endif

    error = analysis_params(argc, argv); if (error) return error;
    if (error) success = false;

    FILE * fout = stdout;
    if (success){
        if (StringNS::string(szfn_out)=="con" || StringNS::string(szfn_out)=="stdout" || StringNS::string(szfn_out)=="screen"){
            fout = stdout;
        } else if (StringNS::string(szfn_out)=="stderr"){
            fout = stderr;
        } else if (szfn_out[0]) {
            fout = fopen(szfn_out, "w");
            if (!fout){ fprintf(stderr, "%s : error : cannot write to %s\n", software_name, szfn_out); success = false; }
        }
    }

    if (success && szfn_top[0]){
        AnalysisTopParameters atp;
        if (szfn_ffpath[0]) atp.init(argc, argv, szfn_ffpath); else atp.init(argc, argv);
        atp.flog = stderr; atp.fout = fout;
        atp.solvent_format = solvent_format;
        atp.reindex_residues = reindex_residues;
        atp.show_atom_spc_name = show_atom_spc_name;
        atp.allow_bond = allow_bond;
        atp.debug_level = debug_level;

        atp.analysis_top(szfn_top, "arg", i_param_szfn_top, nullptr);

        if (fout && fout!=stdout && fout!=stderr) fclose(fout);

        atp.dispose();
    }
    return 0;
}
