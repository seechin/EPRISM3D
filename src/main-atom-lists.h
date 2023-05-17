//>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
//---------------------------------------------------------------------------------
//---------------   System Param: Atom List (Solute and Solvent)   ----------------
//---------------------------------------------------------------------------------
//---------------------------------------------------------------------------------
class AtomSite { public:
  // fundamental parameters
    char name[MAX_NAME]; char mole[MAX_NAME]; const char * nele;
    double mass; double charge; double charge_esp;  // charge for ES field, charge_esp for reporting ES energy
    double sigma, sqrt_sigma; double epsilon, sqrt_epsilon; int iaa;
    int id, grp, multi; bool is_key;
    unsigned int reserved;
  // advanced parameters
    double reverse_rism_factor;
    double ren_charge, ren_bond, ren_dielect;
  // initialization
    void init(int _id, int _grp, char * _mole, int _iaa, char * _name, double _mass, double _charge, double _sigma, double _epsilon){
        memset(name, 0, sizeof(name)); int len = strlen(_name); if (len>MAX_NAME-1) len = MAX_NAME-1; memcpy(name, _name, len);
        memset(mole, 0, sizeof(mole)); len = strlen(_mole); if (len>MAX_NAME-1) len = MAX_NAME-1; memcpy(mole, _mole, len);
        nele = ElementNS::get_atom_element(name);
        iaa = _iaa; mass = _mass; charge = _charge; charge_esp = _charge;
        sigma = _sigma; sqrt_sigma = sqrt(fabs(sigma));
        epsilon = _epsilon<0? -1e-15 : _epsilon; sqrt_epsilon = (_epsilon<0? -1 : 1) * sqrt(fabs(epsilon));
        id = _id; grp = _grp; multi = 1; is_key = false; reserved = 0;
        reverse_rism_factor = -1;
        ren_charge = 0; ren_bond = 1; ren_dielect = 1;
        //printf("atom site %s.%s: %12f %12f %12f %12f\n", mole, name, charge, sigma, epsilon, sqrt_epsilon);
    }
};
const int MAX_BONDS_PER_SOLUTE_ATOM = 6;
class SoluteAtomSite { public:
    int index; char name[8]; int iaa; char mole[8];
    double mass, charge; double sigma, sqrt_sigma, epsilon, sqrt_epsilon;
    int nbond; int ibond[MAX_BONDS_PER_SOLUTE_ATOM];
    int i_sigma_list; int n_sigma_list; int i_epsilon_list; int n_epsilon_list;
    unsigned int reserved;
    void init(int _index, const char * _name, int _iaa, const char * _mole, double _mass, double _charge, double _sigma, double _epsilon){
        index = _index;   memset(name, 0, sizeof(name)); strncpy(name, _name, sizeof(name)-1);
        iaa = _iaa; memset(mole, 0, sizeof(mole)); strncpy(mole, _mole, sizeof(mole)-1);
        mass = _mass; charge = _charge;
        sigma = _sigma; sqrt_sigma = sqrt(fabs(sigma));
        epsilon = _epsilon<0? -1e-15 : _epsilon; sqrt_epsilon = (_epsilon<0? -1 : 1) * sqrt(fabs(epsilon));
        i_sigma_list = -1; n_sigma_list = 0; i_epsilon_list = -1; n_epsilon_list = 0;
        reserved = 0;
    }
};

//>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
//---------------------------------------------------------------------------------
//-----------------------   System Param: Atom Name List   ------------------------
//---------------------------------------------------------------------------------
//---------------------------------------------------------------------------------
class SYSITEM_AtomNameList {
  public:
    char name[MAX_NAME]; char mole[MAX_NAME]; int index; int grp; int iaa;
};
class SYSITEM_PairMapping { // mapping of gvv
  public:
    int grpi; int grpj; int col; int direction;
};
class SYSITEM_BondList {
  public:
    int grpi; int grpj; double bond; double bond_stdev;
    double weight; // weight must be (0,1], and all invalid numbers are recognized as 100%
};
class SYSITEM_ZetaList {
  public:
      //int iaai, iaaj; double zeta0, rc_zeta0;
      int iaai, iaaj; double deltaG_in_Kelvin, rc_zeta0;
};
//-------------------------------------
int search_atom_list(int return_col, SYSITEM_AtomNameList * al, int begin, int nal, StringNS::string mole, StringNS::string atom, int default_ret){
    for (int i=begin; i<nal; i++){
        bool mole_match = mole=="*" || mole==al[i].mole;
        bool atom_match = atom=="*" || atom==al[i].name;
        if (mole_match && atom_match) return return_col==1? al[i].index : return_col==2? al[i].grp : return_col==3? al[i].iaa : i;
    }
    return default_ret;
}
int search_atom_list(int return_col, SYSITEM_AtomNameList * al, int begin, int nal, const char * sz_mole, const char * sz_atom, int default_ret){
    StringNS::string mole = sz_mole? sz_mole : "*";
    StringNS::string atom = sz_atom? sz_atom : "*";
    return search_atom_list(return_col, al, begin, nal, mole, atom, default_ret);
}
int search_atom_list(int return_col, SYSITEM_AtomNameList * al, int begin, int nal, const char * compond, int default_ret){
    char mole_name[MAX_NAME]; char atom_name[MAX_NAME]; memset(mole_name, 0, sizeof(mole_name)); memset(atom_name, 0, sizeof(atom_name));
    for (int i=0; compond[i]; i++){
        if (compond[i]=='.'){
            memcpy(mole_name, compond, i>(MAX_NAME-1)?(MAX_NAME-1):i); if (mole_name[0]==0) strcpy(mole_name, "*");
            const char * compond_next = &compond[i+1]; int compond_next_len = strlen(compond_next);
            memcpy(atom_name, compond_next, compond_next_len>(MAX_NAME-1)?(MAX_NAME-1):compond_next_len); if (atom_name[0]==0) strcpy(atom_name, "*");
            break;
        } else if (compond[i]==':' && compond[i+1]==':'){
            memcpy(mole_name, compond, i>(MAX_NAME-1)?(MAX_NAME-1):i); if (mole_name[0]==0) strcpy(mole_name, "*");
            const char * compond_next = &compond[i+2]; int compond_next_len = strlen(compond_next);
            memcpy(atom_name, compond_next, compond_next_len>(MAX_NAME-1)?(MAX_NAME-1):compond_next_len); if (atom_name[0]==0) strcpy(atom_name, "*");
            break;
        }
    }
    int compond_len = strlen(compond);
    if (!atom_name[0]){ strcpy(mole_name, "*"); memcpy(atom_name, compond, compond_len>(MAX_NAME-1)?(MAX_NAME-1):compond_len); }
    if (atom_name[0]=='*') return default_ret;
    return search_atom_list(return_col, al, begin, nal, mole_name, atom_name, default_ret);
}

/*int search_atom_list_2(int return_col, SYSITEM_AtomNameList * al, int begin, int nal, StringNS::string compound, int default_ret){
    StringNS::string search_strings[2] = { "", "" };
    char sep_sv[4] = { '.', ':', '.', ':' };
    int nw = analysis_general_line(sep_sv, compound, search_strings, 2, false, true);
    if (nw==1){
        int imole = search_atom_list(return_col, al, begin, nal, search_strings[0], "*", -1);
        int iatom = search_atom_list(return_col, al, begin, nal, "*", search_strings[0], -1);
        if (iatom>=0) return iatom;
        else if (imole>=0) return imole;
        else return default_ret;
    } else if (nw==2){
        return search_atom_list(return_col, al, begin, nal, search_strings[0], search_strings[1], default_ret);
    }
    return default_ret;
}*/


int search_atom_list_index(SYSITEM_AtomNameList * al, int begin, int nal, const char * compond, int default_ret){
    return search_atom_list(1, al, begin, nal, compond, default_ret);
}
int search_atom_list_grp(SYSITEM_AtomNameList * al, int begin, int nal, const char * compond, int default_ret){
    return search_atom_list(2, al, begin, nal, compond, default_ret);
}
int search_mole_list(SYSITEM_AtomNameList * al, int begin, int nal, const char * sz_mole, int default_ret){
    if (StringNS::is_string_number(sz_mole)){
        int iaa = atoi(sz_mole);
        return iaa<=0 ? default_ret : iaa-1;
    } else {
        return search_atom_list(3, al, begin, nal, sz_mole, "*", default_ret);
    }
}
