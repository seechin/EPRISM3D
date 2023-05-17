//>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
//---------------------------------------------------------------------------------
//-----------------------------   HI Equation Solver   ----------------------------
//---------------------------------------------------------------------------------
//---------------------------------------------------------------------------------
class HIEquationSolver {
  private:
    bool trim;
    int count;
    double * xa, * ya;
  public:
    double lse_a, lse_b;
  public:
    void set_trim(){ trim = true; }
    void unset_trim(){ trim = false; }
    double f(double x){
        double pre = lse_b * lse_a * exp((1-1/x)/lse_a);
        return ln(x+1e-15) + pre;
    }
    void dispose(){
        if (count>0){ free(xa); free(ya); }
        count = 0; xa = ya = nullptr;
    }
    void prepare(double x_inf, double x_sup){
        for (int i=0; i<=count; i++){
            xa[i] = x_inf + (x_sup - x_inf) * i / (count);
            ya[i] = this->f(xa[i]);
//printf("%12g %12g\n", xa[i], ya[i]);
        }
    }
    void init(int N = 10000){
        lse_a = 0.3; lse_b = 51.8;
        count = N; trim = false;
        xa = (double*) memalloc(sizeof(double) * (count+1) * 2);
        ya = &xa[count+1];
    }
    void set_param(double _lse_a, double _lse_b, double _inf=0, double _sup=2, bool _set_trim=true){
        lse_a = _lse_a; lse_b = _lse_b;
        trim = _set_trim; prepare(_inf, _sup);
    }
    double getx_ds(double y, double xinf, double xsup, double err){
        double inf = xinf; double sup = xsup;
        while (sup-inf>err){
            double xt = (inf + sup) / 2;
            double yt = f(xt);
            if (yt >= y) sup = xt;
            if (yt <= y) inf = xt;
        }
        return (inf + sup) / 2;
    }
    double getx_ds(double y, double err=1e-5){
        double inf = 0; double sup = 0;
        while (f(inf) > y) inf -= 1;
        while (f(sup) < y) sup += 1;
        return getx_ds(y, inf, sup, err);
    }
    double getx_bs(double y, double err=1e-5, bool further_ds=false){
        if (y>ya[count] || y<ya[0]){
            if (!trim) return getx_ds(y, err);
            else return y>ya[count]? xa[count] : xa[0];
        } else {
            int i = 0; int j = count;
            while (j-i>1){
                int k = (i+j)/2;
                if (ya[k]>=y) j = k;
                if (ya[k]<=y) i = k;
            }
            if (j==i) return xa[i];
            if (err > xa[j] - xa[i]){
                return xa[i];
            } else {
                if (further_ds){
                    return getx_ds(y, xa[i], xa[j], err);
                } else {
                    if (ya[i]==ya[j]) return (xa[i] + xa[j])/2;
                    return (xa[i]*(ya[j]-y) + xa[j]*(y-ya[i])) / (ya[j] - ya[i]); //
                }
            }
        }
    }
};

//>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
class RDFGroup { public:
  // index: >0: index; ==0: mathces anything; <0: ignore.
  // mol and atom name: can be string; or leave blank ("") or use "*" to match anything
    int is, iv, grp; // index of solute and solvent. Begin with 1.
    char ms[MAX_NAME], mv[MAX_NAME], as[MAX_NAME], av[MAX_NAME]; // mol and atom name, ignored if
    void init(int _is, int _iv, StringNS::string _ms, StringNS::string _as, StringNS::string _mv, StringNS::string _av, int _grp=-1){
        memset(this, 0, sizeof(RDFGroup));
        is = _is; iv = _iv; grp = _grp;
        memcpy(ms, _ms.text, _ms.length>(MAX_NAME-1)?(MAX_NAME-1):_ms.length);
        memcpy(as, _as.text, _as.length>(MAX_NAME-1)?(MAX_NAME-1):_as.length);
        memcpy(mv, _mv.text, _mv.length>(MAX_NAME-1)?(MAX_NAME-1):_mv.length);
        memcpy(av, _av.text, _av.length>(MAX_NAME-1)?(MAX_NAME-1):_av.length);
    }
};
class RDF_data { public:
    int is, iv; // index of solute and solvent, must be positive
    double * g, * n;
};
class RDF_datas { public:
    int n_rdf_datas[2]; RDF_data * rdf_datas[2];
};
