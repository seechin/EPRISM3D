namespace RISMHI3D_RISMNS {
    //>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    //>>>>>>>                             >>>>>>>>>>>>>>>>>>>>>>>>>
    //>>>>>>>  General Closure Treatment  >>>>>>>>>>>>>>>>>>>>>>>>>
    //>>>>>>>                             >>>>>>>>>>>>>>>>>>>>>>>>>
    //>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

    double perform_closure(IET_Param * sys, size_t N3, size_t i3begin, size_t i3end, __REAL__ **** auuv, __REAL__ **** aulr, __REAL__ **** acuv, __REAL__ **** ahuv, __REAL__ **** aextra, __REAL__ **** ares, __REAL__ **** ahlr, __REAL__ **** aclr, __REAL__ **** add){
            // extra: extra content for e.g. HNCB, RBC
        double ret = 0;
        for (int iv=0; iv<sys->nv; iv++){
            __REAL__ * uuv = &auuv[iv][0][0][0];
            __REAL__ * ulr = &aulr[iv][0][0][0];
            __REAL__ * cuv = &acuv[iv][0][0][0];
            __REAL__ * huv = &ahuv[iv][0][0][0];
            __REAL__ * res = &ares[iv][0][0][0];
            __REAL__ * hlr = ahlr? &ahlr[iv][0][0][0] : nullptr;
            __REAL__ * clr = aclr? &aclr[iv][0][0][0] : nullptr;
            __REAL__ * extra = aextra? &aextra[iv][0][0][0] : nullptr;
            __REAL__ * dd = add? &add[sys->nvm][0][0][0] : nullptr;
            double factor = sys->closure_factors[iv];
            double cutoff = sys->ccutoff;
            double t = 0;
            int closure = sys->closures[iv];
            //if (i3begin<0) i3begin=0; if (i3end>N3) i3end = N3;
            if (i3end>N3) i3end = N3;

            if (closure==CLOSURE_NONE){
                for (size_t i3=i3begin; i3<i3end; i3++) res[i3] = 0;
            } else if (closure==CLOSURE_HNC){   // HNC, can be scaled with -cf
                for (size_t i3=i3begin; i3<i3end; i3++){

                    t = factor*(huv[i3] - cuv[i3]);
                    t = exp(-uuv[i3] + t) - 1;
                    res[i3] = t - huv[i3]; huv[i3] = t;
                }
            } else if (closure==CLOSURE_PLHNC){
                double expC = exp(factor);
                for (size_t i3=i3begin; i3<i3end; i3++){
                    t = -uuv[i3] + (huv[i3] - cuv[i3]);
                    t = (t>factor? t+expC-factor : exp(t)) - 1;
                    res[i3] = t - huv[i3]; huv[i3] = t;
                }
            } else if (closure==CLOSURE_HARDSPHERE){
                for (size_t i3=i3begin; i3<i3end; i3++){
                    t = uuv[i3]>cutoff? uuv[i3] : 0;
                    t = exp(-t + (huv[i3] - cuv[i3])) - 1;
                    res[i3] = t - huv[i3]; huv[i3] = t;
                }
            } else if (closure==CLOSURE_MSA){
                for (size_t i3=i3begin; i3<i3end; i3++){
                    if (uuv[i3]>cutoff){
                        //t = exp(-uuv[i3] + (huv[i3] - cuv[i3])) - 1;
                        t = -1;
                        res[i3] = t - huv[i3]; huv[i3] = t;
                    } else {
                        res[i3] = - uuv[i3] - cuv[i3];
                    }
                }
            } else if (closure==CLOSURE_KGK){
                for (size_t i3=i3begin; i3<i3end; i3++){
                    t = -uuv[i3] + (huv[i3] - cuv[i3]);
                    t = t>-1? t : -1;
                    res[i3] = t - huv[i3]; huv[i3] = t;
                }
            } else if (closure==CLOSURE_KH){
                for (size_t i3=i3begin; i3<i3end; i3++){
                    t = -uuv[i3] + (huv[i3] - cuv[i3]);
                    t = t>0? t : exp(t) - 1;
                    res[i3] = t - huv[i3]; huv[i3] = t;
                }
            } else if (closure==CLOSURE_PSE2){
                for (size_t i3=i3begin; i3<i3end; i3++){
                    t = -uuv[i3] + huv[i3] - cuv[i3];
                    if (t>0){
                        t = t + t*t/2;
                    } else {
                        t = exp(t) - 1;
                    }
                    res[i3] = t - huv[i3]; huv[i3] = t;
                }
            } else if (closure==CLOSURE_PSE3){
                for (size_t i3=i3begin; i3<i3end; i3++){
                    t = -uuv[i3] + huv[i3] - cuv[i3];
                    if (t>0){
                        t = t + t*t/2 + t*t*t/6;
                    } else {
                        t = exp(t) - 1;
                    }
                    res[i3] = t - huv[i3]; huv[i3] = t;
                }
            } else if (closure==CLOSURE_PSE4){
                for (size_t i3=i3begin; i3<i3end; i3++){
                    t = -uuv[i3] + huv[i3] - cuv[i3];
                    if (t>0){
                        t = t + t*t/2 + t*t*t/6 + t*t*t*t/24;
                    } else {
                        t = exp(t) - 1;
                    }
                    res[i3] = t - huv[i3]; huv[i3] = t;
                }
            } else if (closure==CLOSURE_PSE5){
                for (size_t i3=i3begin; i3<i3end; i3++){
                    t = -uuv[i3] + huv[i3] - cuv[i3];
                    if (t>0){
                        t = t + t*t/2 + t*t*t/6 + t*t*t*t/24 + t*t*t*t*t/120;
                    } else {
                        t = exp(t) - 1;
                    }
                    res[i3] = t - huv[i3]; huv[i3] = t;
                }
            } else if (closure==CLOSURE_PSE6){
                for (size_t i3=i3begin; i3<i3end; i3++){
                    t = -uuv[i3] + huv[i3] - cuv[i3];
                    if (t>0){
                        t = t + t*t/2 + t*t*t/6 + t*t*t*t/24 + t*t*t*t*t/120 + t*t*t*t*t*t/720;
                    } else {
                        t = exp(t) - 1;
                    }
                    res[i3] = t - huv[i3]; huv[i3] = t;
                }
            } else if (closure==CLOSURE_PSE7){
                for (size_t i3=i3begin; i3<i3end; i3++){
                    t = -uuv[i3] + huv[i3] - cuv[i3];
                    if (t>0){
                        t = t + t*t/2 + t*t*t/6 + t*t*t*t/24 + t*t*t*t*t/120 + t*t*t*t*t*t/720 + t*t*t*t*t*t*t/5040;
                    } else {
                        t = exp(t) - 1;
                    }
                    res[i3] = t - huv[i3]; huv[i3] = t;
                }
            } else if (closure==CLOSURE_PSE8){
                for (size_t i3=i3begin; i3<i3end; i3++){
                    t = -uuv[i3] + huv[i3] - cuv[i3];
                    if (t>0){
                        t = t + t*t/2 + t*t*t/6 + t*t*t*t/24 + t*t*t*t*t/120 + t*t*t*t*t*t/720 + t*t*t*t*t*t*t/5040 + t*t*t*t*t*t*t*t/40320;
                    } else {
                        t = exp(t) - 1;
                    }
                    res[i3] = t - huv[i3]; huv[i3] = t;
                }
            } else if (closure==CLOSURE_PSE9){
                for (size_t i3=i3begin; i3<i3end; i3++){
                    t = -uuv[i3] + huv[i3] - cuv[i3];
                    if (t>0){
                        t = t + t*t/2 + t*t*t/6 + t*t*t*t/24 + t*t*t*t*t/120 + t*t*t*t*t*t/720 + t*t*t*t*t*t*t/5040 + t*t*t*t*t*t*t*t/40320 + t*t*t*t*t*t*t*t*t/362880;
                    } else {
                        t = exp(t) - 1;
                    }
                    res[i3] = t - huv[i3]; huv[i3] = t;
                }
            } else if (closure==CLOSURE_PSE10){
                for (size_t i3=i3begin; i3<i3end; i3++){
                    t = -uuv[i3] + huv[i3] - cuv[i3];
                    if (t>0){
                        t = t + t*t/2 + t*t*t/6 + t*t*t*t/24 + t*t*t*t*t/120 + t*t*t*t*t*t/720 + t*t*t*t*t*t*t/5040 + t*t*t*t*t*t*t*t/40320 + t*t*t*t*t*t*t*t*t/362880 + t*t*t*t*t*t*t*t*t*t/3628800;
                    } else {
                        t = exp(t) - 1;
                    }
                    res[i3] = t - huv[i3]; huv[i3] = t;
                }
            } else if (closure==CLOSURE_PSE){
                for (size_t i3=i3begin; i3<i3end; i3++){
                    t = -uuv[i3] + huv[i3] - cuv[i3];
                    if (t>0){ double s = t;
                        t += ((factor-1)>1?1:(factor-1)<0?0:(factor-1)) * s*s/2;
                        t += ((factor-2)>1?1:(factor-2)<0?0:(factor-2)) * s*s*s/6;
                        t += ((factor-3)>1?1:(factor-3)<0?0:(factor-3)) * s*s*s*s/24;
                        t += ((factor-4)>1?1:(factor-4)<0?0:(factor-4)) * s*s*s*s*s/120;
                        t += ((factor-5)>1?1:(factor-5)<0?0:(factor-5)) * s*s*s*s*s*s/720;
                        t += ((factor-6)>1?1:(factor-6)<0?0:(factor-6)) * s*s*s*s*s*s*s/5040;
                        t += ((factor-7)>1?1:(factor-7)<0?0:(factor-7)) * s*s*s*s*s*s*s*s/40320;
                        t += ((factor-8)>1?1:(factor-8)<0?0:(factor-8)) * s*s*s*s*s*s*s*s*s/362880;
                        t += ((factor-9)>1?1:(factor-9)<0?0:(factor-9)) * s*s*s*s*s*s*s*s*s*s/3628800;
                    } else {
                        t = exp(t) - 1;
                    }
                    res[i3] = t - huv[i3]; huv[i3] = t;
                }
            } else if (closure==CLOSURE_PY){
                for (size_t i3=i3begin; i3<i3end; i3++){

                    t = exp(-uuv[i3])*(1 + huv[i3] - cuv[i3]) - 1;
                    res[i3] = t - huv[i3]; huv[i3] = t;
                }
            } else if (closure==CLOSURE_HNCB){
                for (size_t i3=i3begin; i3<i3end; i3++){

                    t = huv[i3] - cuv[i3];
                    double hh = extra? extra[i3] : 0;
                    t = exp(-uuv[i3] + t - hh*hh/2) - 1;
                    res[i3] = t - huv[i3]; huv[i3] = t;
                }
            } else if (closure==CLOSURE_D2){
                for (size_t i3=i3begin; i3<i3end; i3++){

                    t = huv[i3] - cuv[i3];
                    t = exp(-uuv[i3] + t - t*t/2) - 1;
                    res[i3] = t - huv[i3]; huv[i3] = t;
                }
            } else if (closure==CLOSURE_MS){
                for (size_t i3=i3begin; i3<i3end; i3++){

                    t = huv[i3] - cuv[i3];
                    t = exp(-uuv[i3] + sqrt((1+2*t)>0?(1+2*t):0)-1) - 1;
                    res[i3] = t - huv[i3]; huv[i3] = t;
                }
            } else if (closure==CLOSURE_MSHNC){
                for (size_t i3=i3begin; i3<i3end; i3++){

                    t = huv[i3] - cuv[i3];
                    if (t<=0){
                        t = exp(-uuv[i3] + t) - 1;
                        res[i3] = t - huv[i3]; huv[i3] = t;
                    } else {
                        t = exp(-uuv[i3] + sqrt(1+2*t)-1) - 1;
                        res[i3] = t - huv[i3]; huv[i3] = t;
                    }
                }
            } else if (closure==CLOSURE_BPGG){
                for (size_t i3=i3begin; i3<i3end; i3++){

                    t = huv[i3] - cuv[i3];
                    t = exp(-uuv[i3] + pow((1+factor*t)>0?(1+factor*t):0, 1.0/factor)-1) - 1;
                    res[i3] = t - huv[i3]; huv[i3] = t;
                }
            } else if (closure==CLOSURE_BPGGHNC){
                for (size_t i3=i3begin; i3<i3end; i3++){

                    t = huv[i3] - cuv[i3];
                    if (t<=0){
                        t = exp(-uuv[i3] + t) - 1;
                        res[i3] = t - huv[i3]; huv[i3] = t;
                    } else {
                        t = exp(-uuv[i3] + pow(1+factor*t, 1.0/factor)-1) - 1;
                        res[i3] = t - huv[i3]; huv[i3] = t;
                    }
                }
            } else if (closure==CLOSURE_VM){
                for (size_t i3=i3begin; i3<i3end; i3++){

                    t = huv[i3] - cuv[i3];
                    t = exp(-uuv[i3] + t - t*t/2/(1+factor*t)) - 1;
                    res[i3] = t - huv[i3]; huv[i3] = t;
                }
            } else if (closure==CLOSURE_MHNC){
                for (size_t i3=i3begin; i3<i3end; i3++){

                    t = huv[i3] - cuv[i3];
                    double E_mhnc = - t*t/(2 + 2*0.8*t);
                    t = exp(-uuv[i3] + t + E_mhnc) - 1;
                    res[i3] = t - huv[i3]; huv[i3] = t;
                }
            } else if (closure==CLOSURE_MP){
                for (size_t i3=i3begin; i3<i3end; i3++){
                    t = huv[i3] - cuv[i3];
                    t = exp(-uuv[i3])*((1+factor)*exp(t/(1+factor)) - factor) - 1;
                    res[i3] = t - huv[i3]; huv[i3] = t;
                }
            } else if (closure==CLOSURE_RBC_HNC){
                for (size_t i3=i3begin; i3<i3end; i3++){
                    t = huv[i3] - cuv[i3];
                    t = exp(-uuv[i3] + factor*t) * (extra? extra[i3] : 1) - 1;
                    res[i3] = t - huv[i3]; huv[i3] = t;
                }
            } else if (closure==CLOSURE_RBC_KH){
                for (size_t i3=i3begin; i3<i3end; i3++){
                    t = -uuv[i3] + huv[i3] - cuv[i3];
                    t = (t>0? t+1 +ln(extra?fabs(extra[i3]+MACHINE_REASONABLE_ERROR):1) : exp(t)*(extra?extra[i3]:1)) - 1;
                    res[i3] = t - huv[i3]; huv[i3] = t;
                }
            } else if (closure==CLOSURE_USER1){
                for (size_t i3=i3begin; i3<i3end; i3++){
                    t = exp(-uuv[i3] + huv[i3] - cuv[i3]) - 1;
                    res[i3] = t - huv[i3]; huv[i3] = t;
                }
            } else if (closure==CLOSURE_USER2){
                for (size_t i3=i3begin; i3<i3end; i3++){
                    t = exp(-uuv[i3] + huv[i3] - cuv[i3]) - 1;
                    res[i3] = t - huv[i3]; huv[i3] = t;
                }
            } else if (closure==CLOSURE_USER3){
                for (size_t i3=i3begin; i3<i3end; i3++){
                    t = exp(-uuv[i3] + huv[i3] - cuv[i3]) - 1;
                    res[i3] = t - huv[i3]; huv[i3] = t;
                }
            } else { // default: HNC (no scaling) or experimental
                for (size_t i3=i3begin; i3<i3end; i3++){
                    t = (huv[i3] - cuv[i3]);
                    t = exp(-uuv[i3] + t) - 1;
                    res[i3] = t - huv[i3]; huv[i3] = t;
                }
            }

            for (size_t i3=i3begin; i3<i3end; i3++) ret += res[i3]*res[i3];

          // closure enhancement
            if (sys->closure_enhance_level==2){
                for (size_t i3=i3begin; i3<i3end; i3++){
                    //double g = 1 + huv[i3]; res[i3] /= sqrt(1 + g*g);
                    //res[i3] /= sqrt(1 + huv[i3]*huv[i3]);
                    res[i3] /= (1 + huv[i3]*huv[i3]);
                }
            } else if (sys->closure_enhance_level==1){
                for (size_t i3=i3begin; i3<i3end; i3++){
                    res[i3] /= sqrt(1 + huv[i3]*huv[i3]);
                }
            } else if (sys->closure_enhance_level>MACHINE_REASONABLE_ERROR){
                for (size_t i3=i3begin; i3<i3end; i3++){
                    res[i3] /= pow(1 + huv[i3]*huv[i3], sys->closure_enhance_level/2);
                }
            }
          // closure cutoff
            if (sys->closure_enhance_cutoff>0){
                for (size_t i3=i3begin; i3<i3end; i3++){
                    if (res[i3]>sys->closure_enhance_cutoff){
                        res[i3] = sys->closure_enhance_cutoff;
                    } else if (res[i3]<-sys->closure_enhance_cutoff){
                        res[i3] = -sys->closure_enhance_cutoff;
                    }
                }
            }
        }
        return ret;
    }

    double perform_closure(IET_Param * sys, IET_arrays * arr, int id){
        size_t N3 = arr->nx*arr->ny*arr->nz;
      #ifdef _LOCALPARALLEL_
        size_t idjamin = N3 / sys->nt * id; size_t idjamax = N3 / sys->nt * (id+1); if (id+1 >= sys->nt) idjamax = N3;
      #else
        size_t idjamin = 0; size_t idjamax = N3;
      #endif
        return perform_closure(sys, N3, idjamin, idjamax, arr->uuv, arr->ulr, arr->cuv, arr->huv, arr->rismhi_cache[2], arr->res, arr->hlr, arr->clr, arr->dd);
    }
    double perform_closure(IET_Param * sys, IET_arrays * arr){
      #ifdef _LOCALPARALLEL_
        for (int i=1; i<sys->nt; i++) __mp_tasks[i] = MPTASK_RISM_CLOSURE;
        __mp_returns[0] = perform_closure(sys, arr, 0);
        wait_subroutines(sys);
        double ret = 0; for (int i=0; i<sys->nt; i++) ret += __mp_returns[i];
        return ret;
      #else
        size_t N3 = arr->nx*arr->ny*arr->nz;
        return perform_closure(sys, N3, 0, N3, arr->uuv, arr->ulr, arr->cuv, arr->huv, arr->rismhi_cache[2], arr->res, arr->hlr, arr->clr, arr->dd);
      #endif
    }


    //>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    //>>>>>>>                            >>>>>>>>>>>>>>>>>>>>>>>>>>
    //>>>>>>>     Original RISM Main     >>>>>>>>>>>>>>>>>>>>>>>>>>
    //>>>>>>>                            >>>>>>>>>>>>>>>>>>>>>>>>>>
    //>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    void perform_rism_equation(IET_Param * sys, IET_arrays * arr, __REAL__ **** cuv, __REAL__ **** huv, __REAL__ **** res, __REAL__ **** dd, int nx, int ny, int nz, int nv, Vector box, __REAL__ *** wvv, double dk_wvv, int n_wvv, __REAL__ *** nhkvv, double dk_nhkvv, int n_nhkvv, __REAL__ *** xvv, double dk_xvv, int n_xvv, double *** fftin, double *** fftout, fftw_plan & planf, fftw_plan & planb){
            size_t N3 = nx * ny * nz; size_t N4 = N3 * nv;
            if (xvv && (!dd||sys->use_homogeneous_rism)){
                perform_3rx1k_convolution(&arr->fftw_mp, cuv, nx, ny, nz, box, nv, nv, xvv, dk_wvv, sys->xvv_k_shift, n_xvv, huv, fftin, fftout, planf, planb, true);
            } else {
              // 1. C(r) * w(r) -> huv
                if (wvv){
                    perform_3rx1k_convolution(&arr->fftw_mp, cuv, nx, ny, nz, box, nv, nv, wvv, dk_wvv, sys->xvv_k_shift, n_wvv, huv, fftin, fftout, planf, planb, true);
                }
              // 2. C(r)n(r) -> res, C(r)n(r) * (w(r) + nhkvv(r)) -> huv
                if (nhkvv){
                    if (sys->use_homogeneous_rism){
                        for (int iv=0; iv<nv; iv++){
                            int ivm = sys->av[iv].iaa; //fprintf(stderr, "%d -(mol)-> %d : %s\n", iv, ivm, sys->av[iv].name);
                            for (size_t i3=0; i3<N3; i3++) res[iv][0][0][i3] = cuv[iv][0][0][i3] * sys->nbulk[ivm];
                        }
                        perform_3rx1k_convolution(&arr->fftw_mp, res, nx, ny, nz, box, nv, nv, nhkvv, dk_nhkvv, sys->xvv_k_shift, n_nhkvv, huv, fftin, fftout, planf, planb, false);
                        //perform_3rx1k_convolution(&arr->fftw_mp, cuv, nx, ny, nz, box, nv, nv, nhkvv, dk_nhkvv, sys->xvv_k_shift, n_nhkvv, huv, fftin, fftout, planf, planb, wvv?false:true);
                    } else {
                        for (int iv=0; iv<nv; iv++){
                            int ivm = sys->av[iv].iaa; //fprintf(stderr, "%d -(mol)-> %d : %s\n", iv, ivm, sys->av[iv].name);
                            if (dd){
                                for (size_t i3=0; i3<N3; i3++) res[iv][0][0][i3] = cuv[iv][0][0][i3] * dd[ivm][0][0][i3];
                            } else {
                                for (size_t i3=0; i3<N3; i3++) res[iv][0][0][i3] = cuv[iv][0][0][i3] * sys->nbulk[ivm];
                            }
                        }
                        perform_3rx1k_convolution(&arr->fftw_mp, res, nx, ny, nz, box, nv, nv, nhkvv, dk_nhkvv, sys->xvv_k_shift, n_nhkvv, huv, fftin, fftout, planf, planb, false);
                    }
                }
            }
            lap_timer_fftw();
    }
    void perform_rism_equation(IET_Param * sys, IET_arrays * arr, __REAL__ **** cuv, __REAL__ **** huv, __REAL__ **** res, __REAL__ *** wvv=nullptr, __REAL__ *** nhkvv=nullptr, __REAL__ *** xvv=nullptr){
        if (!wvv) wvv = arr->convolution_wvv; if (!nhkvv) nhkvv = arr->convolution_nhkvv;
        if (!xvv) xvv = arr->xvv;
        perform_rism_equation(sys, arr, cuv, huv, res, arr->dd, arr->nx, arr->ny, arr->nz, arr->nv, arr->box, wvv, arr->dk_wvv, arr->n_wvv, nhkvv, arr->dk_nhkvv, arr->n_nhkvv, xvv, arr->dk_wvv, arr->n_xvv, arr->fftin, arr->fftout, arr->planf, arr->planb);
    }
    void perform_rism_equation_without_hi(IET_Param * sys, IET_arrays * arr, __REAL__ **** cuv, __REAL__ **** huv, __REAL__ **** res, __REAL__ *** wvv=nullptr, __REAL__ *** nhkvv=nullptr, __REAL__ *** xvv=nullptr){
        if (!wvv) wvv = arr->convolution_wvv; if (!nhkvv) nhkvv = arr->convolution_nhkvv;
        if (!xvv) xvv = arr->xvv;
        perform_rism_equation(sys, arr, cuv, huv, res, nullptr, arr->nx, arr->ny, arr->nz, arr->nv, arr->box, wvv, arr->dk_wvv, arr->n_wvv, nhkvv, arr->dk_nhkvv, arr->n_nhkvv, xvv, arr->dk_wvv, arr->n_xvv, arr->fftin, arr->fftout, arr->planf, arr->planb);
    }
    //>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    //>>>>>>>                              >>>>>>>>>>>>>>>>>>>>>>>>
    //>>>>>>>    Closure preparation       >>>>>>>>>>>>>>>>>>>>>>>>
    //>>>>>>>                              >>>>>>>>>>>>>>>>>>>>>>>>
    //>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    void prepare_for_closure(IET_Param * sys, IET_arrays * arr){
        size_t N3 = arr->nx * arr->ny * arr->nz; int N4 = N3 * sys->nv;
        double beta = sys->default_temperature / sys->temperature;

        bool enable_HNCB = false;       // HNCB: extra = h*hbar
        bool enable_RBC = false;        // RBC: extra = exp(-ljr) * wvv
        bool enable_prepare = false;


        for (int i=0; i<sys->nv; i++){
            if (sys->closures[i] == CLOSURE_HNCB) enable_HNCB = enable_prepare = true;
            else if (sys->closures[i] == CLOSURE_RBC_HNC) enable_RBC = enable_prepare = true;
            else if (sys->closures[i] == CLOSURE_RBC_KH) enable_RBC = enable_prepare = true;
            //else if (sys->closures[i] == CLOSURE_RCCA2) prepare_cmds[n_prepare++] = CLOSURE_RCCA2;
        }
        if (enable_prepare){
            __REAL__ **** extra = arr->rismhi_cache[2];
            __REAL__ **** temp = arr->rismhi_cache[3];
            clear_tensor4d(extra, N4);
            if (enable_HNCB){
                perform_rism_equation(sys, arr, arr->huv, temp, arr->res, arr->convolution_wvv, arr->convolution_nhkvv);
                for (int iv=0; iv<sys->nv; iv++) if (sys->closures[iv]==CLOSURE_HNCB){
                    cp_tensor3d(temp[iv], extra[iv], N3);
                }
            } else if (enable_RBC && arr->uljr){
                for (size_t i4=0; i4<N4; i4++) extra[0][0][0][i4] = exp(-arr->uljr[0][0][0][i4]*beta);
                perform_rism_equation(sys, arr, extra, temp, arr->res, arr->convolution_wvv, nullptr);
                clear_tensor4d(extra, N4);
                for (int iv=0; iv<sys->nv; iv++) if (sys->closures[iv]==CLOSURE_RBC_HNC || sys->closures[iv]==CLOSURE_RBC_KH){
                    cp_tensor3d(temp[iv], extra[iv], N3);
                }
            }

            /*
            for (int ip=0; ip<n_prepare; ip++){
                if (prepare_cmds[ip] == CLOSURE_HNCB){
                    perform_rism_equation(sys, arr, arr->huv, temp, arr->res, arr->convolution_wvv, arr->convolution_nhkvv);
                }
                for (int iv=0; iv<sys->nv; iv++) if (sys->closures[iv]==prepare_cmds[ip]){
                    cp_tensor3d(temp[iv], extra[iv], N3);
                }
            }
            */
        }
    }

    //>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    //>>>>>>>                              >>>>>>>>>>>>>>>>>>>>>>>>
    //>>>>>>>    Unrenormalized OZ Main    >>>>>>>>>>>>>>>>>>>>>>>>
    //>>>>>>>                              >>>>>>>>>>>>>>>>>>>>>>>>
    //>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    double perform_3drism_loop_ssoz(IET_Param * sys, IET_arrays * arr){
      // 1. RISM equation
        perform_rism_equation(sys, arr, arr->cuv, arr->huv, arr->res, arr->convolution_wvv, arr->convolution_nhkvv);
      // 2. closure: ΔH -> res
        prepare_for_closure(sys, arr);
        double ret = perform_closure(sys, arr);
      // x. end
        lap_timer_rism();
        return ret;
    }

    //>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    //>>>>>>>                              >>>>>>>>>>>>>>>>>>>>>>>>
    //>>>>>>>    Renormalized RISM Main    >>>>>>>>>>>>>>>>>>>>>>>>
    //>>>>>>>                              >>>>>>>>>>>>>>>>>>>>>>>>
    //>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    double perform_3drism_loop_rrism(IET_Param * sys, IET_arrays * arr){
        size_t N3 = arr->nx * arr->ny * arr->nz; size_t N4 = N3 * sys->nv; double beta = sys->default_temperature / sys->temperature;
      // 1. hsr -> arr->huv
        perform_rism_equation(sys, arr, arr->cuv, arr->huv, arr->res, arr->convolution_wvv, arr->convolution_nhkvv);
      // 2. h = hsr + hlr
      /*
        for (size_t i4=0; i4<N4; i4++){
            arr->huv[0][0][0][i4] += arr->hlr[0][0][0][i4];
        }
      */
        for (int iv=0; iv<sys->nv; iv++){
            for (size_t i3=0; i3<N3; i3++){
              // note: enable_dilute_rism only perform on closure, so here huv is similar to disable_dilute_rism
                arr->huv[iv][0][0][i3] = arr->huv[iv][0][0][i3] + arr->hlr[iv][0][0][i3];
            }
        }
      // 3. sr closure
        prepare_for_closure(sys, arr);
        double ret = perform_closure(sys, arr);
        lap_timer_rism();
        return ret;
    }

    //>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    //>>>>>>>                                  >>>>>>>>>>>>>>>>>>>>
    //>>>>>>>  General Renormalization Scheme  >>>>>>>>>>>>>>>>>>>>
    //>>>>>>>                                  >>>>>>>>>>>>>>>>>>>>
    //>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    void prepare_3d_rism(IET_Param * sys, IET_arrays * arr){
      // prepare for the short-long range separation
        if (sys->ietal==IETAL_SSOZ){
        } else {
          // long range total correlation: arr->hlr = (β U_{CoulLR+SR} q / εr) ρv * χvv
            size_t N3 = arr->nx * arr->ny * arr->nz; size_t N4 = N3 * sys->nv;
            double beta = sys->default_temperature / sys->temperature;
          // 1. arr->clr = clr = - ulr
            for (size_t i4=0; i4<N4; i4++) arr->clr[0][0][0][i4] = - arr->ulr[0][0][0][i4];
          // 2. arr->hlr = clr ρv * χvv
            __REAL__ *** xvv = arr->xvv;
            __REAL__ *** wvv = arr->wvv;
            __REAL__ *** nhkvv = arr->nhkvv;
            if (sys->hlr_no_hi){
                perform_rism_equation_without_hi(sys, arr, arr->clr, arr->hlr, arr->res, wvv, nhkvv, xvv);
            } else {
                perform_rism_equation(sys, arr, arr->clr, arr->hlr, arr->res, wvv, nhkvv, xvv);
            }
        }
    }

    //>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    //>>>>>>>                                  >>>>>>>>>>>>>>>>>>>>
    //>>>>>>>         report functions         >>>>>>>>>>>>>>>>>>>>
    //>>>>>>>                                  >>>>>>>>>>>>>>>>>>>>
    //>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    void fill_closure_additional_buffer(IET_Param * sys, int ic, char * additional_buffer, size_t sizeof_additional_buffer){
        int closure= sys->closures[ic];
        double closure_factor = sys->closure_factors[ic];
        if (closure==CLOSURE_PLHNC){
            snprintf(additional_buffer, sizeof_additional_buffer, "(%g)", closure_factor);
        } else {
            bool display_additional = false;
            for (int i=0; !display_additional&&i<sizeof(closures_needto_display_factor)/sizeof(closures_needto_display_factor[0]); i++) if (closure==closures_needto_display_factor[i]) display_additional = true;
            if (display_additional){
                snprintf(additional_buffer, sizeof_additional_buffer, "(%g)", closure_factor);
            } else strncpy(additional_buffer, "", sizeof_additional_buffer);
        }
    }
    bool fill_closure_buffer(IET_Param * sys, char * closure_buffer, size_t sizeof_closure_buffer){
        memset(closure_buffer, 0, sizeof_closure_buffer);
        bool multiple_closure = false; for (int iv=1; iv<sys->nv; iv++) if (sys->closures[iv] != sys->closures[0]) multiple_closure = true;

        char additional_buffer[64]; fill_closure_additional_buffer(sys, 0, additional_buffer, sizeof(additional_buffer));
        if (multiple_closure){
            int closure_fill[MAX_SOL]; int n_closure_fill[MAX_SOL]; int n_fill = 0;
            for (int iv=0; iv<sys->nv; iv++){
                if (n_fill<=0||sys->closures[iv]!=closure_fill[n_fill-1]){
                    n_fill++; closure_fill[n_fill-1] = sys->closures[iv]; n_closure_fill[n_fill-1] = 1;
                } else {
                    n_closure_fill[n_fill-1]++;
                }
            }
            for (int iv=0; iv<n_fill; iv++){
                if (iv<4){
                    char buffer[128]; memset(buffer, 0, sizeof(buffer));
                    fill_closure_additional_buffer(sys, iv, additional_buffer, sizeof(additional_buffer));
                    if (n_closure_fill[iv]==1){
                        snprintf(buffer, sizeof(buffer), "%s%s%s", iv>0?",":"", CLOSURE_name[closure_fill[iv]], additional_buffer);
                    } else {
                        snprintf(buffer, sizeof(buffer), "%s%sx%d%s", iv>0?",":"", CLOSURE_name[closure_fill[iv]], n_closure_fill[iv], additional_buffer);
                    }
                    strncat(closure_buffer, buffer, sizeof_closure_buffer-strlen(buffer));
                } else if (iv==4) strncat(closure_buffer, ",...", sizeof_closure_buffer-4);
            }

            /*for (int iv=0; iv<sys->nv; iv++){
                if (iv<4){
                    if (iv>0) strcat(closure_buffer, ","); strcat(closure_buffer, CLOSURE_name[sys->closures[iv]]);
                } else if (iv==4) strcat(closure_buffer, ",...");
            }*/
        } else snprintf(closure_buffer, sizeof_closure_buffer, "%s%s", CLOSURE_name[sys->closures[0]], additional_buffer); // strncpy(closure_buffer, CLOSURE_name[sys->closures[0]], sizeof_closure_buffer);
        return multiple_closure;
    }

    bool report_3drism_convergence(FILE * flog, IET_Param * sys, IET_arrays * arr, double err, bool & success, const char * closure_buffer, int istep, bool is_scf=true){
        bool is_out_tty = isatty(fileno(flog));
      // 3. report and loop control
        bool stop_loop = false;
            if (err <= sys->errtolrism){ success = true; stop_loop = true; }
            if (err < 0 || err > 1e7){ success = false; stop_loop = true; } // error diverge
            //if (diverge_counter>10){ success = false; stop_loop = true; } // stop converging
            if (istep+1 >= sys->stepmax_rism) stop_loop = true;
        //if (flog && (sys->detail_level>=1 || istep+1>=sys->stepmax_rism || stop_loop)){
        if (flog && (sys->detail_level>=1 || stop_loop)){
            fprintf(flog, "  %s-%s step %d %s ", IETAL_name[sys->ietal], closure_buffer, istep+1, is_scf?"stdev":"change");
            //fprintf(flog, (err<sys->errtolrism||istep+1>=sys->stepmax_rism)? "%g" : err<1e-3?"%.1e":"%.4f", err); //fprintf(flog, (" %.2g"), err);
            fprintf(flog, err>0.1?"%.5g":"%.5e", err);
            if (sys->ndiis_rism>1 && arr->diis_rism.ndiis>0) fprintf(flog, " (DIISx%d)", arr->diis_rism.ndiis);
            if (sys->debug_level>=3||sys->detail_level>=3){ char time_buffer[64]; fprintf(flog, " (%s)", get_current_time_text(time_buffer)); }
            fprintf(flog, "  %s", stop_loop? "\n" : (istep<sys->stepmax_rism&&is_out_tty&&sys->detail_level==1)?"  \r":"  \n");
            //fflush(flog);
        }
        fflush(flog);
        //if (!stop_loop && sys->detail_level==1 && flog && flog!=stderr && flog!=stdout) fseek(flog, flog_pointer, SEEK_SET);
        return stop_loop;
    }

    //>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    //>>>>>>>                             >>>>>>>>>>>>>>>>>>>>>>>>>
    //>>>>>>>      Main RISM Control      >>>>>>>>>>>>>>>>>>>>>>>>>
    //>>>>>>>                             >>>>>>>>>>>>>>>>>>>>>>>>>
    //>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    bool perform_3drism(FILE * flog, IET_Param * sys, IET_arrays * arr, bool guess_init_n=true){
        bool success = false; double err = 1; double err_last = 1; int diverge_counter = 0;
        size_t N3 = arr->nx*arr->ny*arr->nz; size_t N4 = sys->nv * N3;

        char closure_buffer[128];
        if (fill_closure_buffer(sys, closure_buffer, sizeof(closure_buffer))){
            if (sys->debug_level>=1){
                fprintf(flog, "debug:: multiple closures:");
                for (int iv=0; iv<sys->nv; iv++){
                    fprintf(flog, " %d/%s.%s:%s", iv+1, sys->av[iv].mole, sys->av[iv].name, CLOSURE_name[sys->closures[iv]]);
                }
                fprintf(flog, "\n");
            }
        }

      // step 1: other preparations
        prepare_3d_rism(sys, arr);
      // step 2: RISM loop
        //size_t flog_pointer = 0;
        for (int istep=0; istep<sys->stepmax_rism; istep++){
            arr->istep = istep;
            while (sys->suspend_calculation){ usleep(100); sys->is_suspend_calculation = true; continue; }
            sys->is_suspend_calculation = false;
            //if (flog && flog!=stderr && flog!=stdout) flog_pointer = ftell(flog);
          // 1. Perform RISM calculation
            double err_original = 0;
            if (sys->ietal==IETAL_SSOZ){            // SSOZ
                err_original = perform_3drism_loop_ssoz(sys, arr);
            } else if (sys->ietal==IETAL_RRISM){    // RRISM
                err_original = perform_3drism_loop_rrism(sys, arr);
            } else {
                return false;
            }
            err_original = sqrt(err_original/N4);
            lap_timer_rism();
          // 2. DIIS step in
            if (sys->_n_hold_list>0){
                for (int i=0; i<sys->_n_hold_list; i++){
                    int ivh = sys->_hold_list[i]-1; if (ivh>=0 && ivh<sys->nv) clear_tensor3d(arr->res[ivh], N3);
                }
            }
            if (sys->ndiis_rism<=1){
                for (size_t i4=0; i4<N4; i4++){
                    double res = arr->res[0][0][0][i4]; err += res * res;
                    arr->cuv[0][0][0][i4] += res * sys->delrism;
                }
                err = sqrt(fabs(err)/N4);
            } else {
                arr->diis_current = &arr->diis_rism;
                err = arr->diis_rism.advance(sys, &arr->cuv[0][0][0][0], &arr->res[0][0][0][0], sys->delrism, true);
            }
            if (err >= err_last*(1-sqrt(MACHINE_REASONABLE_ERROR))) diverge_counter ++; err_last = err;
            lap_timer_diis();
          // 3. report and loop control
            bool stop_loop = report_3drism_convergence(flog, sys, arr, sys->err_from_unscaled_res? err_original : err, success, closure_buffer, istep);
            if (stop_loop) break;
        }

        return success;
    }

}

bool perform_3d_iet(FILE * flog, IET_Param * sys, IET_arrays * arr, bool guess_init_n=true){
    if (sys->ietal==IETAL_SSOZ){
        if (sys->debug_level>=2) fprintf(sys->log(), "DEBUG:: perform_3drism(algorithm=SSOZ)\n");
        return RISMHI3D_RISMNS::perform_3drism(flog, sys, arr, guess_init_n);
    } else if (sys->ietal==IETAL_RRISM){
        if (sys->debug_level>=2) fprintf(sys->log(), "DEBUG:: perform_3drism(algorithm=RISM)\n");
        return RISMHI3D_RISMNS::perform_3drism(flog, sys, arr, guess_init_n);
    } else {
        return false;
    }
}
