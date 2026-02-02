#include "Utils.h"
#include "hungarian.h"
#include "matching.h"

#include "ODESolverCVODE.h"
#include "ODESolverDP5.h"

#include<iostream>
#include<cmath>
#include<iomanip>
#include<complex>
#include<limits>
#include<limits>
#include<cstdlib>
#include<fstream>
#include<algorithm>
#include<chrono>
#include<assert.h>
using namespace std;


//Prototipi 
//INTEGRATOR 
//void dYdt(double ,double *,double *);
void dYdt_reduced(double, const double*, double*, double*, double*);
//void RK4Step(double t,double *Y, double h,int Neq,void(*RHS_Func)(double,double*,double*))
//example
//void RK4Step(double ,double *, double ,int,void(*)(double,double*,double*));
void RK4StepWithC(double, double*, double, int, void(*)(double, double*, double*, double*, double*), double*, double*);
//REDUCED Competitive Inhibitor 
// === FUNCTION PROTOTYPES ===
double cubic_root_clamped_scalar(double a3, double a2, double a1, double a0,
                                 double lo, double hi, double x0,
                                 double ftol, double xtol, int maxit);

double solve_quadratic_in_interval(double a2, double a1, double a0, double lo, double hi);

void phy_coeffs(double k1, double k2, double E, double S1, double S2, double coeffs[4]);

double select_root_scalar(double k1, double k2, double E, double S1, double S2, 
                         double x_prev, double rtol, double atol);

struct Solution {
    double C1, C2, C3;
    bool success;
};

// === FUNCTION PROTOTYPES ===
bool solveLinear3x3(double A[3][3], double b[3], double x[3]);

void computeEquations(double C1, double C2, double C3,
                      double E_T, double I_hat, double S_hat,
                      double a1, double a2, double i, double k1, double d1, double d2, double i_rev,
                      double F[3]);

void computeJacobian(double C1, double C2, double C3,
                     double E_T, double I_hat, double S_hat,
                     double a1, double a2, double i, double k1, double d1, double d2, double i_rev,
                     double J[3][3]);

Solution solve_Cs_allosteric(double E_T, double I_hat, double S_hat,
                             double a1, double a2, double i, double k1,
                             double d1, double d2, double i_rev,
                             bool init_p, double init_v[3], bool log_i);

double C_root(double, double, double);
struct UserData {
    int rhs_evals;
    Real C_evaluated_at_t;
    double* C;
    double* params;
    ofstream* pfdata;
};
bool derivative(OdeReal t, const OdeReal* y, OdeReal* ydot, void* user_data)
{
    UserData* ud = reinterpret_cast<UserData*>(user_data);
    dYdt_reduced(t, y, ydot, ud->C, ud->params);
    ud->C_evaluated_at_t = t;
    ud->rhs_evals++;
    return true;
}
void integration_step_cb(OdeReal t, const OdeReal* y, void* user_data)
{
    UserData* ud = reinterpret_cast<UserData*>(user_data);
    (*ud->pfdata) << t;
    for (int i = 0; i < 4; i++) {
        (*ud->pfdata) << "\t" << y[i];
    }
    for (int i = 0; i < 13; i++) {
        (*ud->pfdata) << "\t" << ud->C[i];
    }
    (*ud->pfdata) << endl;
}


int main(){
    // we do not need to define 
    //double Y[9]; // Y[8]
    //storeC previous iteration for numerical solver 
    //initialize arrays
    OdeVectorReal Y(4);
    double C[14] = { 0.0, 0.0, 0.0,0.0, 0.0, 0.0,0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,0.0 };
    double params[34];
    //parameters
    double step;



    double k22a = 0.106, km22a = 0.02385, k20 = 0.66;  // E+S ↔ ES → E+pS
    double k22b = 0.106, km22b = 0.0159, k220 = 0.0159;
    double k1a = 0.106, km1a = 0.02385, k2a = 0.66;   // E+S ↔ ES → E+pS
    double k1b = 0.106, km1b = 0.0159, k2b = 0.0159;  //  pS+p1 ↔ pSp1 → S+p1
    double k1c = 0.106, km1c = 0.01184, k2c = 0.0242; // pS+S1 ↔ pSS1 → pS+pS1
    double k1d = 0.106, km1d = 0.0159, k2d = 0.0159;   //  pS1+p2 ↔ pS1p2 → S1+p2
    double k1e = 0.106, km1e = 0.000421, k2e = 0.000639; // pS1+S2 ↔ pS1S2 → pS1+pS2
    double k1f = 0.106, km1f = 0.0159, k2f = 0.0159;     //  pS2+p3 ↔ pS2p3 → S2+p3

    double ki_EGF = 0.106, kminusi_EGF = 0.00005936;
    double ki_S = 0.106, kminusi_S = 0.00005936;
    double ki_S1 = 13.06, kminusi_S1 = 0.0012296; //changed to scale trametinib on MLD scale 
    double ki_S2 = 0.106, kminusi_S2 = 0.0002967;

    double KMegf = (km22a + k2a) / k22a;
    double KMa = (km1a + k2a) / k1a;
    double KMb = (km1b + k2b) / k1b;
    double KMc = (km1c + k2c) / k1c;
    double KMd = (km1d + k2d) / k1d;
    double KMe = (km1e + k2e) / k1e;
    double KMf = (km1f + k2f) / k1f;

    double KI_e = kminusi_EGF / ki_EGF;
    double KI_S = ((kminusi_S / ki_S)) * 10.0;
    double KI_S1 = kminusi_S1 / ki_S1;
    double KI_S2 = kminusi_S2 / ki_S2;

    // enzyme pools (constants)
    double E_T = 0.003;
    double RAFtotal = 0.0003;
    double phosphoRAF = 0.00003;
    double p1_T = 0.003;
    double p2_T = 0.002;
    double p3_T = 0.002;

    double sTOTAL = 0.02;
    double s1TOTAL = 0.02;
    double s2TOTAL = 0.01;

    double I_EGFR_T = 0.0;
    double I_S_T = 0.5;
    double I_S1_T = 0.5;
    double I_S2_T = 0.0;

    // Populate params array for dYdt_reduced
    // params layout: EGFR_T, RAFtotal, sTOTAL, s1TOTAL, s2TOTAL, I_EGFR_T, I_S_T, I_S1_T, I_S2_T,
    //                phospRAF, p1_T, p2_T, p3_T, KMegf, KI_e, k20, k220,
    //                KMa, KI_S, k1c, ki_S1, k2c, km1c, kminusi_S1, KMb, KMe, KI_S2, KMd, KMf,
    //                k2a, k2b, k2d, k2e, k2f
    params[0] = E_T;           // EGFR_T
    params[1] = RAFtotal;
    params[2] = sTOTAL;
    params[3] = s1TOTAL;
    params[4] = s2TOTAL;
    params[5] = I_EGFR_T;
    params[6] = I_S_T;
    params[7] = I_S1_T;
    params[8] = I_S2_T;
    params[9] = phosphoRAF;    // phospRAF
    params[10] = p1_T;
    params[11] = p2_T;
    params[12] = p3_T;
    params[13] = KMegf;
    params[14] = KI_e;
    params[15] = k20;
    params[16] = k220;
    params[17] = KMa;
    params[18] = KI_S;
    params[19] = k1c;
    params[20] = ki_S1;
    params[21] = k2c;
    params[22] = km1c;
    params[23] = kminusi_S1;
    params[24] = KMb;
    params[25] = KMe;
    params[26] = KI_S2;
    params[27] = KMd;
    params[28] = KMf;
    params[29] = k2a;
    params[30] = k2b;
    params[31] = k2d;
    params[32] = k2e;
    params[33] = k2f;

    long int numero_iterazioni;
    //condizioni iniziali 



    Y[0] = 0.0;
    Y[1] = 0.0;
    Y[2] = 0.0;
    Y[3] = 0.0;





    step = 0.01;
    numero_iterazioni = 12000000;
    int write_every = 10000;  // Write only every 10000 iterations to reduce I/O overhead

    ofstream fdata;
    fdata.open("solution_no_inhibitors.dat");
    fdata << setiosflags(ios::scientific);


    int k = 0;
    //RK4Step
    double t = 0.0;

    // Start timing
    auto start_time = chrono::high_resolution_clock::now();

#if 0
    for (int i = 0; i <= numero_iterazioni; i++) {

        RK4StepWithC(t, Y, step, 4, dYdt_reduced, C, params);
        t += step;

        // Write only every write_every iterations to reduce file I/O overhead
        if (i % write_every == 0) {
            fdata << t << " " << Y[0] << " " << Y[1] << " " << Y[2] << " " << Y[3] << " "
                << C[0] << " " << C[1] << " " << C[2] << " " << C[3] << " " << C[4] << " " << C[5] << " "
                << C[6] << " " << C[7] << " " << C[8] << " " << C[9] << " " << C[10] << " " << C[11] << " " << C[12] << endl;
        }
    }
#else
    bcm3::logger->SetLogToConsole(bcm3::Logger::Info);

    std::unique_ptr<ODESolver> solver = std::make_unique<ODESolverCVODE>();
    //std::unique_ptr<ODESolver> solver = std::make_unique<ODESolverDP5>();

    UserData ud;
    ud.rhs_evals = 0;
    ud.C = C;
    ud.params = params;
    ud.pfdata = &fdata;

    solver->Initialize(4, reinterpret_cast<void*>(&ud));
    solver->SetTolerance(1e-8, 1e-10);
    //solver->SetSolverParameter("max_dt", 0, 100.0);
    solver->SetDerivativeFunction(&derivative);
    //solver->SetIntegrationStepCallback(&integration_step_cb);

    //OdeVectorReal tp(1000);
    //for (int i = 0; i < 1000; i++) {
    //    tp(i) = i * 120;
    //}

    MatrixReal output;

#if 1
    for (int i = 0; i < 1000; i++) {
        solver->SolveStoreIntegrationPoints(Y, 120000, false);
    }
#else
    solver->SolveStoreIntegrationPoints(Y, 120000, false);
#endif

#endif

    // End timing
    auto end_time = chrono::high_resolution_clock::now();
    auto duration = chrono::duration_cast<chrono::milliseconds>(end_time - start_time);

    cout << "Execution time: " << duration.count() / 1000.0 << " seconds" << endl;

    fdata << endl << endl;

    fdata.close();



    return 0;
}




void dYdt_reduced(double t, const double* Y, double* R, double* C, double* params) {

    // remeber to define callback for cached values and first call 
    // EGFR_T,RAFtotal,sTOTAL,s1TOTAL,s2TOTAL,I_EGFR_T,I_S_T,I_S1_T,I_S2_T,phospRAF,p1_T,p2_T,p3_T,KMegf,KI_e,k20,k220, 
    //    KMa,KI_S,k1c,ki_S1,k2c,   km1c      ,kminusi_S1,KMb,KMe,KI_S2,KMd,KMf,k2a,k2b,k2d,k2e,k2f
    double pRAF_tot = Y[0], pS_tot = Y[1], pS1_tot = Y[2], pS2_tot = Y[3];
    //paramters
    double EGFR_T = params[0], RAFtotal = params[1], sTOTAL = params[2], s1TOTAL = params[3], s2TOTAL = params[4], I_EGFR_T = params[5];
    double I_S_T = params[6], I_S1_T = params[7], I_S2_T = params[8], phospRAF = params[9], p1_T = params[10], p2_T = params[11], p3_T = params[12], KMegf = params[13];
    double  KI_e = params[14], k20 = params[15], k220 = params[16], KMa = params[17], KI_S = params[18], k1c = params[19], ki_S1 = params[20];
    double k2c = params[21], km1c = params[22], kminusi_S1 = params[23], KMb = params[24], KMe = params[25];
    double KI_S2 = params[26], KMd = params[27], KMf = params[28], k2a = params[29], k2b = params[30], k2d = params[31], k2e = params[32], k2f = params[33];
    // intermediate numerical  
    //double c_all_startz,Cs_all;
    // physical
    double C1_EGFR, Cinh1_EGFR, C1_startz, Cinh1_startz, C1, Cinh1, C2, Cinh12, C3, C4, C5, Cinh3, C6, Cinh1_all, Cinh12_all, C2RAF;

    double v0, v01, v1, v2, v3, v4, v5, v6;


    //cahed values 
    double C1_EGFR_prev, Cinh1_EGFR_prev, C1_prev, Cinh1_prev, C5_prev, Cinh3_prev;
    // apply conservation laws to get free enzyme concentrations
    double RAF_tot, S_tot, S1_tot, S2_tot;

    RAF_tot = RAFtotal - pRAF_tot;
    S_tot = sTOTAL - pS_tot;
    S1_tot = s1TOTAL - pS1_tot;
    S2_tot = s2TOTAL - pS2_tot;

    //calculate complexes

    if (t == 0) {
        C1_EGFR_prev = -1.0, Cinh1_EGFR_prev = -1.0, C1_prev = -1.0, Cinh1_prev = -1.0, C5_prev = -1.0, Cinh3_prev = -1.0;

    }
    else {
        C1_EGFR_prev = C[0], Cinh1_EGFR_prev = C[12], C1_prev = C[2], Cinh1_prev = C[9], C5_prev = C[7], Cinh3_prev = C[10];
    }


    C1_EGFR = select_root_scalar(KMegf, KI_e, EGFR_T, RAF_tot, I_EGFR_T, C1_EGFR_prev, 1e-8, 1e-10);
    Cinh1_EGFR = select_root_scalar(KI_e, KMegf, EGFR_T, I_EGFR_T, RAF_tot, Cinh1_EGFR_prev, 1e-8, 1e-10);

    C1_startz = select_root_scalar(KMa, KI_S, pRAF_tot, S_tot, I_S_T, C1_prev, 1e-8, 1e-10);
    Cinh1_startz = select_root_scalar(KI_S, KMa, pRAF_tot, I_S_T, S_tot, Cinh1_prev, 1e-8, 1e-10);

    C2RAF = C_root(phospRAF, pRAF_tot - C1_startz - Cinh1_startz, KMb);

    C1 = select_root_scalar(KMa, KI_S, pRAF_tot - C2RAF, S_tot, I_S_T, C1_prev, 1e-8, 1e-10);
    Cinh1 = select_root_scalar(KI_S, KMa, pRAF_tot - C2RAF, I_S_T, S_tot, Cinh1_prev, 1e-8, 1e-10);



    Solution c_all_startz;
    if (t == 0) {
        double init[3] = { 0.0,0.0,0.0 };
        c_all_startz = solve_Cs_allosteric(pS_tot, I_S1_T, S1_tot,
            k1c, ki_S1, ki_S1, k2c, km1c, kminusi_S1, kminusi_S1,
            true, init, false);

    }
    else {
        double init[3] = { C[4], C[5], C[11] };
        c_all_startz = solve_Cs_allosteric(pS_tot, I_S1_T, S1_tot,
            k1c, ki_S1, ki_S1, k2c, km1c, kminusi_S1, kminusi_S1,
            false, init, false);

    }
    //TEsT CROOT 
    C2 = C_root(p1_T, pS_tot - c_all_startz.C1 - c_all_startz.C2 - c_all_startz.C3, KMb);
    Solution Cs_all;
    if (t == 0) {
        double init[3] = { 0.0,0.0,0.0 };
        Cs_all = solve_Cs_allosteric(pS_tot - C2, I_S1_T, S1_tot,
            k1c, ki_S1, ki_S1, k2c, km1c, kminusi_S1, kminusi_S1,
            true, init, false);

    }
    else {
        double init[3] = { C[4], C[5], C[11] };
        Cs_all = solve_Cs_allosteric(pS_tot - C2, I_S1_T, S1_tot,
            k1c, ki_S1, ki_S1, k2c, km1c, kminusi_S1, kminusi_S1,
            false, init, false);

    }



    C3 = Cs_all.C1, Cinh1_all = Cs_all.C2, Cinh12_all = Cs_all.C3;



    // how we have to define the last module ? 
    C4 = C_root(p2_T, pS1_tot - select_root_scalar(KMe, KI_S2, pS1_tot, S2_tot, I_S2_T, C5_prev, 1e-8, 1e-10) - select_root_scalar(KI_S2, KMe, pS1_tot, I_S2_T, S2_tot, Cinh3_prev, 1e-8, 1e-10), KMd);


    //LAYER3
    C5 = select_root_scalar(KMe, KI_S2, pS1_tot - C4, S2_tot, I_S2_T, C5_prev, 1e-8, 1e-10);
    Cinh3 = select_root_scalar(KI_S2, KMe, pS1_tot - C4, I_S2_T, S2_tot, Cinh3_prev, 1e-8, 1e-10);
    C6 = C_root(p3_T, pS2_tot, KMf);

    //C1_EGFR_prev= C[0],Cinh1_EGFR_prev=C[12],C1_prev=C[2],Cinh1_prev=C[9],C5_prev=C[7],Cinh3_prev = C[10]
    C[0] = C1_EGFR;
    C[1] = C2RAF;
    C[2] = C1;
    C[3] = C2;
    C[4] = C3;
    C[5] = Cinh1_all;
    C[6] = C4;
    C[7] = C5;
    C[8] = C6;
    C[9] = Cinh1;
    C[10] = Cinh3;
    C[11] = Cinh12_all;
    C[12] = Cinh1_EGFR;

    v0 = k20 * C1_EGFR;
    v01 = k220 * C2RAF;
    v1 = k2a * C1;
    v2 = k2b * C2;
    v3 = k2c * C3;
    v4 = k2d * C4;
    v5 = k2e * C5;
    v6 = k2f * C6;

    R[0] = v0 - v01;
    R[1] = v1 - v2;
    R[2] = v3 - v4;
    R[3] = v5 - v6;
}








static const int max_neq = 8;

// RK4Step that passes C array
void RK4StepWithC(double t, double* Y, double h, int Neq,
    void(*RHS_Func)(double, double*, double*, double*, double*), double* C, double* params) {
    assert(Neq <= max_neq);
    int i;
    double Y1[max_neq], k1[max_neq], k2[max_neq], k3[max_neq], k4[max_neq];

    RHS_Func(t, Y, k1, C, params);
    for (i = 0; i < Neq; i++) {
        Y1[i] = Y[i] + 0.5 * h * k1[i];
    }

    RHS_Func(t + 0.5 * h, Y1, k2, C, params);
    for (i = 0; i < Neq; i++) {
        Y1[i] = Y[i] + 0.5 * h * k2[i];
    }

    RHS_Func(t + 0.5 * h, Y1, k3, C, params);
    for (i = 0; i < Neq; i++) {
        Y1[i] = Y[i] + h * k3[i];
    }

    RHS_Func(t + h, Y1, k4, C, params);
    for (i = 0; i < Neq; i++) {
        Y[i] += h / 6.0 * (k1[i] + 2 * k2[i] + 2 * k3[i] + k4[i]);
    }
}


////////////////////////////

double C_root(double E_tot, double S_tot, double KM) {
    double disc = (E_tot + S_tot + KM) * (E_tot + S_tot + KM) - 4.0 * E_tot * S_tot;
    disc = max(disc, 0.0);  // guard against tiny negatives
    return 0.5 * ((E_tot + S_tot + KM) - sqrt(disc));
}


//////////

// Halley's method with bisection fallback (mimics Python implementation)
double cubic_root_clamped_scalar(double a3, double a2, double a1, double a0,
    double lo, double hi, double x0,
    double ftol = 1e-12, double xtol = 1e-12, int maxit = 12) {
    double x = min(max(x0, lo), hi);

    // Quick accept
    double f = ((a3 * x + a2) * x + a1) * x + a0;
    if (abs(f) <= ftol) {
        return x;
    }

    bool have_bracket = false;
    double flo = 0.0, fhi = 0.0;

    for (int iter = 0; iter < maxit; ++iter) {
        // Halley's method: x_n+1 = x_n - (2*f*f') / (2*f'^2 - f*f'')
        double fp = (3.0 * a3 * x + 2.0 * a2) * x + a1;
        double fpp = 6.0 * a3 * x + 2.0 * a2;
        double denom = 2.0 * fp * fp - f * fpp;

        double xn = (denom != 0.0) ? x - (2.0 * f * fp) / denom : x;

        // Check if Halley step is out of bounds or invalid
        if (!(lo <= xn && xn <= hi) || !isfinite(xn)) {
            // Try bisection
            if (!have_bracket) {
                flo = ((a3 * lo + a2) * lo + a1) * lo + a0;
                fhi = ((a3 * hi + a2) * hi + a1) * hi + a0;
                have_bracket = (flo * fhi) <= 0.0;
            }

            if (have_bracket) {
                if (flo * f < 0.0) {
                    hi = x;
                    fhi = f;
                }
                else {
                    lo = x;
                    flo = f;
                }
                xn = 0.5 * (lo + hi);
            }
            else {
                // Newton fallback
                xn = (fp != 0.0) ? x - f / fp : 0.5 * (lo + hi);
                xn = min(max(xn, lo), hi);
            }
        }

        // Check convergence
        if (abs(xn - x) <= xtol * max(1.0, abs(x))) {
            return xn;
        }

        x = xn;
        f = ((a3 * x + a2) * x + a1) * x + a0;

        if (abs(f) <= ftol) {
            return x;
        }
    }

    return x;
}

// Solve quadratic in interval [lo, hi]
double solve_quadratic_in_interval(double a2, double a1, double a0, double lo, double hi) {
    if (a2 == 0.0) {
        if (a1 == 0.0) return 0.0;
        double x = -a0 / a1;
        return min(max(x, lo), hi);
    }

    double disc = a1 * a1 - 4.0 * a2 * a0;
    if (disc < 0.0) {
        // No real roots, return endpoint with smaller residual
        double flo = a0;
        double fhi = a2 * hi * hi + a1 * hi + a0;
        return (abs(flo) < abs(fhi)) ? lo : hi;
    }

    double s = sqrt(disc);
    double q = -0.5 * (a1 + copysign(s, a1));
    double r1 = q / a2;
    double r2 = (q != 0.0) ? a0 / q : (-a1 - copysign(s, a1)) / (2.0 * a2);

    // Check if roots are in interval
    if (lo <= r1 && r1 <= hi) return r1;
    if (lo <= r2 && r2 <= hi) return r2;

    // Clamp and return closest
    double r1c = min(max(r1, lo), hi);
    double r2c = min(max(r2, lo), hi);
    return (abs(r1c - r1) <= abs(r2c - r2)) ? r1c : r2c;
}

void phy_coeffs(double k1, double k2, double E, double S1, double S2, double coeffs[4]) {
    coeffs[0] = -(k1 - k2);
    coeffs[1] = (E + k1 + S1) * (k1 - k2) - (S1 * k2 + S2 * k1);
    coeffs[2] = (-E * (k1 - k2) + (S1 * k2 + S2 * k1) + k2 * (E + k1)) * S1;
    coeffs[3] = -E * k2 * S1 * S1;
}

// Main root selector with warm-start and adaptive tolerances (mimics Python)
double select_root_scalar(double k1, double k2, double E, double S1, double S2,
    double x_prev = -1.0, double rtol = 1e-8, double atol = 1e-10) {
    double bnd = min(E, S1);
    if (bnd <= 0.0) {
        return 0.0;
    }

    // Get polynomial coefficients
    double dk = k1 - k2;
    double a3 = -dk;
    double a2 = (E + k1 + S1) * dk - (S1 * k2 + S2 * k1);
    double a1 = (-E * dk + (S1 * k2 + S2 * k1) + k2 * (E + k1)) * S1;
    double a0 = -E * k2 * S1 * S1;

    // Check for degeneracy
    if (abs(a3) < 1e-14) {
        return solve_quadratic_in_interval(a2, a1, a0, 0.0, bnd);
    }

    // Scale to y ∈ [0,1] for numerical stability
    double invb = 1.0 / bnd;
    double A3 = a3;
    double A2 = a2 * invb;
    double invb2 = invb * invb;
    double A1 = a1 * invb2;
    double A0 = a0 * invb2 * invb;

    // Adaptive tolerances
    double xtol_x = max(1e-12, rtol * bnd + atol);
    double xtol_y = xtol_x * invb;
    double ftol_y = 1e-12;

    // Warm start
    double x0 = 0.5 * bnd;
    if (x_prev >= 0.0 && isfinite(x_prev)) {
        x0 = x_prev;
    }
    double y0 = min(max(x0 * invb, 0.0), 1.0);

    // Solve in scaled space and convert back
    double y = cubic_root_clamped_scalar(A3, A2, A1, A0, 0.0, 1.0, y0, ftol_y, xtol_y, 12);
    return y * bnd;
}


//////////////////////////////////////////////////////////////



bool solveLinear3x3(double A[3][3], double b[3], double x[3]) {
    constexpr double EPS = 1e-12;

    // Forward elimination
    for (int i = 0; i < 3; ++i) {
        int maxRow = i;
        for (int k = i + 1; k < 3; ++k)
            if (fabs(A[k][i]) > fabs(A[maxRow][i]))
                maxRow = k;

        // Swap rows
        for (int j = 0; j < 3; ++j)
            swap(A[i][j], A[maxRow][j]);
        swap(b[i], b[maxRow]);

        if (fabs(A[i][i]) < EPS) return false;

        // Eliminate
        for (int k = i + 1; k < 3; ++k) {
            double f = A[k][i] / A[i][i];
            for (int j = i; j < 3; ++j)
                A[k][j] -= f * A[i][j];
            b[k] -= f * b[i];
        }
    }

    // Back substitution
    for (int i = 2; i >= 0; --i) {
        x[i] = b[i];
        for (int j = i + 1; j < 3; ++j)
            x[i] -= A[i][j] * x[j];
        x[i] /= A[i][i];
    }

    return true;
}

void computeEquations(double C1, double C2, double C3,
    double E_T, double I_hat, double S_hat,
    double a1, double a2, double i, double k1, double d1, double d2, double i_rev,
    double F[3]) {
    F[0] = E_T * S_hat * a1
        - E_T * C1 * a1
        - E_T * C3 * a1
        - I_hat * C1 * i
        - S_hat * C1 * a1
        - S_hat * C2 * a1
        - S_hat * C3 * a1
        + C1 * C1 * a1
        + C1 * C2 * a1
        + C1 * C2 * i
        + 2 * C1 * C3 * a1
        + C1 * C3 * i
        - C1 * d1
        - C1 * k1
        + C2 * C3 * a1
        + C3 * C3 * a1
        + C3 * i_rev;

    F[1] = E_T * I_hat * a2
        - E_T * C2 * a2
        - E_T * C3 * a2
        - I_hat * C1 * a2
        - I_hat * C2 * a2
        - I_hat * C3 * a2
        + C1 * C2 * a2
        + C1 * C3 * a2
        + C2 * C2 * a2
        + 2 * C2 * C3 * a2
        - C2 * d2
        + a2 * C3 * C3;

    F[2] = I_hat * C1 * i
        - C1 * C2 * i
        - C1 * C3 * i
        - C3 * i_rev;
}

void computeJacobian(double C1, double C2, double C3,
    double E_T, double I_hat, double S_hat,
    double a1, double a2, double i, double k1, double d1, double d2, double i_rev,
    double J[3][3]) {
    J[0][0] = -E_T * a1 - I_hat * i - S_hat * a1 + 2.0 * C1 * a1 + a1 * C2 + C2 * i + 2.0 * C3 * a1 + C3 * i - d1 - k1;
    J[0][1] = -S_hat * a1 + C1 * a1 + C1 * i + C3 * a1;
    J[0][2] = -E_T * a1 - S_hat * a1 + 2.0 * C1 * a1 + C1 * i + C2 * a1 + 2.0 * C3 * a1 + i_rev;

    J[1][0] = -I_hat * a2 + C2 * a2 + C3 * a2;
    J[1][1] = -E_T * a2 - I_hat * a2 + C1 * a2 + 2.0 * C2 * a2 + 2.0 * C3 * a2 - d2;
    J[1][2] = -E_T * a2 - I_hat * a2 + C1 * a2 + 2.0 * C2 * a2 + 2.0 * C3 * a2;

    J[2][0] = I_hat * i - C2 * i - C3 * i;
    J[2][1] = -C1 * i;
    J[2][2] = -C1 * i - i_rev;
}

Solution solve_Cs_allosteric(double E_T, double I_hat, double S_hat, double a1, double a2, double i, double k1, double d1, double d2, double i_rev, bool init_p, double init_v[3], bool log_i) {
    constexpr int max_iter = 100;
    constexpr double tol = 1e-12;

    double best_sol[3] = { 0.0, 0.0, 0.0 };
    bool found = false;
    // Needs proper bounds 
    // Log inputs
    if (log_i) {
        cout << "solve_Cs_allosteric called: E_T=" << E_T
            << ", I_hat=" << I_hat << ", S_hat=" << S_hat
            << ", a1=" << a1 << ", a2=" << a2 << ", i=" << i
            << ", k1=" << k1 << ", d1=" << d1 << ", d2=" << d2
            << ", i_rev=" << i_rev << endl;
    }

    for (double x = init_p ? 0.0 : init_v[0]; x + init_v[0] <= 0.5; x += init_p ? 0.01 : 1) {// std::min(E_T,S_hat))
        for (double y = init_p ? 0.0 : init_v[1]; y + init_v[1] <= 0.5; y += init_p ? 0.01 : 1) {//std::min(E_T,I_hat))
            for (double z = init_p ? 0.0 : init_v[2]; z + init_v[2] <= 0.5; z += init_p ? 0.01 : 1) {//std::min({E_T,I_hat,S_hat}))

                double C[3] = { x, y, z };
                // solver steps 
                for (int iter = 0; iter < max_iter; ++iter) {
                    double F[3], J[3][3];
                    computeEquations(C[0], C[1], C[2], E_T, I_hat, S_hat, a1, a2, i, k1, d1, d2, i_rev, F);
                    computeJacobian(C[0], C[1], C[2], E_T, I_hat, S_hat, a1, a2, i, k1, d1, d2, i_rev, J);

                    double dx[3];
                    double J_copy[3][3], F_copy[3];
                    for (int ii = 0; ii < 3; ++ii) {
                        F_copy[ii] = F[ii];
                        for (int jj = 0; jj < 3; ++jj)
                            J_copy[ii][jj] = J[ii][jj];
                    }

                    if (!solveLinear3x3(J_copy, F_copy, dx)) break;

                    for (int i = 0; i < 3; ++i) {
                        C[i] -= dx[i];
                        C[i] = std::max(0.0, C[i]);
                    }
                    //if (log_i) {
                    //            cout << "  residual : [" << fabs(dx[0]) + fabs(dx[1]) + fabs(dx[2])<<" " <<C[0]<<" "<<C[1]<<" "<<C[2]<<" " << endl;
                    //        }
                    if (fabs(dx[0]) + fabs(dx[1]) + fabs(dx[2]) < tol) {
                        if (C[0] >= 0 && C[1] >= 0 && C[2] >= 0) {
                            best_sol[0] = C[0];
                            best_sol[1] = C[1];
                            best_sol[2] = C[2];
                            found = true;
                            if (log_i) {
                                cout << "  ✓ Solution found: [" << C[0] << ", " << C[1] << ", " << C[2] << "]" << endl;
                            }
                        }
                        break;  // replaces goto: once found, stop iterating for this initial guess  you break the solver because tollerance reached 
                    }
                }

                if (found) break;
            }
            if (found) break;
        }
        if (found) break;
        // early exit outer loops once found
    }
    if (log_i && !found) {
        cout << "  ✗ No solution found!" << endl;
    }

    return { best_sol[0], best_sol[1], best_sol[2], found };
}

#if 0

int main(int argc, char* argv[])
{
	MatrixReal test(3, 3);
	test.row(0) << 2.2,  5.4, 1.1;
	test.row(1) << -2.5, 3.3, 8.0;
	test.row(2) << 7.0,  0.5, 0.1;

	MatrixReal matching(3, 3);
	bool result = HungarianMatching(test, matching, HUNGARIAN_MATCH_MAX);

	std::cout << matching;

	std::vector<WeightedBipartiteEdge> edges;
	for (int i = 0; i < 3; i++) {
		for (int j = 0; j < 3; j++) {
			edges.push_back(WeightedBipartiteEdge(i, j, -test(i, j)));
		}
	}

	std::vector<int> matching2 = hungarianMinimumWeightPerfectMatching(3, edges, 9);

	return 0;
}

#endif