#include <stdio.h>

#define NMAX 1000

// Fortran subroutine declaration
extern void dsimthnf_(int *n, double *wth, double *ha, double *s0, double *pr,
                      double *psi1, double *psi2, double *psi3, double *psi4,
                      double *p0, double *hg, double *ec, double *beta,
                      double *alphae, double *rd, double *lambda, double *mi,
                      double *rho1, double *rho2, double *rho3, double *rho4, double *rhof,
                      double *k1, double *k2, double *k3, double *k4, double *kf,
                      double *cp1, double *cp2, double *cp3, double *cp4, double *cpf,
                      double *sig1, double *sig2, double *sig3, double *sig4, double *sigf,
                      double *tol, int *maxit, double *W,
                      double *f, double *p, double *ht, double *g, double *s,
                      double *q1, double *q2);

int main() {
    int N = 101;
    double WTH = 1.0;

    double HA = 1.0;
    double S0 = 0.5;
    double PR = 6.2;

    // Volume fractions (phy1 to phy4)
    double PSI1 = 0.01; // Silver (Ag)
    double PSI2 = 0.01; // Gold (Au)
    double PSI3 = 0.01; // Tantalum (Ta)
    double PSI4 = 0.01; // Copper (Cu)

    // Other physical parameters
    double P0 = 1.0;
    double HG = 0.5;
    double EC = 0.1;
    double BETA = 0.5;
    double ALPHAE = 0.1;
    double RD = 0.1;
    double LAMBDA = 0.1;
    double MI = 0.1;

    // Densities (ρ)
    double RHO1 = 10500.0;   // Ag
    double RHO2 = 19320.0;   // Au
    double RHO3 = 16650.0;   // Ta
    double RHO4 = 8933.0;    // Cu
    double RHOF = 1063.0;    // Blood

    // Thermal conductivities (k)
    double K1 = 429.0;
    double K2 = 314.0;
    double K3 = 0.52;
    double K4 = 401.0;
    double KF = 0.492;

    // Specific heat capacities (Cp)
    double CP1 = 235.0;
    double CP2 = 129.0;
    double CP3 = 686.2;
    double CP4 = 385.0;
    double CPF = 3594.0;

    // Electrical conductivities (σ)
    double SIG1 = 6.30e7;
    double SIG2 = 4.10e7;
    double SIG3 = 7.70e6;
    double SIG4 = 5.96e7;
    double SIGF = 6.67e-1;

    // Solver controls
    double TOL = 1e-7;
    int MAXIT = 100000;
    double W = 0.95;

    // Output arrays
    double F[NMAX], P[NMAX], HT[NMAX], G[NMAX], S[NMAX];

    // Output values for skin friction and Nusselt number
    double Q1, Q2;

    // Call Fortran subroutine
    dsimthnf_(&N, &WTH, &HA, &S0, &PR,
              &PSI1, &PSI2, &PSI3, &PSI4,
              &P0, &HG, &EC, &BETA, &ALPHAE, &RD, &LAMBDA, &MI,
              &RHO1, &RHO2, &RHO3, &RHO4, &RHOF,
              &K1, &K2, &K3, &K4, &KF,
              &CP1, &CP2, &CP3, &CP4, &CPF,
              &SIG1, &SIG2, &SIG3, &SIG4, &SIGF,
              &TOL, &MAXIT, &W,
              F, P, HT, G, S,
              &Q1, &Q2);

    // Display results
    printf("Skin Friction Coefficient (Q1): %lf\n", Q1);
    printf("Nusselt Number (Q2): %lf\n", Q2);

    return 0;
}
