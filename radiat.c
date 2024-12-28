/* ///////////////////////////////////////////////////////////////////// */
/*! 
  \file  
  \brief Compute right hand side for Tabulated cooling

  \authors A. Mignone (mignone@ph.unito.it)\n
           M. Sormani\n

 \b References

  \date   Apr 09, 2019
*/
/* ///////////////////////////////////////////////////////////////////// */
#include "pluto.h"
#include "local_pluto.h"

extern double gCooling_x1, gCooling_x2, gCooling_x3;

/* ***************************************************************** */
void Radiat (double *v, double *rhs)
/*!
 *   Provide r.h.s. for tabulated cooling.
 * 
 ******************************************************************* */
{
  int    klo, khi, kmid;
  static int ntab;
  double  mu, muH, mue, T, Tmid, scrh, dT, prs;
  static double *L_tab, *T_tab, E_cost;
  double dummy[4];
  
/* -------------------------------------------
        Read tabulated cooling function
   ------------------------------------------- */

  if (T_tab == NULL){
    FILE *fcool;
    printLog (" > Reading cooltable.dat from disk...\n");
    fcool = fopen("cooltable.dat","r");
    if (fcool == NULL){
      printLog ("! Radiat: cooltable.dat could not be found.\n");
      QUIT_PLUTO(1);
    }
    L_tab = ARRAY_1D(20000, double);
    T_tab = ARRAY_1D(20000, double);

    ntab = 0;
    while (fscanf(fcool, "%lf  %lf\n", T_tab + ntab, 
                                       L_tab + ntab)!=EOF) {
      ntab++;
    }
    E_cost = UNIT_LENGTH/UNIT_DENSITY/pow(UNIT_VELOCITY, 3.0);
    fclose(fcool);
  }

/* ---------------------------------------------
            Get pressure and temperature 
   --------------------------------------------- */

  prs = v[RHOE]*(g_gamma-1.0);
  if (prs < 0.0) {
    prs     = g_smallPressure;
    v[RHOE] = prs/(g_gamma - 1.0);
  }

  mu  = MeanMolecularWeight(v, dummy);
  muH = dummy[2];
  mue = dummy[0];
  T   = prs/v[RHO]*KELVIN*mu;

  if (T != T || T==0.){
    printLog ("! Radiat(): Nan found in f_m_cl 1: rho = %12.6e, prs = %12.6e\n",v[RHO], prs);
    QUIT_PLUTO(1);
  }

  double T_cl = 1.e4;
  double f_m_cl = v[TRC];
 // printLog ("TRC = %12.6e\n", v[TRC]);
  const double TOLERANCE = 1.e-6;

  if (fabs(f_m_cl - 1.0) <= TOLERANCE) {
  //if (f_m_cl == 1.0){
  rhs[RHOE] = 0.0;
  }
  else {
  double rho_cl = prs*muH*CONST_mp*UNIT_VELOCITY*UNIT_VELOCITY/(T_cl * CONST_kB);  // cold gas density
  double f_v_h = 1 - v[RHO]*f_m_cl/rho_cl; // hot gas volume fraction
  double f_v_cl = 1 - f_v_h;  // cold gas volume fraction
  double rho_h = (v[RHO] - f_v_cl*rho_cl)/f_v_h;
  T = rho_cl*T_cl/rho_h;
  if (T != T || T==0.){

    printLog ("! Radiat(): Nan found: rho = %12.6e, rho_cl = %12.6e, f_v_h = %12.6e, f_m_cl = %12.6e, rho_h = %12.6e\n",v[RHO], rho_cl, f_v_h, f_m_cl, rho_h);
    QUIT_PLUTO(1);
  }

/* ----------------------------------------------
        Table lookup by binary search  
   ---------------------------------------------- */

  klo = 0;
  khi = ntab - 1;

  if (T > T_tab[khi] || T < T_tab[klo]){
    rhs[RHOE] = 0.0;
    return;
    QUIT_PLUTO(1);
  }

  while (klo != (khi - 1)){
    kmid = (klo + khi)/2;
    Tmid = T_tab[kmid];
    if (T <= Tmid){
      khi = kmid;
    }else if (T > Tmid){
      klo = kmid;
    }
  }

/* -----------------------------------------------
    Compute r.h.s
   ----------------------------------------------- */
  double nH = rho_h*UNIT_DENSITY/(muH*CONST_amu); //UNIT_DENSITY/CONST_amu*H_MASS_FRAC/CONST_AH*v[RHO];
  double ne = rho_h*UNIT_DENSITY/(mue*CONST_amu); //nH*(1.0 + 0.5*CONST_AZ*FRAC_Z);
  double n  = rho_h*UNIT_DENSITY/(mu*CONST_amu);
  dT        = T_tab[khi] - T_tab[klo];
  scrh      = L_tab[klo]*(T_tab[khi] - T)/dT + L_tab[khi]*(T - T_tab[klo])/dT;
  rhs[RHOE] = -nH*nH*scrh*E_cost*f_v_h;
  }
  
/* ----------------------------------------------
    Temperature cutoff
   ---------------------------------------------- */

  //rhs[RHOE] *= 1.0 - 1.0/cosh( pow( T/g_minCoolingTemp, 12));
}
