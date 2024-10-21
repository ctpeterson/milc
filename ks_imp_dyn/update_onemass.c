/******************************************************************** 

Rewrite of update_onemass.c because it needed it. Includes 2nd-order
Omelyan integrator. Does not include R-algorithm because nobody
uses it anymore.

Note: "steps" in the MILC convetion means "number of force evaluations",
not the number of full Omelyan steps

(Rewrite) author: Curtis Taylor Peterson
Original author(s): Unknown :/

********************************************************************/

#include "ks_imp_includes.h"

int update(){
  int step,iters = 0;
  double d_action();
  double startaction,endaction,xrandom;
  Real final_rsq;
  Real lmbda = 0.1931833275037836;
  imp_ferm_links_t** fn;

  /* ---- Utility functions for HMC ---- */

  // Get Hasenbusch phi field
  void grsource_imp_hasenbusch(){
    // det(M)/det(H)
    grsource_imp(F_OFFSET(phi), mass, EVENANDODD, fn[0]);
    iters += ks_congrad( 
      F_OFFSET(phi), F_OFFSET(xxx), hmass, 
		  aniter, anrestart, arsqmin, MILC_PRECISION, 
		  EVENANDODD, &final_rsq, fn[0] 
    );
    dslash_site(F_OFFSET(xxx), F_OFFSET(hphi), EVEN, fn[0]);
    scalar_mult_add_latvec(
      F_OFFSET(hphi), 
      F_OFFSET(xxx),
      2.0*hmass,
      F_OFFSET(hphi), 
      EVEN
    );
    clear_latvec(F_OFFSET(phi), EVENANDODD);
    clear_latvec(F_OFFSET(xxx), EVENANDODD);
    // det(H)
    grsource_imp(F_OFFSET(phi), hmass, EVEN, fn[0]);
  }

  // Momentum/fermion heatbath
  void heatbath(){
    ranmom(); 
    #ifdef HASENBUSCH
      grsource_imp_hasenbusch();
    #else
      grsource_imp(F_OFFSET(phi), mass, EVEN, fn[0]); 
    #endif
  }

  // Fermion solve D = -1.714319e-03
  void solve(int niter, int nrestart, Real rsqmin){
    #ifdef HASENBUSCH
      // det(M)/det(H)
      iters += ks_congrad( 
        F_OFFSET(hphi), F_OFFSET(hxxx), mass, 
        niter, nrestart, rsqmin, MILC_PRECISION, 
        EVEN, &final_rsq, fn[0] 
      );
      // det(H)
      iters += ks_congrad( 
        F_OFFSET(phi), F_OFFSET(xxx), hmass, 
        niter, nrestart, rsqmin, MILC_PRECISION, 
        EVEN, &final_rsq, fn[0] 
      );
    #else
      iters += ks_congrad( 
        F_OFFSET(phi), F_OFFSET(xxx), mass, 
        niter, nrestart, rsqmin, MILC_PRECISION, 
        EVEN, &final_rsq, fn[0] 
      );
    #endif
  }

  // Smear fields
  void smear(){
    restore_fermion_links_from_site(fn_links, MILC_PRECISION);
	  fn = get_fm_links(fn_links);
  }

  // Momentum update wrapper
  void update_v(int step, double dtau){
    // CG solve
    smear(); solve(fniter,fnrestart,frsqmin);
    #ifdef HASENBUSCH
      dslash_site(F_OFFSET(hxxx), F_OFFSET(hxxx), ODD, fn[0]);
    #endif
    dslash_site(F_OFFSET(xxx), F_OFFSET(xxx), ODD, fn[0]);

    // Integrator step
    #ifdef INT_LEAPFROG
      update_h(dtau);
    #elif defined INT_OMELYAN
      update_h(dtau);
    #elif defined INT_OMELYAN_3G1F
      update_h_fermion(dtau);
    #endif
  }

  // Nested gauge update
  void update_u_gauge(int step, int substep, Real dtau){
    if (step == 1){update_u(lmbda*dtau);}
    else {
      if (substep == 1){update_u(2.0*lmbda*dtau);}
      else {update_u(lmbda*dtau);}
    }
    update_h_gauge(0.5*dtau);
    update_u((1.0-2.0*lmbda)*dtau);
    update_h_gauge(0.5*dtau);
    if (step == steps){update_u(lmbda*dtau);}
    else {if (substep != 3){update_u(lmbda*dtau);}}
  }

  /* ---- Run HMC trajectory ---- */

  // Do heatbath & calculate initial Hamiltonian
  smear(); heatbath(); 
  solve(aniter,anrestart,arsqmin); 
  startaction = d_action();

  // Save backup gauge field
	gauge_field_copy(F_OFFSET(link[0]), F_OFFSET(old_link[0]));

  // Molecular dynamics trajectory
  for(step=1; step <= steps; step++){
    #ifdef INT_LEAPFROG
      if (step == 1){update_u(0.5*epsilon);}
      else {update_u(epsilon);}
      update_v(step,epsilon);
      if (step == steps){update_u(0.5*epsilon);}
    #elif defined INT_OMELYAN
      if (step == 1){update_u(lmbda*epsilon);}
      else {update_u(2.0*lmbda*epsilon);}
      update_v(step,0.5*epsilon);
      update_u((1.0-2.0*lmbda)*epsilon);
      update_v(step,0.5*epsilon);
      if (step == steps){update_u(lmbda*epsilon);}
    #elif defined INT_OMELYAN_3G1F
      update_u_gauge(step,1,lmbda*epsilon); 
      update_v(step,0.5*epsilon);
      update_u_gauge(step,2,(1.0-2.0*lmbda)*epsilon); 
      update_v(step,0.5*epsilon);
      update_u_gauge(step,3,lmbda*epsilon); 
    #endif
  }

  // Calculate final Hamiltonian
  smear(); solve(aniter,anrestart,arsqmin); 
  endaction = d_action();

  // Metropolis
  if(this_node==0)xrandom = myrand(&node_prn);
  broadcast_float(&xrandom);
  if( exp( (double)(startaction-endaction) ) < xrandom ){
	  if(steps > 0)
	    gauge_field_copy(F_OFFSET(old_link[0]), F_OFFSET(link[0]));
        #ifdef FN
	        invalidate_fermion_links(fn_links);
        #endif
	    node0_printf("REJECT: delta S = %e\n", (double)(endaction-startaction));
  }
  else {node0_printf("ACCEPT: delta S = %e\n", (double)(endaction-startaction));}

  if(steps > 0)return (iters/steps);
  else return(-99);
}