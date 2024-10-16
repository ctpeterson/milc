/******************************************************************** 

Rewrite of update_onemass.c because it needed it. Includes 2nd-order
Omelyan integrator. I tried to keep the code as close to the original
structure as I could, even though I don't like how it was 
originally written. Does not include R-algorithm because nobody
uses it anymore.

Note: "steps" in the MILC convetion means "number of force evaluations",
not the number of full Omelyan steps

(Rewrite) author: Curtis Taylor Peterson
Original author(s): Unknown :/

********************************************************************/

#include "ks_imp_includes.h"


#ifdef HASENBUSCH
  // "phi" field for Hasenbusch; based on legacy code & 
  // David Schaich's KS_nHYP_FA code
  /*
  int gsource_imp_hasenbusch(
      field_offset dest, 
      Real mass1, 
      Real mass2,
      Real rsqmin,
      int niter, 
      int nrestart,
      int parity,  
      imp_ferm_links_t *fn
    ){
    int iters = 0; 
    Real final_rsq;
    quark_invert_control qic;
    su3_vector *g_rand = create_v_field();
    su3_vector *d1g_rand = create_v_field();
    su3_vector *d2id1g_rand = create_v_field();

    // Set up quark invert control for field solve
    qic.prec = MILC_PRECISION;
    qic.min = 0;
    qic.max = iters;
    qic.nrestart = nrestart;
    qic.parity = EVENANDODD;
    qic.start_flag = 0;
    qic.nsrc = 1;
    qic.resid = sqrt(rsqmin);
    qic.relresid = 0;

    // D1*g
    grsource_plain_field(g_rand, EVENANDODD);
    ks_dirac_adj_op(g_rand, d1g_rand, mass1, EVENANDODD, fn);

    // D2^{-1}*D1*g
    iters = ks_congrad_field(d1g_rand, d2id1g_rand, &qic, mass2, fn);

    // Clear temporary fields from memory
    destroy_v_field(g_rand);
    destroy_v_field(d1g_rand);
    destroy_v_field(d2id1g_rand);

    // Create source
    FORALLSITES(i,s){
      for(j=0;j<3;j++){
        //s->xxx.c[j] = complex_gaussian_rand_no(&node_prn);
        s->xxx.c[j] = complex_gaussian_rand_no(&(s->site_prn));
      }
    }
    grsource_plain_field(g_rand, EVENANDODD);
    copy_site_member_from_v_field(F_OFFSET(xxx), g_rand);

    // Hit source w/ D(m)^{dag}
    dslash_site(F_OFFSET(xxx), F_OFFSET(hxxx), EVENANDODD, fn);
    scalar_mult_latvec(F_OFFSET(hxxx), -1.0, F_OFFSET(hxxx), EVENANDODD);
    scalar_mult_add_latvec( 
      F_OFFSET(hxxx), F_OFFSET(xxx), 
      2.0*mass1, F_OFFSET(hxxx), EVENANDODD 
    );

    // Hit D(m)^{dag}*"source" w/ D(m_h)^{-1}
    iters = ks_congrad( 
      F_OFFSET(hxxx), F_OFFSET(hpsi), mass2, 
      niter, nrestart, rsqmin, MILC_PRECISION, 
      EVENANDODD, &final_rsq, fn 
    );
    dslash_site(F_OFFSET(hpsi), dest, parity, fn);
    scalar_mult_add_latvec(dest, F_OFFSET(hpsi), dest, 2.0*mass2, parity);

    // Clear lattice vectors for safety
    //clear_latvec(F_OFFSET(xxx), EVENANDODD); 
    //clear_latvec(F_OFFSET(hxxx), EVENANDODD);
    //clear_latvec(F_OFFSET(hpsi), EVENANDODD);
    destroy_v_field(g_rand);
    

    // Return CG iterations
    return iters;
  }
  */
#endif

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
		  fniter, fnrestart, frsqmin, MILC_PRECISION, 
		  EVENANDODD, &final_rsq, fn[0] 
    );
    dslash_site(F_OFFSET(xxx), F_OFFSET(hphi), EVEN, fn[0]);
    scalar_mult_add_latvec(
      F_OFFSET(hphi), F_OFFSET(xxx),
      2.0*hmass,
      F_OFFSET(hphi), 
      EVEN
    );
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

  // Fermion solve
  void solve(){
    #ifdef HASENBUSCH
      // det(M)/det(H)
      iters += ks_congrad( 
        F_OFFSET(hphi), F_OFFSET(hxxx), mass, 
        aniter, anrestart, arsqmin, MILC_PRECISION, 
        EVEN, &final_rsq, fn[0] 
      );
      // det(H)
      iters += ks_congrad( 
        F_OFFSET(phi), F_OFFSET(xxx), hmass, 
        aniter, anrestart, arsqmin, MILC_PRECISION, 
        EVEN, &final_rsq, fn[0] 
      );
    #else
      iters += ks_congrad( 
        F_OFFSET(phi), F_OFFSET(xxx), mass, 
        aniter, anrestart, arsqmin, MILC_PRECISION, 
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
  void update_v(int step, double dtau, double offset){
    // CG solve
    smear();
    #ifdef HASENBUSCH
      // det(M)/det(H)
      iters += ks_congrad( 
        F_OFFSET(hphi), F_OFFSET(hxxx), mass, 
        fniter, fnrestart, frsqmin, MILC_PRECISION, 
        EVEN, &final_rsq, fn[0] 
      );
      dslash_site(F_OFFSET(hxxx), F_OFFSET(hxxx), ODD, fn[0]);
      // det(H)
      iters += ks_congrad( 
        F_OFFSET(phi), F_OFFSET(xxx), hmass, 
        fniter, fnrestart, frsqmin, MILC_PRECISION, 
        EVEN, &final_rsq, fn[0] 
      );
    #else
      iters += ks_congrad( 
        F_OFFSET(phi), F_OFFSET(xxx), mass, 
        fniter, fnrestart, frsqmin, MILC_PRECISION, 
        EVEN, &final_rsq, fn[0] 
      );
    dslash_site(F_OFFSET(xxx), F_OFFSET(xxx), ODD, fn[0]);
    #endif

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
  smear(); heatbath(); solve(); startaction = d_action();

  // Save backup gauge field
	gauge_field_copy(F_OFFSET(link[0]), F_OFFSET(old_link[0]));

  // Molecular dynamics trajectory
  for(step=1; step <= steps; step++){
    #ifdef INT_LEAPFROG
      if (step == 1){update_u(0.5*epsilon);}
      else {update_u(epsilon);}
      update_v(step,epsilon,0.5);
      if (step == steps){update_u(0.5*epsilon);}
    #elif defined INT_OMELYAN
      if (step == 1){update_u(lmbda*epsilon);}
      else {update_u(2.0*lmbda*epsilon);}
      update_v(step,0.5*epsilon,0.333);
      update_u((1.0-2.0*lmbda)*epsilon);
      update_v(step,0.5*epsilon,0.667);
      if (step == steps){update_u(lmbda*epsilon);}
    #elif defined INT_OMELYAN_3G1F
      update_u_gauge(step,1,lmbda*epsilon); 
      update_v(step,0.5*epsilon,0.333);
      update_u_gauge(step,2,(1.0-2.0*lmbda)*epsilon); 
      update_v(step,0.5*epsilon,0.667);
      update_u_gauge(step,3,lmbda*epsilon); 
    #endif
  }

  // Calculate final Hamiltonian
  smear(); solve(); endaction = d_action();

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