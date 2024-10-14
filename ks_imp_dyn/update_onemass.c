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

int update(){
    void predict_next_xxx(Real *oldtime,Real *newtime,Real *nexttime);
    int step,iters = 0;
    double d_action();
    double startaction,endaction,xrandom;
    Real lmbda = 0.1931833275037836;
    Real final_rsq;
    Real cg_time = 0.0;
    Real old_cg_time = -1.0e6;
    Real next_cg_time = -1.0e6;
    imp_ferm_links_t** fn;

    // Momentum update wrapper
    void update_v(int step, double dtau, double offset){
      next_cg_time = ((Real)step-offset)*epsilon;
	  predict_next_xxx(&old_cg_time,&cg_time,&next_cg_time);
      restore_fermion_links_from_site(fn_links, MILC_PRECISION);
	  fn = get_fm_links(fn_links);
      iters += ks_congrad( 
        F_OFFSET(phi), F_OFFSET(xxx), mass, 
		niter, nrestart, rsqmin, MILC_PRECISION, 
		EVEN, &final_rsq, fn[0] 
      );
	  dslash_site( F_OFFSET(xxx), F_OFFSET(xxx), ODD, fn[0]);
	  cg_time = ((Real)step-offset)*epsilon;
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

    // Define values for first & last gauge update step (minor optimization)
    #ifdef INT_LEAPFROG
      node0_printf("using leapfrog integrator\n");
    #elif defined INT_OMELYAN
      node0_printf("using 2nd-order Omelyan integrator\n");
    #elif defined INT_OMELYAN_3G1F
      node0_printf("using \"3G1F\" nested 2nd-order Omelyan integrator\n");
    #else
      node0_printf("No integration algorithm, or unknown one\n");
      terminate(1);
    #endif

    // First link smearing
    restore_fermion_links_from_site(fn_links, MILC_PRECISION);
	fn = get_fm_links(fn_links);

    // Momentum/fermion heatbath
    ranmom(); grsource_imp(F_OFFSET(phi), mass, EVEN, fn[0]); 
    
    // Calculate initial Hamiltonian
	iters += ks_congrad( 
        F_OFFSET(phi), F_OFFSET(xxx), mass, 
		niter, nrestart, rsqmin, MILC_PRECISION, 
		EVEN, &final_rsq, fn[0] 
    );
    startaction = d_action();

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
      }
    #endif

    // Smear one last time & calculate final Hamiltonian
    next_cg_time = steps*epsilon;
    predict_next_xxx(&old_cg_time,&cg_time,&next_cg_time);
    restore_fermion_links_from_site(fn_links, MILC_PRECISION);
    fn = get_fm_links(fn_links);
    iters += ks_congrad( 
        F_OFFSET(phi), F_OFFSET(xxx), mass,
		niter, nrestart, rsqmin, MILC_PRECISION, 
		EVEN, &final_rsq, fn[0] 
    );
    cg_time = steps*epsilon;
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

void predict_next_xxx(Real *oldtime,Real *newtime,Real *nexttime)
{
register int i;
register site *s;
register Real x;
su3_vector tvec;
    if( *newtime != *oldtime ) x = (*nexttime-*newtime)/(*newtime-*oldtime);
    else x = 0.0;
    if( *oldtime < 0.0 ){
        FOREVENSITES(i,s){
	    s->old_xxx = s->xxx;
        }
    }
    else  {
        FOREVENSITES(i,s){
            sub_su3_vector( &(s->xxx), &(s->old_xxx), &tvec);
	    s->old_xxx = s->xxx;
            scalar_mult_add_su3_vector( &(s->xxx), &tvec,x, &(s->xxx) );
        }
    }
    *oldtime = *newtime;
    *newtime = *nexttime;
}