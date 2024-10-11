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
    double d_action();
    int step,iters = 0;
    double startaction,endaction,xrandom;
    Real lmbda = 0.8;//0.1931833275037836;
    Real final_rsq;
    Real cg_time = 0.0;
    Real old_cg_time = -1.0e6;
    Real next_cg_time = -1.0e6;
    imp_ferm_links_t** fn;

    // Momnetum update wrapper
    void update_v(int step, double offset){
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
        update_h(epsilon);
    }

    // First link smearing
    restore_fermion_links_from_site(fn_links, MILC_PRECISION);
	fn = get_fm_links(fn_links);

    // Momentum/fermion heatbath
    ranmom();
    grsource_imp(F_OFFSET(phi), mass, EVEN, fn[0]); 
    
    // Calculate initial Hamiltonian
	iters += ks_congrad( 
        F_OFFSET(phi), F_OFFSET(xxx), mass, 
		niter, nrestart, rsqmin, MILC_PRECISION, 
		EVEN, &final_rsq, fn[0] 
    );
    startaction = d_action();

    // Save backup gauge field
	gauge_field_copy(F_OFFSET(link[0]), F_OFFSET(old_link[0]));

    // MILC's convention for a 2nd-order Omelyan integrator
    for(step=2; step <= steps; step+=2){
        if (step == 2){update_u(0.5*lmbda*epsilon);}
        else {update_u(lmbda*epsilon);}
        update_v(step,0.333);
        update_u((2.0-lmbda)*epsilon);
        update_v(step,0.667);
        if (step == steps){update_u(0.5*lmbda*epsilon);}
    }

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