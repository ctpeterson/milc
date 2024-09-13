/****** update_h.c  -- ******************/
/* updates momentum matrices for improved action */
/* D.T. & J.H., naik term    8/96
*  D.T., fat link fermion term 5/97
*  D.T. general quark action 1/98
*  D.T. two types of quarks 3/99
*  T.D. and A.H. improved gauge updating spliced in 5/97
*
* MIMD version 7 */

#include "ks_imp_includes.h"	/* definitions files and prototypes */

#ifdef ONEMASS
void eo_fermion_force_oneterm( Real eps, Real weight, su3_vector *temp_x,
			  int prec, fermion_links_t *fl )
{
  su3_vector *xxx[1] = { temp_x };
  int nterms = 1;
  Real residues[1] = { weight };

  eo_fermion_force_multi( eps, residues, xxx, nterms, prec, fl );
  
}

void eo_fermion_force_oneterm_site( Real eps, Real weight, field_offset x_off,
				    int prec, fermion_links_t *fl)
{
  su3_vector *temp_x = create_v_field_from_site_member(x_off);
  su3_vector *xxx[1] = { temp_x };
  int nterms = 1;
  Real residues[1] = { weight };

  eo_fermion_force_multi( eps, residues, xxx, nterms, prec, fl );

  destroy_v_field(temp_x);
}
#else
void eo_fermion_force_twoterms( Real eps, Real weight1, Real weight2, 
				su3_vector *x1_off, su3_vector *x2_off,
				int prec, fermion_links_t *fl ){
  su3_vector *xxx[2] = { x1_off, x2_off };
  int nterms = 2;
  Real residues[2] = { weight1, weight2 };

  eo_fermion_force_multi( eps, residues, xxx, nterms, prec, fl );
}

void eo_fermion_force_twoterms_site( Real eps, Real weight1, Real weight2,
				     field_offset x1_off, field_offset x2_off,
				     int prec, fermion_links_t *fl)
{
  su3_vector *x1 = create_v_field_from_site_member(x1_off);
  su3_vector *x2 = create_v_field_from_site_member(x2_off);
  su3_vector *xxx[2] = { x1, x2 };
  int nterms = 2;
  Real residues[2] = { weight1, weight2 };

  eo_fermion_force_multi( eps, residues, xxx, nterms, prec, fl );

  destroy_v_field(x1);
  destroy_v_field(x2);
}
#endif

void update_h( Real eps ){
  int ff_prec = MILC_PRECISION;  /* Just use prevailing precision for now */
#ifdef FN
//    free_fn_links(&fn_links);
//    free_fn_links(&fn_links_dmdu0);
  invalidate_fermion_links(fn_links);
#endif
    /* gauge field force */
    rephase(OFF);
    imp_gauge_force(eps,F_OFFSET(mom));
    rephase(ON);
    /* fermionic force */
    /* First compute M*xxx in temporary vector xxx_odd */
    /* See long comment at end of file */
	/* The diagonal term in M doesn't matter */
    restore_fermion_links_from_site(fn_links, MILC_PRECISION);

#ifdef ONEMASS
    eo_fermion_force_oneterm_site( eps, ((Real)nflavors)/4., F_OFFSET(xxx),
				   ff_prec, fn_links );
#else
/**
   eo_fermion_force_oneterm( eps, ((Real)nflavors1)/4., F_OFFSET(xxx1),
        ff_prec, fn_links );
    eo_fermion_force_oneterm( eps, ((Real)nflavors2)/4., F_OFFSET(xxx2),
        ff_prec, fn_links );
**/
/**/
    eo_fermion_force_twoterms_site( eps, ((Real)nflavors1)/4., 
				    ((Real)nflavors2)/4., F_OFFSET(xxx1), 
				    F_OFFSET(xxx2), ff_prec, fn_links );
    /**/
#endif
} /* update_h */


