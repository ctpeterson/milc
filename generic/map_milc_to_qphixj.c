/*********************** map_milc_to_qphixjmilc.c **********************/
/* Functions for mapping MILC data layouts to raw QPHIXJ layouts       */
/* C. DeTar 10/19/2005
   C. DeTar 4/29/07    add DiracFermion = wilson_fermion
*/

#include "generic_includes.h"
#include "../include/generic_qphixj.h"
#include "../include/generic_qphixjmilc.h"
#include "../include/openmp_defs.h"

/* Copy types with possible conversion */

/* Convert (or copy) MILC types between specific and prevailing precision */

#if (PRECISION==1)

static void 
f2p_mat(su3_matrix *dest, fsu3_matrix *src){
  memcpy((void *)dest, (void *)src, sizeof(fsu3_matrix));
}

static void 
p2f_mat(fsu3_matrix *dest, su3_matrix *src){
  memcpy((void *)dest, (void *)src, sizeof(fsu3_matrix));
}

static void 
d2p_mat(su3_matrix *dest, dsu3_matrix *src){
  int i,j;
  
  for(i = 0; i < 3; i++)for(j = 0; j < 3; j++){
    dest->e[i][j].real = src->e[i][j].real;
    dest->e[i][j].imag = src->e[i][j].imag;
  }
}

static void 
p2d_mat(dsu3_matrix *dest, su3_matrix *src){
  int i,j;
  
  for(i = 0; i < 3; i++)for(j = 0; j < 3; j++){
    dest->e[i][j].real = src->e[i][j].real;
    dest->e[i][j].imag = src->e[i][j].imag;
  }
}

static void 
f2p_vec(su3_vector *dest, fsu3_vector *src){
  memcpy((void *)dest, (void *)src, sizeof(fsu3_vector));
}

static void 
p2f_vec(fsu3_vector *dest, su3_vector *src){
  memcpy((void *)dest, (void *)src, sizeof(fsu3_vector));
}

static void 
d2p_vec(su3_vector *dest, dsu3_vector *src){
  int i;
  
  for(i = 0; i < 3; i++){
    dest->c[i].real = src->c[i].real;
    dest->c[i].imag = src->c[i].imag;
  }
}

static void 
p2d_vec(dsu3_vector *dest, su3_vector *src){
  int i;
  
  for(i = 0; i < 3; i++){
    dest->c[i].real = src->c[i].real;
    dest->c[i].imag = src->c[i].imag;
  }
}

#ifndef HAVE_NO_QPHIXJ_D

static void 
f2p_wvec(wilson_vector *dest, fwilson_vector *src){
  memcpy((void *)dest, (void *)src, sizeof(fwilson_vector));
}

static void 
p2f_wvec(fwilson_vector *dest, wilson_vector *src){
  memcpy((void *)dest, (void *)src, sizeof(fwilson_vector));
}

static void 
d2p_wvec(wilson_vector *dest, dwilson_vector *src){
  int i,j;
  
  for(j = 0; j < 4; j++)for(i = 0; i < 3; i++){
    dest->d[j].c[i].real = src->d[j].c[i].real;
    dest->d[j].c[i].imag = src->d[j].c[i].imag;
  }
}

static void 
p2d_wvec(dwilson_vector *dest, wilson_vector *src){
  int i,j;
  
  for(j = 0; j < 4; j++)for(i = 0; i < 3; i++){
    dest->d[j].c[i].real = src->d[j].c[i].real;
    dest->d[j].c[i].imag = src->d[j].c[i].imag;
  }
}

#endif

#else

static void 
f2p_mat(su3_matrix *dest, fsu3_matrix *src){
  int i,j;
  
  for(i = 0; i < 3; i++)for(j = 0; j < 3; j++){
    dest->e[i][j].real = src->e[i][j].real;
    dest->e[i][j].imag = src->e[i][j].imag;
  }
}

/* Convert (or copy) su3_matrix from prevailing to single precision */
static void 
p2f_mat(fsu3_matrix *dest, su3_matrix *src){
  int i,j;
  
  for(i = 0; i < 3; i++)for(j = 0; j < 3; j++){
    dest->e[i][j].real = src->e[i][j].real;
    dest->e[i][j].imag = src->e[i][j].imag;
  }
}

static void 
d2p_mat(su3_matrix *dest, dsu3_matrix *src){
  memcpy((void *)dest, (void *)src, sizeof(dsu3_matrix));
}

static void 
p2d_mat(dsu3_matrix *dest, su3_matrix *src){
  memcpy((void *)dest, (void *)src, sizeof(dsu3_matrix));
}

static void 
f2p_vec(su3_vector *dest, fsu3_vector *src){
  int i;
  
  for(i = 0; i < 3; i++){
    dest->c[i].real = src->c[i].real;
    dest->c[i].imag = src->c[i].imag;
  }
}

/* Convert (or copy) su3_vector from prevailing to single precision */
static void 
p2f_vec(fsu3_vector *dest, su3_vector *src){
  int i;
  
  for(i = 0; i < 3; i++){
    dest->c[i].real = src->c[i].real;
    dest->c[i].imag = src->c[i].imag;
  }
}

static void 
d2p_vec(su3_vector *dest, dsu3_vector *src){
  memcpy((void *)dest, (void *)src, sizeof(dsu3_vector));
}

static void 
p2d_vec(dsu3_vector *dest, su3_vector *src){
  memcpy((void *)dest, (void *)src, sizeof(dsu3_vector));
}

#ifndef HAVE_NO_QPHIXJ_D

static void 
f2p_wvec(wilson_vector *dest, fwilson_vector *src){
  int i,j;
  
  for(j = 0; j < 4; j++)for(i = 0; i < 3; i++){
    dest->d[j].c[i].real = src->d[j].c[i].real;
    dest->d[j].c[i].imag = src->d[j].c[i].imag;
  }
}

/* Convert (or copy) wilson_vector from prevailing to single precision */
static void 
p2f_wvec(fwilson_vector *dest, wilson_vector *src){
  int i,j;
  
  for(j = 0; j < 4; j++)for(i = 0; i < 3; i++){
    dest->d[j].c[i].real = src->d[j].c[i].real;
    dest->d[j].c[i].imag = src->d[j].c[i].imag;
  }
}

static void 
d2p_wvec(wilson_vector *dest, dwilson_vector *src){
  memcpy((void *)dest, (void *)src, sizeof(dwilson_vector));
}

static void 
p2d_wvec(dwilson_vector *dest, wilson_vector *src){
  memcpy((void *)dest, (void *)src, sizeof(dwilson_vector));
}

#endif
#endif

/********************************************************************/
#include "map_milc_to_qphixj_all.c"
/********************************************************************/

/********************************************************************/
/* Procedures involving MILC and raw fields                         */
/********************************************************************/


/* Hash tables for optimizing data mapping */

make_create_hash(G);
make_create_hash(F);
make_create_hash(V);
make_create_hash(D);


/* Storage for raw gauge field */

make_create_raw4(F, G, fsu3_matrix);
make_create_raw4(D, G, dsu3_matrix);

make_destroy_raw4(F, G, fsu3_matrix);
make_destroy_raw4(D, G, dsu3_matrix);

#if 0

/* Storage for raw gauge momentum */

make_create_raw4(F, F, fsu3_matrix);
make_create_raw4(D, F, dsu3_matrix);

make_destroy_raw4(F, F, fsu3_matrix);
make_destroy_raw4(D, F, dsu3_matrix);

#endif

/* Storage for raw su3 vector field */

make_create_raw(F, V, fsu3_vector);
make_create_raw(D, V, dsu3_vector);

make_destroy_raw(F, V, fsu3_vector);
make_destroy_raw(D, V, dsu3_vector);

/* Storage for raw wilson vector field */

make_create_raw(F, D, fwilson_vector);
make_create_raw(D, D, dwilson_vector);

make_destroy_raw(F, D, fwilson_vector);
make_destroy_raw(D, D, dwilson_vector);

/* Map gauge field from site to raw */

make_create_raw4_from_site(F, G, fsu3_matrix, su3_matrix);
make_create_raw4_from_site(D, G, dsu3_matrix, su3_matrix);

/* Map gauge field from field to raw */

make_create_raw4_from_field(F, G, fsu3_matrix, su3_matrix);
make_create_raw4_from_field(D, G, dsu3_matrix, su3_matrix);

#if 0

/* Map gauge momentum from site to raw */

make_create_raw4_from_site(F, F, fsu3_matrix, anti_hermitmat);
make_create_raw4_from_site(D, F, dsu3_matrix, anti_hermitmat);

/* Map gauge momentum from field to raw */

make_create_raw4_from_field(F, F, fsu3_matrix, anti_hermitmat);
make_create_raw4_from_field(D, F, dsu3_matrix, anti_hermitmat);

#endif

/* Map color vector from site to raw */

make_create_raw_from_site(F, V, fsu3_vector, su3_vector);
make_create_raw_from_site(D, V, dsu3_vector, su3_vector);

/* Map color vector from field to raw */

make_create_raw_from_field(F, V, fsu3_vector, su3_vector);
make_create_raw_from_field(D, V, dsu3_vector, su3_vector);


/* Map wilson vector from site to raw */

make_create_raw_from_site(F, D, fwilson_vector, wilson_vector);
make_create_raw_from_site(D, D, dwilson_vector, wilson_vector);

/* Map wilson vector from field to raw */

make_create_raw_from_field(F, D, fwilson_vector, wilson_vector);
make_create_raw_from_field(D, D, dwilson_vector, wilson_vector);

/* Map gauge field from raw to site */

make_unload_raw4_to_site(F, G, su3_matrix, fsu3_matrix);
make_unload_raw4_to_site(D, G, su3_matrix, dsu3_matrix);

/* Map gauge field from raw to field */

make_unload_raw4_to_field(F, G, su3_matrix, fsu3_matrix);
make_unload_raw4_to_field(D, G, su3_matrix, dsu3_matrix);

#if 0

/* Map gauge momentum from raw to site */

make_unload_raw4_to_site(F, F, anti_hermitmat, fsu3_matrix);
make_unload_raw4_to_site(D, F, anti_hermitmat, dsu3_matrix);

/* Map gauge momentum from raw to field */

make_unload_raw4_to_field(F, F, anti_hermitmat, fsu3_matrix);
make_unload_raw4_to_field(D, F, anti_hermitmat, dsu3_matrix);

#endif

/* Map color vector from raw to site */

make_unload_raw_to_site(F, V, su3_vector, fsu3_vector);
make_unload_raw_to_site(D, V, su3_vector, dsu3_vector);

/* Map color vector from raw to field */

make_unload_raw_to_field(F, V, su3_vector, fsu3_vector);
make_unload_raw_to_field(D, V, su3_vector, dsu3_vector);

#ifndef HAVE_NO_QPHIXJ_D

/* Map wilson vector from raw to site */

make_unload_raw_to_site(F, D, wilson_vector, fwilson_vector);
make_unload_raw_to_site(D, D, wilson_vector, dwilson_vector);

/* Map wilson vector from raw to field */

make_unload_raw_to_field(F, D, wilson_vector, fwilson_vector);
make_unload_raw_to_field(D, D, wilson_vector, dwilson_vector);

#endif

/********************************************************************/
/* Procedures involving MILC and QPHIXJ fields                         */
/********************************************************************/

#if 0

make_create_from_site4(F, G, QPHIXJ_F3_GaugeField, fsu3_matrix, float);
make_create_from_site4(D, G, QPHIXJ_D3_GaugeField, dsu3_matrix, double);

/* Map gauge field from field to QPHIXJ */

make_create_from_field4(F, G, QPHIXJ_F3_GaugeField, fsu3_matrix, su3_matrix, float);
make_create_from_field4(D, G, QPHIXJ_D3_GaugeField, dsu3_matrix, su3_matrix, double);

/* Map gauge momentum from site to QPHIXJ */

make_create_from_site4(F, F, QPHIXJ_F3_Force, fsu3_matrix, float);
make_create_from_site4(D, F, QPHIXJ_D3_Force, dsu3_matrix, double);

/* Map gauge momentum from field to QPHIXJ */

make_create_from_field4(F, F, QPHIXJ_F3_Force, fsu3_matrix, anti_hermitmat, float);
make_create_from_field4(D, F, QPHIXJ_D3_Force, dsu3_matrix, anti_hermitmat, double);

/* Map color vector from site to QPHIXJ */

make_create_from_site(F, V, QPHIXJ_F3_ColorVector, fsu3_vector, float);
make_create_from_site(D, V, QPHIXJ_D3_ColorVector, dsu3_vector, double);

/* Map color vector from field to QPHIXJ */

make_create_from_field(F, V, QPHIXJ_F3_ColorVector, fsu3_vector, su3_vector, float);
make_create_from_field(D, V, QPHIXJ_D3_ColorVector, dsu3_vector, su3_vector, double);

/* Map wilson vector from site to QPHIXJ */

make_create_from_site(F, D, QPHIXJ_F3_DiracFermion, fwilson_vector, float);
make_create_from_site(D, D, QPHIXJ_D3_DiracFermion, dwilson_vector, double);

#endif

/* Map wilson vector from field to QPHIXJ */

make_create_from_field(F, D, QPHIXJ_F3_DiracFermion, fwilson_vector, wilson_vector, float);
make_create_from_field(D, D, QPHIXJ_D3_DiracFermion, dwilson_vector, wilson_vector, double);

#if 0
/* Map QPHIXJ color vector field to site */

make_unload_to_site(F, V, QPHIXJ_F3_ColorVector, fsu3_vector, float);
make_unload_to_site(D, V, QPHIXJ_D3_ColorVector, dsu3_vector, double);

#endif

/* Map QPHIXJ Dirac fermion field to field */

make_unload_to_field(F, D, QPHIXJ_F3_DiracFermion, fwilson_vector, wilson_vector, float);
make_unload_to_field(D, D, QPHIXJ_D3_DiracFermion, dwilson_vector, wilson_vector, double);


#if 0

/* Map QPHIXJ wilson vector field to site */

make_unload_to_site(F, D, QPHIXJ_F3_DiracFermion, fwilson_vector, float);
make_unload_to_site(D, D, QPHIXJ_D3_DiracFermion, dwilson_vector, double);

/* Map QPHIXJ wilson vector field to field */

make_unload_to_field(F, D, QPHIXJ_F3_DiracFermion, fwilson_vector, wilson_vector, float);
make_unload_to_field(D, D, QPHIXJ_D3_DiracFermion, dwilson_vector, wilson_vector, double);

/* Map QPHIXJ gauge force to site */

make_unload_to_site4(F, F, QPHIXJ_F3_Force, fsu3_matrix, float);
make_unload_to_site4(D, F, QPHIXJ_D3_Force, dsu3_matrix, double);

/* Map QPHIXJ gauge force to field */

make_unload_to_field4(F, F, QPHIXJ_F3_Force, fsu3_matrix, anti_hermitmat, float);
make_unload_to_field4(D, F, QPHIXJ_D3_Force, dsu3_matrix, anti_hermitmat, double);

/* Map QPHIXJ gauge field to site */

make_unload_to_site4(F, G, QPHIXJ_F3_GaugeField, fsu3_matrix, float);
make_unload_to_site4(D, G, QPHIXJ_D3_GaugeField, dsu3_matrix, double);

/* Map QPHIXJ gauge field to field */

make_unload_to_field4(F, G, QPHIXJ_F3_GaugeField, fsu3_matrix, su3_matrix, float);
make_unload_to_field4(D, G, QPHIXJ_D3_GaugeField, dsu3_matrix, su3_matrix, double);


/* Map MILC gauge links in site plus raw clover terms to QPHIX fermion links */

make_create_L_from_sites(F, QPHIXJ_F3_FermionLinksWilson, fsu3_matrix, float);
make_create_L_from_sites(D, QPHIXJ_D3_FermionLinksWilson, dsu3_matrix, double);

#endif

/* Map MILC gauge links in site plus raw backward links plus raw clover terms to QPHIX fermion links */

make_create_L_from_fieldback_and_sites(F, QPHIXJ_F3_FermionLinksWilson, fsu3_matrix, float);
make_create_L_from_fieldback_and_sites(D, QPHIXJ_D3_FermionLinksWilson, dsu3_matrix, double);

#if 0

/* Map preconstructed MILC fat and long link fields to QPHIXJ fermion links */

make_create_L_from_fields(F, QPHIXJ_F3_FermionLinksAsqtad, fsu3_matrix, su3_matrix, float);
make_create_L_from_fields(D, QPHIXJ_D3_FermionLinksAsqtad, dsu3_matrix, su3_matrix, double);

#endif

#ifndef HAVE_NO_QPHIXJ_D

/* Unload fat and long links from QPHIXJ fermion links to MILC fields */

#if 0
make_unload_L_to_fields(F, QPHIXJ_F3_FermionLinksAsqtad, fsu3_matrix, su3_matrix, float );

make_unload_L_to_fields(D, QPHIXJ_D3_FermionLinksAsqtad, dsu3_matrix, su3_matrix, double );
#endif

#endif