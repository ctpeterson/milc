/******************  rand_ahmat.c  (in su3.a) ***************************
*									*
* void random_anti_hermitian( anti_hermitmat *mat_antihermit, passthru *prn_pt)*
* Creates gaussian random anti-hermitian matrices			*
* Normalization is < |m01|^2 > = 1, or < m01.real*m01.real > = 1/2	*
* The argument "prn_pt" is a pointer to be passed to gaussian_rand_no() *
* RS6000 may choke on void *						*
*/
#include "../include/config.h"
#include <math.h>
#include "../include/complex.h"
#include "../include/su3.h"
#include <stdio.h>

//#ifdef OLD_GAUSSRAND

void random_anti_hermitian( anti_hermitmat *mat_antihermit, double_prn *prn_pt) {
Real r3,r8;
Real sqrt_third;

	//printf("YES LEGACY\n");

	sqrt_third = sqrt( (double)(1.0/3.0) );
        r3=gaussian_rand_no(prn_pt);
	r8=gaussian_rand_no(prn_pt);
	mat_antihermit->m00im=r3+sqrt_third*r8;
	mat_antihermit->m11im= -r3+sqrt_third*r8;
	mat_antihermit->m22im= -2.0*sqrt_third*r8;
	mat_antihermit->m01.real=gaussian_rand_no(prn_pt);
	mat_antihermit->m02.real=gaussian_rand_no(prn_pt);
	mat_antihermit->m12.real=gaussian_rand_no(prn_pt);
	mat_antihermit->m01.imag=gaussian_rand_no(prn_pt);
	mat_antihermit->m02.imag=gaussian_rand_no(prn_pt);
	mat_antihermit->m12.imag=gaussian_rand_no(prn_pt);
	
	
	//if (this_node == 0)
	//{
	/*
	printf(
	  "%0.17lf %0.17lf %0.17lf %0.17lf %0.17lf %0.17lf \n",
	  mat_antihermit->m01.real,
	  mat_antihermit->m02.real,
	  mat_antihermit->m12.real,
	  mat_antihermit->m01.imag,
	  mat_antihermit->m02.imag,
	  mat_antihermit->m12.imag
	);
	*/
	/*
	printf("%0.17lf %0.17lf \n", mat_antihermit->m01.real,mat_antihermit->m01.imag);
	printf("%0.17lf %0.17lf \n", -mat_antihermit->m01.real,mat_antihermit->m01.imag);
	printf("%0.17lf %0.17lf \n", mat_antihermit->m02.real,mat_antihermit->m02.imag);
	printf("%0.17lf %0.17lf \n", -mat_antihermit->m02.real,mat_antihermit->m02.imag);
	printf("%0.17lf %0.17lf \n", mat_antihermit->m12.real,mat_antihermit->m12.imag);
	printf("%0.17lf %0.17lf \n", -mat_antihermit->m12.real,mat_antihermit->m12.imag);
	//}
	exit(0);
	*/
	

}/*random_anti_hermitian_*/

/*
#else

void random_anti_hermitian( anti_hermitmat *mat_antihermit, double_prn *prn_pt) {
complex r38;
Real sqrt_third;

	//printf("NOT LEGACY\n");

	sqrt_third = sqrt( (double)(1.0/3.0) );
	r38 = complex_gaussian_rand_no(prn_pt);
	mat_antihermit->m00im=r38.real+sqrt_third*r38.imag;
	mat_antihermit->m11im= -r38.real+sqrt_third*r38.imag;
	mat_antihermit->m22im= -2.0*sqrt_third*r38.imag;
	mat_antihermit->m01=complex_gaussian_rand_no(prn_pt);
	mat_antihermit->m02=complex_gaussian_rand_no(prn_pt);
	mat_antihermit->m12=complex_gaussian_rand_no(prn_pt);

}// random_anti_hermitian_

#endif
*/
