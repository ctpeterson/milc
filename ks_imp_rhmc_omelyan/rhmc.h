/** 
 * Copyright (c) 1998-2018 MILC Developers and Contributors

 * Permission is hereby granted, free of charge, to any person obtaining
 * a copy of this software and associated documentation files (the
 * "Software"), to deal in the Software without restriction, including
 * without limitation the rights to use, copy, modify, merge, publish,
 * distribute, sublicense, and/or sell copies of the Software, and to
 * permit persons to whom the Software is furnished to do so, subject to
 * the following conditions:

 * The above copyright notice, this permission notice and the following
 * disclaimers shall be included in all copies or substantial portions of
 * the Software.

 * Neither the name of the MILC collaboration, nor the names of its
 * developers or contributors may be used to endorse or promote products
 * derived from this Software without specific prior written permission.

 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
 * EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
 * MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
 * NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE
 * LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION
 * OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION
 * WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE. 
 */
/**
 * @brief Defines data structures and methods for (Rational)HMC sampling
 * @author Curtis Taylor Peterson <curtistaylorpetersonwork@gmail.com>
 */

#include "ks_imp_includes.h"
#include "integrators.h"

#if defined(USE_FF_GPU)
#include "../include/generic_quda.h"
#define special_alloc qudaAllocateManaged
#define special_free qudaFreeManaged
#else
#define special_alloc malloc
#define special_free free
#endif

typedef enum 
{
  /**
   * @typedef @enum FermionAction
   * @brief Enumerate supported fermion actions
   */
  ActionHISQ, /**< Standard HISQ action */
  ActionHypISQ /**< HISQ with hypercubic smearing */
} FermionAction;

typedef struct
{
  /**
   * @typedef @struct HamiltonianMonteCarlo
   * @brief Data structure representing Hamiltonian Monte Carlo
   */
  int nnaiks; /**< Number of outer naik terms */
  int nmultix; /**< Number of multi-shift solution vectors */
  int cgIters; /**< Number of multi-shift solution vectors */
  Real hi,hf; /**< Initial & final value of Hamiltonian */
  su3_vector **multi_x; /**< Mult-shift solution vectors */
  su3_vector *sumvec; /**< Auxiliary SU(3) vector */
  imp_ferm_links_t** fn; /**< Smeared fermion links */
  FermionAction fermionAction; /**< Fermion action */
  UpdateOrganism integrator; /**< Integrator for fermion & gauge */
} HamiltonianMonteCarlo;

void newHamiltonianMonteCarlo(
  HamiltonianMonteCarlo *hmc,
  FermionAction fact,
  BaseIntegrator fermionIntegrator,
  BaseIntegrator gaugeIntegrator,
  Real trajectoryLength,
  int fermionSteps, 
  int gaugeSteps
){
  /**
   * @brief Constructs Hamiltonian Monte Carlo data structure
   * @param fact Fermion action
   * @param fermionIntegrator Named integrator for outer fermion update
   * @param gaugeIntegrator Named integrator for inner gauge update
   * @param trajectoryLength Integrator trajectory length
   * @param fermionSteps Number of outer fermion steps
   * @param gaugeSteps Number of inner gauge steps per fermion gauge update
   * @return HamiltonainMonteCarlo data structure
   */
  int i,j;
  
  // Initialize a few attributes
  hmc->fermionAction = fact;
  hmc->cgIters = 0;

  // Prepare multi-shift solution vectors
  hmc->nmultix = max_rat_order;
  for(j = 0, i = 0; i < n_pseudo; i++){j += rparam[i].MD.order;}
  if(j > hmc->nmultix) hmc->nmultix = j;
  hmc->multi_x = (su3_vector **)malloc(hmc->nmultix*sizeof(su3_vector *));
  if(hmc->multi_x == NULL){printf("update: No room for multi_x\n"); terminate(1);}
  for(i = 0; i < hmc->nmultix; i++){
    hmc->multi_x[i] = (su3_vector *)special_alloc( sizeof(su3_vector)*sites_on_node );
    if(hmc->multi_x[i] == NULL){printf("update: No room for multi_x\n"); terminate(1);}
  }

  // Prepare auxiliary SU(3) vector
  hmc->sumvec = (su3_vector *)malloc( sizeof(su3_vector)*sites_on_node );
  if(hmc->sumvec == NULL){printf("update: No room for sumvec\n"); terminate(1);}

  // Prepare fermion action
  switch (hmc->fermionAction)
  {
    case ActionHISQ:
    case ActionHypISQ: hmc->nnaiks = fermion_links_get_n_naiks(fn_links); break;
  }

  // Prepare integrator
  newIntegrator(
    &hmc->integrator,
    trajectoryLength,
    fermionIntegrator,
    gaugeIntegrator,
    fermionSteps,
    gaugeSteps
  );
}

void heatbath(HamiltonianMonteCarlo *hmc)
{
  /**
   * @brief Performs heatbath for momenta and pseudofermion fields
   * @param hmc HamiltonainMonteCarlo data structure
   */
  int iphi = 0;
  ranmom();
  for (int inaik = 0; inaik < hmc->nnaiks; inaik++) 
  {
    for(int jphi = 0; jphi < n_pseudo_naik[inaik]; jphi++) 
    {
      restore_fermion_links_from_site(fn_links, prec_gr[iphi]);
      hmc->fn = get_fm_links(fn_links);
      grsource_imp_rhmc(
        F_OFFSET(phi[iphi]), 
        &(rparam[iphi].GR), 
        EVEN,
		    hmc->multi_x, 
        hmc->sumvec, 
        rsqmin_gr[iphi], 
        niter_gr[iphi],
		    prec_gr[iphi], 
        hmc->fn[inaik], 
        inaik, 
		    rparam[iphi].naik_term_epsilon
      );
      iphi++;
    }
  }
}

void metropolis(HamiltonianMonteCarlo *hmc)
{
  /**
   * @brief Performs Metropolis accept/reject step
   * @param hmc HamiltonianMonteCarlo data structure
   */
  Real xrandom;
  double dh = (double)(hmc->hf-hmc->hi);
  if(this_node == 0)xrandom = myrand(&node_prn);
  broadcast_float(&xrandom);
  switch (xrandom <= exp(-dh))
  {
    case true: 
      node0_printf("ACCEPT: delta S = %e\n",dh); 
      reunitarize_ks(); 
      break;
    case false:
      node0_printf("REJECT: delta S = %e\n",dh);
      gauge_field_copy(F_OFFSET(old_link[0]),F_OFFSET(link[0]));
      invalidate_fermion_links(fn_links);
      break;
  }
}

void molecularDynamics(HamiltonianMonteCarlo *hmc)
{
  /**
   * @brief Wraps call to integrator for molecular dynamics
   * @param hmc Hamiltonian Monte Carlo data structure
   */
  gauge_field_copy(F_OFFSET(link[0]),F_OFFSET(old_link[0]));
  hmc->hi = d_action_rhmc(hmc->multi_x,hmc->sumvec);
  hmc->cgIters += integrate(&hmc->integrator,hmc->multi_x);
  hmc->hf = d_action_rhmc(hmc->multi_x,hmc->sumvec);
}

void freeHamiltonianMonteCarlo(HamiltonianMonteCarlo *hmc)
{
  /**
   * @brief Frees memory of multi-shift vectors & auxiliary vector
   * @param hmc HamiltonianMonteCarlo data structure
   */
  for(int i = 0; i < hmc->nmultix; i++){special_free(hmc->multi_x[i]);}
  free(hmc->sumvec);
  freeUpdateOrganism(&hmc->integrator);
}