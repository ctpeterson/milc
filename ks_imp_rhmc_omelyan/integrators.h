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
 * @brief Defines data structures and methods for molecular dynamics
 * @author Curtis Taylor Peterson <curtistaylorpetersonwork@gmail.com>
 * @details
 * - Sexton:1992nu
 *   - Title: "Hamiltonian evolution for the hybrid Monte Carlo algorithm"
 *   - J.C. Sexton, D.H. Weingarten
 *   - doi: 0550-3213(92)90263-B
 * - Omelyan:2003bjg
 *   - Title: "Symplectic analytically integrable decomposition algorithms:
 *     classification, derivation, and application to molecular dynamics,
 *     quantum and celestial mechanics simulations"
 *   - I.P. Omelyan, I.M. Mryglod, R. Folk
 *   - doi: S0010-4655(02)00754-3
 */

#ifndef INTEGRATORS_H
#define INTEGRATORS_H

#include "ks_imp_includes.h"
#include <stdbool.h>

typedef enum {
  /**
   * @typedef @enum Update
   * @brief Enumerate supported integrator updates
   */
  UpdateT, /**< Gauge (a.k.a., "position" or "time") update */
  UpdateV  /**< Momentum update */
} Update;

typedef enum
{
  /** 
   * @typedef @enum BaseIntegrator
   * @brief Enumerate supported integrators
   * @details
   */
  LeapFrog, /**< Standard TVT leapfrog */
  SextonWeingarten, /**< Sexton:1992nu, Eqn. (6.1) */
  Omelyan2MN, /**< Omelyan:2003bjg, Eqn. (31) */
  Omelyan4MN4FP, /**< Omelyan:2003bjg, Eqn. (58) & (62) */
  Omelyan4MN5FP /**< Omelyan:2003bjg, Eqn. (72) & (80) */
} BaseIntegrator;

#ifndef FERMION_INT_ALG
#define FERMION_INT_ALG Omelyan2MN
#endif 

#ifndef GAUGE_INT_ALG
#define GAUGE_INT_ALG Omelyan2MN
#endif 

typedef struct {
  /**
   * @typedef @struct UpdateAtom
   * @brief Single-update data structure
   */
  Real stepSize; /**< Update step size */
  Update update; /**< Update type (Update) */
} UpdateAtom;

typedef struct
{
  /**
   * @typedef @struct UpdateMolecule
   * @brief List of single-update data structures (UpdateAtom)
   */
  UpdateAtom *updates; /**< "List" of updates */
  size_t size,capacity; /**< "List" properties */
} UpdateMolecule;

typedef struct 
{
  /**
   * @typedef @struct UpdateProtein
   * @brief Represents sequence of updates (UpdateMolecule)
   */
  int steps; /**< Total number of steps */
  Real lmbda,rho,theta,vartheta; /**< Integrator parameters (Omelyan:2003bjg) */
  Real stepSize; /**< trajectoryLength/steps */
  UpdateMolecule molecule; /**< "List" of updates at each step */
} UpdateProtein;

typedef struct
{
  /**
   * @typedef @struct UpdateOrganism
   * @brief Represents full symplectic integrator (UpdateProtein)
   */
  bool nested; /**< Flag telling integrator if gauge update nested */
  UpdateProtein gauge; /**< Gauge integrator properties */
  UpdateProtein fermion; /**< Fermion integrator properties */
} UpdateOrganism;

char integratorToChar(BaseIntegrator integrator)
{
  switch (integrator)
  {
    case LeapFrog: return "leap frog"; break;
    case Omelyan2MN: return "2nd-order Omelyan"; break;
    case Omelyan4MN4FP: return "4th-order Omelyan w/ 4 force evals."; break;
    case Omelyan4MN5FP: return "4th-order Omelyan w/ 5 force evals."; break;
  }
}

UpdateAtom newUpdateAtom(Real stepSize, Update update)
{
  /**
   * @brief UpdateAtom constructor
   * @param stepSize Update step size
   * @param update Update type
   * @return UpdateAtom struct
   */
  UpdateAtom atom;
  atom.stepSize = stepSize;
  atom.update = update;
  return atom;
}

void initUpdateMolecule(UpdateMolecule *molecule, size_t capacity)
{
  /**
   * @brief Initialize "list" of integrator updates
   * @return UpdateMolecule struct
   */
  molecule->updates = (UpdateAtom *)malloc(capacity*sizeof(UpdateAtom));
  molecule->size = capacity;
  molecule->capacity = capacity;
}

void assign(UpdateMolecule *molecule, int index, Update update, Real stepSize)
{
  /**
   * @brief Sets element of UpdateMolecule type to UpdateAtom type at "index"
   * @param molecule Pointer to "list" of integrator updates
   * @param index Location of assignment
   * @param update Update type
   * @param stepSize Update step size
   */
  UpdateAtom updateAtom = newUpdateAtom(stepSize,update);
  molecule->updates[index] = updateAtom;
}

int size(UpdateMolecule *molecule)
{
  /**
   * @brief Gets number of updates in "list" of updates
   * @param molecule "List" of integrator updates
   * @return Number of updates in "list" of updates
   */
  return molecule->size;
}

float getStepSize(UpdateMolecule *molecule, int update)
{
  /**
   * @brief Grab update step size from "list" of integrator updates
   * @param molecule "List" of integrator updates
   * @param int Index of integrator update in "list" of integrator updates
   * @return Update step size from "list" of integrator updates
   */
  return molecule->updates[update].stepSize;
}

Update getUpdate(UpdateMolecule *molecule, int update)
{
  /**
   * @brief Grab update type (Update) from "list" of integrator updates
   * @param molecule "List" of integrator updates
   * @param int Index of integrator update in "list" of integrator updates
   * @return Update enum member
   */
  return molecule->updates[update].update;
}

void freeUpdateMolecule(UpdateMolecule *molecule)
{
  /**
   * @brief Frees "list" of updates in molecule from memory
   * @param molecule Pointer to "list" of integrator updates
   */
  free(molecule->updates);
  molecule->updates = NULL;
  molecule->size = molecule->capacity = 0;
}

// Used to access "freeUpdateMolecule"
void freeUpdateProtein(UpdateProtein *protein)
{freeUpdateMolecule(&protein->molecule);}
void freeUpdateOrganism(UpdateOrganism *organism)
{freeUpdateProtein(&organism->fermion);freeUpdateProtein(&organism->gauge);}

void newUpdateMolecule(UpdateProtein *protein, BaseIntegrator scheme){
  /**
   * @brief Create "list" of integrator updates
   * @param protein Represents sequence of "lists" of integrator updates
   * @param scheme Type of symplectic integrator
   */
  Real lmbda = protein->lmbda;
  Real rho = protein->rho;
  Real theta = protein->theta;
  Real vartheta = protein->vartheta;
  Real stepSize = protein->stepSize;
  switch (scheme)
  {
    case LeapFrog:
      initUpdateMolecule(&protein->molecule,3);
      assign(&protein->molecule,0,UpdateT,0.5*stepSize);
      assign(&protein->molecule,1,UpdateV,stepSize);
      assign(&protein->molecule,2,UpdateT,0.5*stepSize);
      break;
    case Omelyan2MN:
    case SextonWeingarten:
      initUpdateMolecule(&protein->molecule,5);
      assign(&protein->molecule,0,UpdateT,lmbda*stepSize);
      assign(&protein->molecule,1,UpdateV,0.5*stepSize);
      assign(&protein->molecule,2,UpdateT,(1.0-2.0*lmbda)*stepSize);
      assign(&protein->molecule,3,UpdateV,0.5*stepSize);
      assign(&protein->molecule,4,UpdateT,lmbda*stepSize);
      break;
    case Omelyan4MN4FP:
      initUpdateMolecule(&protein->molecule,9);
      assign(&protein->molecule,0,UpdateT,rho*stepSize);
      assign(&protein->molecule,1,UpdateV,lmbda*stepSize);
      assign(&protein->molecule,2,UpdateT,theta*stepSize);
      assign(&protein->molecule,3,UpdateV,(0.5-lmbda)*stepSize);
      assign(&protein->molecule,4,UpdateT,(1.0-2.0*theta-2.0*rho)*stepSize);
      assign(&protein->molecule,5,UpdateV,(0.5-lmbda)*stepSize);
      assign(&protein->molecule,6,UpdateT,theta*stepSize);
      assign(&protein->molecule,7,UpdateV,lmbda*stepSize);
      assign(&protein->molecule,8,UpdateT,rho*stepSize);
      break;
    case Omelyan4MN5FP:
      initUpdateMolecule(&protein->molecule,11);
      assign(&protein->molecule,0,UpdateT,rho*stepSize);
      assign(&protein->molecule,1,UpdateV,vartheta*stepSize);
      assign(&protein->molecule,2,UpdateT,theta*stepSize);
      assign(&protein->molecule,3,UpdateV,lmbda*stepSize);
      assign(&protein->molecule,4,UpdateT,(0.5-theta-rho)*stepSize);
      assign(&protein->molecule,5,UpdateV,(1.0-2.0*lmbda-2.0*vartheta)*stepSize);
      assign(&protein->molecule,6,UpdateT,(0.5-theta-rho)*stepSize);
      assign(&protein->molecule,7,UpdateV,lmbda*stepSize);
      assign(&protein->molecule,8,UpdateT,theta*stepSize);
      assign(&protein->molecule,9,UpdateV,vartheta*stepSize);
      assign(&protein->molecule,10,UpdateT,rho*stepSize);
      break;
    default: break;
  }
}

void newUpdateProtein(
  UpdateProtein *protein,
  BaseIntegrator scheme,
  Real stepSize,
  int steps
){
  /**
   * @brief Creates single integrator for either gauge for fermion field
   * @param scheme Type of symplectic integrator
   * @param stepSize "trajectory length"/"total number of steps"
   * @param steps Number of steps per call (of integrator)
   * @return UpdateProtein struct
   */
  protein->stepSize = stepSize;
  protein->steps = steps;
  node0_printf("Integrator:\n");
  switch (scheme)
  {
    case LeapFrog: break;
    case SextonWeingarten: 
      node0_printf("  -leapfrog\n");
      protein->lmbda = 1.0/6.0; break;
    case Omelyan2MN: 
      node0_printf("  -2nd-order Omelyan\n");
      protein->lmbda = 0.1931833275037836; break;
    case Omelyan4MN4FP:
      node0_printf("  -4th-order Omelyan w/ 4 force evaluations\n");
      protein->rho = 0.1786178958448091;
      protein->theta = -0.06626458266981843;
      protein->lmbda = 0.7123418310626056;
      break;
    case Omelyan4MN5FP:
      node0_printf("  -4th-order Omelyan w/ 5 force evaluations\n");
      protein->rho = 0.2750081212332419;
      protein->theta = -0.1347950099106792;
      protein->vartheta = -0.08442961950707149;
      protein->lmbda = 0.3549000571574260;
      break;
  }
  node0_printf("  -steps = %i\n",protein->steps);
  newUpdateMolecule(protein,scheme);
}

void newIntegrator(
  UpdateOrganism *organism,
  Real trajectoryLength,
  BaseIntegrator fermionIntegrator,
  BaseIntegrator gaugeIntegrator,
  int fermSteps, 
  int gaugeSteps
){
  /**
   * @brief Creates full integrator for gauge & fermion fields
   * @param trajectoryLength Integrator trajectory length
   * @param fermionIntegrator Type of integrator for fermion fields,
   *        both fermion & Hasenbusch
   * @param gaugeIntegrator Type of integrator for gauge field
   * @param fermSteps Number of fermion steps per call to full integrator
   * @param gaugeSteps Number of gauge steps per fermion (outer) gauge update
   *        if integrator is nested (gaugeSteps > 0). For non-nested integrator,
   *        choose gaugeSteps = -1
   * @return UpdateOrganism struct
   */
  Real fermStepSize = trajectoryLength/(Real)fermSteps;
  Real gaugeStepSize = 1.0/(Real)gaugeSteps;
  newUpdateProtein(&organism->fermion,fermionIntegrator,fermStepSize,fermSteps);
  newUpdateProtein(&organism->gauge,gaugeIntegrator,gaugeStepSize,gaugeSteps);
  organism->nested = !((fermionIntegrator == gaugeIntegrator) && (gaugeSteps == -1));
}

void gaugeUpdateInternal(UpdateOrganism *organism, float outerStepSize)
{
  /**
   * @brief Runs through one step of inner "nested" gauge integrator
   * @param organism Data structure representing full integrator
   * @param outerStepSize "Trajectory length" for this nested gauge update
   * @details
   * Runs through a full cycle (= "updates") of one step of the nested 
   * gauge integrator. Value for each step size is the outer trajectory 
   * length ("outerStepSize") multiplied by the step size for the nested 
   * gauge integrator (= "integrator step factor" x 1/"number of internal 
   * gauge steps")
   */
  for (int update = 0; update < size(&organism->gauge.molecule); update++)
  {
    float innerStepSize = getStepSize(&organism->gauge.molecule,update);
    switch (getUpdate(&organism->gauge.molecule,update))
    {
      case UpdateT: update_u(innerStepSize*outerStepSize); break;
      case UpdateV: update_h_gauge(innerStepSize*outerStepSize); break;
    }
  }
}

int fermionUpdateNested(UpdateOrganism *organism, su3_vector **multi_x)
{
  /**
   * @brief Runs through one step of "outer" fermion integrator
   * @param organism Data structure representing full integrator
   * @param multi_x List of vectors for multi-shift solver
   * @details
   * Runs through a full cycle (= "updates") of outer fermion integrator.
   * At each gauge update, a full cycle through the gauge integrator
   * is performed.
   * @return Number of CG iterations from fermion/Hasenbusch solve
   */
  int iters = 0;
  for (int update = 0; update < size(&organism->fermion.molecule); update++)
  {
    float stepSize = getStepSize(&organism->fermion.molecule,update);
    switch (getUpdate(&organism->fermion.molecule,update))
    {
      case UpdateT:
        for (int step = 0; step < organism->gauge.steps; step++)
        {gaugeUpdateInternal(organism,stepSize);}
        break;
      case UpdateV: iters += update_h_fermion(stepSize,multi_x); break;
    }
  }
  return iters;
}

int collectiveUpdate(UpdateOrganism *organism, su3_vector **multi_x)
{
  /**
   * @brief Runs through one step of gauge + fermion integrator
   * @param organism Data structure representing full integrator
   * @param multi_x List of vectors for multi-shift solver
   * @return Number of CG iterations from fermion/Hasenbusch solve
   */
  int iters = 0;
  int updates = size(&organism->fermion.molecule);
  UpdateMolecule molecule = organism->fermion.molecule;
  for (int update = 0; update < updates; update++)
  {
    float stepSize = getStepSize(&organism->fermion.molecule,update);
    switch (getUpdate(&organism->fermion.molecule,update))
    {
      case UpdateT: update_u(stepSize); break;
      case UpdateV: iters += update_h_rhmc(stepSize,multi_x); break;
    }
  }
  return iters;
}

int integrate(UpdateOrganism *organism, su3_vector **multi_x)
{
  /**
   * @brief Constructs full integrator & runs molecular dynamics update with it
   * @param organism Full gauge + fermion integrator
   * @param multi_x List of vectors for multi-shift solver
   * @return Number of CG iterations from fermion/Hasenbusch solve
   */
  int iters = 0;
  int steps = organism->fermion.steps;
  for (int step = 0; step < steps; step++)
  {
    switch (organism->nested)
    {
      case true: iters += fermionUpdateNested(organism, multi_x); break;
      case false: iters += collectiveUpdate(organism, multi_x); break;
    }
  }
  return iters;
}

#endif