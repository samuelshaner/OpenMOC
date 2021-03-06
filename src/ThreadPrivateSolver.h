/**
 * @file ThreadPrivateSolver.h
 * @brief The ThreadPrivateSolver class.
 * @date May 28, 2013
 * @author William Boyd, MIT, Course 22 (wboyd@mit.edu)
 */


#ifndef THREADPRIVATESOLVER_H_
#define THREADPRIVATESOLVER_H_

#ifdef __cplusplus
#define _USE_MATH_DEFINES
#include <math.h>
#include <omp.h>
#include <stdlib.h>
#include "CPUSolver.h"
#endif

/** Indexing scheme for the thread private FSR scalar flux for each thread */
#define _thread_flux(tid,r,e) (_thread_flux[(tid)][(r)*_num_groups+(e)])

/** Indexing scheme for the thread private Cmfd Mesh surface currents for each 
 * thread in each FSR and energy group */
#define _thread_currents(tid,r,e) (_thread_currents[(tid)*_num_mesh_cells*8*_cmfd->getNumCmfdGroups() + (r)*_cmfd->getNumCmfdGroups() + std::min((e) / _cmfd->getCmfdGroupWidth(), _cmfd->getNumCmfdGroups()-1)])


/**
 * @class ThreadPrivateSolver ThreadPrivateSolver.h "openmoc/src/ThreadPrivateSolver.h"
 * @brief This is a subclass of the CPUSolver which uses thread private
 *        arrays for the FSR scalar fluxes to minimize OpenMPC atomics.
 * @details Since this class stores a separate copy of the FSR scalar
 *          fluxes for each OMP thread, the memory requirements are greater
 *          than for the CPUSolver, but the parallel performance and scaling
 *          are much better.
 */
class ThreadPrivateSolver : public CPUSolver {

protected:

  /** An array for the FSR scalar fluxes for each thread */
  FP_PRECISION** _thread_flux;

  /** An array for the CMFD Mesh surface currents for each thread */
  FP_PRECISION* _thread_currents;

  void initializeFluxArrays();
  void initializeCmfd();

  void flattenFSRFluxes(FP_PRECISION value);
  void zeroSurfaceCurrents();
  void scalarFluxTally(segment* curr_segment, int azim_index,
                       FP_PRECISION* track_flux,
                       FP_PRECISION* fsr_flux, bool fwd);
  void reduceThreadScalarFluxes();
  void reduceThreadSurfaceCurrents();
  void transportSweep();

public:
  ThreadPrivateSolver(Geometry* geometry=NULL,
                      TrackGenerator* track_generator=NULL,
                      Cmfd* cmfd=NULL);
  virtual ~ThreadPrivateSolver();
};


#endif /* THREADPRIVATESOLVER_H_ */
