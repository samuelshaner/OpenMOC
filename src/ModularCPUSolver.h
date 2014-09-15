/**
 * @file ModularCPUSolver.h
 * @brief The ModularCPUSolver class.
 * @date August 21, 2014
 * @author Sam Shaner, MIT, Course 22 (shaner@mit.edu)
 */


#ifndef MODULARCPUSOLVER_H_
#define MODULARCPUSOLVER_H_

#ifdef __cplusplus
#define _USE_MATH_DEFINES
#include <math.h>
#include <omp.h>
#include <stdlib.h>
#include "CPUSolver.h"
#include "ModularTrackGenerator.h"
#endif

/**
 * @class ModularCPUSolver ModularCPUSolver.h "src/ModularCPUSolver.h"
 * @brief This a subclass of the Solver class for multi-core CPUs using
 *        OpenMP multi-threading.
 * @details This Solver subclass uses OpenMP's multi-threading directives,
 *          including the mutual exclusion locks. Although the algorithm is
 *          more memory efficient than its ThreadPrivateSolver subclass, its
 *          parallel performance scales very poorly. As a result, this class
 *          is not recommended for general use and is primarily intended to
 *          expose the limitations of OpenMP's mutual exclusion formulation.
 */
class ModularCPUSolver : public CPUSolver {

protected:

  ModularTrackGenerator* _modular_track_generator;
  std::vector< std::vector< std::vector<Track*> > > _modular_tracks;
 
  bool _half_update;

  virtual void initializeFluxArrays();
  virtual void initializeFSRs();
  virtual void transportSweep();
  virtual void checkTrackSpacing();
  void transferBoundaryFluxModular(int track_id, int azim_index,
                                    bool direction,
                                   FP_PRECISION* track_flux, Track* track);
  void tallyCurrent(int track_id, int azim_index,
                    bool direction,
                    FP_PRECISION* track_flux, Track* track);
  void transferTrackFluxes();
  virtual FP_PRECISION computeFSRSources();
  
public:
  ModularCPUSolver(Geometry* geometry=NULL, TrackGenerator* track_generator=NULL);
  virtual ~ModularCPUSolver();

  void setHalfUpdate(bool half_update);

};


#endif /* MODULARCPUSOLVER_H_ */
