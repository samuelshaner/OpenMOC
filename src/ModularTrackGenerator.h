/**
 * @file ModularTrackGenerator.h
 * @brief The ModularTrackGenerator class.
 * @date August 20, 2014
 * @author Samuel Shaner, MIT, Course 22 (shaner@mit.edu)
 */


#ifndef MODULARTRACKGENERATOR_H_
#define MODULARTRACKGENERATOR_H_

#ifdef __cplusplus
#define _USE_MATH_DEFINES
#include <iostream>
#include <fstream>
#include <sstream>
#include <unistd.h>
#include <omp.h>
#include "Track.h"
#include "Geometry.h"
#include "TrackGenerator.h"
#endif


/**
 * @class ModularTrackGenerator ModularTrackGenerator.h 
 *        "src/ModularTrackGenerator.h"
 * @brief The ModularTrackGenerator is a subclass of TrackGenerator and
 *        generates tracks in a two-step procedure.
 * @details The ModularTrackGenerator creates Tracks and initializes boundary
 *          conditions (vacuum or reflective) for each Track. Tracks are first
 *          generated for the entire geometry and then broken up into subtracks
 *          by chopping them on a coarse cartesian mesh. The tracks are then
 *          segmented independently for each cartesian mesh cell/domain.
 */
class ModularTrackGenerator : public TrackGenerator {

private:

  /** Number of azimuthal angles in \f$ [0, \pi] \f$ */
  int _num_cells;
  int _cx;
  int _cy;

  vector<vector<vector<Track*> > > _modular_tracks;

  void computeEndPoint(Point* start, Point* end,  const double phi,
                       const double width, const double height);

public:
  ModularTrackGenerator(Geometry* geometry, int num_azim, double spacing);
  virtual ~ModularTrackGenerator();

};

#endif /* MODULARTRACKGENERATOR_H_ */
