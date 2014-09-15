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
 * @struct macro_track
 * @brief A macro_track represents the first and last track of 
 *        subtrack train.
 */
struct macro_track {

  /** A pointer to the start and end track */
  Track* _start;
  Track* _end;

};


/** Indexing macro for the scalar flux in each FSR and energy group */
#define curr_track(i,j) (_modular_track_map[&(_tracks[(i)][(j)])])

/** Indexing macro for the scalar flux in each FSR and energy group */
#define refl_track(i,j) (_modular_track_map[&(_tracks[_num_azim - (i) - 1][(j)])])


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

  int _cx;
  int _cy;
  Lattice* _lattice;
  
  std::vector< std::vector< std::vector<Track*> > > _modular_tracks;
  
  std::map<Track*, macro_track> _modular_track_map;

  /** The number of shared memory OpenMP threads */
  int _num_threads;


public:
  ModularTrackGenerator(Geometry* geometry, int num_azim, double spacing);
  virtual ~ModularTrackGenerator();

  int getNumX();
  int getNumY();
  void setNumX(int cx);
  void setNumY(int cy);
  void setLatticeStructure(int cx, int cy);
  int getNumCells();
  virtual void generateTracks();
  void decomposeTracks();
  void setEndPoint(Track* track, double dist);
  int findNextCell(Point* point, double angle);
  double nextCellDist(Point* point, double angle);
  virtual void initializeBoundaryConditions();
  virtual void segmentize();
  virtual void dumpTracksToFile();
  virtual void retrieveSegmentCoords(double* coords, int num_segments);
  virtual void initializeTrackFileDirectory();
  virtual void retrieveTrackCoords(double* coords, int num_tracks);
  int findDomainCell(LocalCoords* coords);
  void setNumThreads(int num_threads);

  std::vector< std::vector< std::vector<Track*> > > getModularTracks();
  std::map<Track*, macro_track> getModularTrackMap();
};

#endif /* MODULARTRACKGENERATOR_H_ */
