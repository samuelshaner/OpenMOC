/**
 * @file TrackGenerator.h
 * @brief The TrackGenerator class.
 * @date January 23, 2012
 * @author William Boyd, MIT, Course 22 (wboyd@mit.edu)
 */


#ifndef TRACKGENERATOR_H_
#define TRACKGENERATOR_H_

#ifdef __cplusplus
#define _USE_MATH_DEFINES
#include <iostream>
#include <fstream>
#include <sstream>
#include <unistd.h>
#include <omp.h>
#include "Track.h"
#include "Geometry.h"
#endif


/**
 * @class TrackGenerator TrackGenerator.h "src/TrackGenerator.h"
 * @brief The TrackGenerator is dedicated to generating and storing Tracks
 *        which cyclically wrap across the Geometry.
 * @details The TrackGenerator creates Track and initializes boundary
 *          conditions (vacuum or reflective) for each Track.
 */
class TrackGenerator {

protected:

  /** Number of azimuthal angles in \f$ [0, \pi] \f$ */
  int _num_azim;

  /** The track spacing (cm) */
  double _spacing;

  /** An integer array of the number of Tracks for each azimuthal angle */
  int* _num_tracks;

  /** The total number of Tracks for all azimuthal angles */
  int _tot_num_tracks;

  /** An integer array of the number of segments per Track  */
  int* _num_segments;

  /** The total number of segments for all Tracks */
  int _tot_num_segments;

  /** An integer array of the number of Tracks starting on the x-axis for each
   *  azimuthal angle */
  int* _num_x;

  /** An integer array of the number of Tracks starting on the y-axis for each
   *  azimuthal angle */
  int* _num_y;

  /** An array of the azimuthal angle quadrature weights */
  FP_PRECISION* _azim_weights;

  /** A 2D ragged array of Tracks */
  Track** _tracks;

  /** Pointer to the Geometry */
  Geometry* _geometry;

  /** Boolean for whether to use Track input file (true) or not (false) */
  bool _use_input_file;

  /** Filename for the *.tracks input / output file */
  std::string _tracks_filename;

  /** Boolean whether the Tracks have been generated (true) or not (false) */
  bool _contains_tracks;

  void computeEndPoint(Point* start, Point* end,  const double phi,
                       const double width, const double height);

  virtual void initializeTrackFileDirectory();
  void initializeTracks();
  void recalibrateTracksToOrigin();
  virtual void initializeBoundaryConditions();
  virtual void segmentize();
  void dumpTracksToFile();
  bool readTracksFromFile();

public:
  TrackGenerator(Geometry* geometry, int num_azim, double spacing);
  virtual ~TrackGenerator();

  /* Get parameters */
  int getNumAzim();
  double getTrackSpacing();
  Geometry* getGeometry();
  int getNumTracks();
  int* getNumTracksArray();
  int getNumSegments();
  int* getNumSegmentsArray();
  Track** getTracks();
  FP_PRECISION* getAzimWeights();

  /* Set parameters */
  void setNumAzim(int num_azim);
  void setTrackSpacing(double spacing);
  void setGeometry(Geometry* geometry);

  /* Worker functions */
  bool containsTracks();
  virtual void retrieveTrackCoords(double* coords, int num_tracks);
  virtual void retrieveSegmentCoords(double* coords, int num_segments);
  virtual void generateTracks();
};

#endif /* TRACKGENERATOR_H_ */
