#include "ModularTrackGenerator.h"


/**
 * @brief Constructor for the TrackGenerator assigns default values.
 * @param geometry a pointer to a Geometry object
 * @param num_azim number of azimuthal angles in \f$ [0, 2\pi] \f$
 * @param spacing track spacing (cm)
 */
ModularTrackGenerator::ModularTrackGenerator(Geometry* geometry, 
                                             const int num_azim,
                                             const double spacing) 
    :TrackGenerator(geometry, num_azim, spacing) {

  _cx = 1;
  _cy = 1;
}


/**
 * @brief Destructor frees memory for all Tracks.
 */
ModularTrackGenerator::~ModularTrackGenerator() {

}


int ModularTrackGenerator::getNumX() {
  return _cx;
}


void ModularTrackGenerator::setNumX(int cx) {
  _cx = cx;
}


int ModularTrackGenerator::getNumY() {
  return _cy;
}


void ModularTrackGenerator::setNumY(int cy) {
  _cy = cy;
}


int ModularTrackGenerator::getNumCells() {
  return _cx*_cy;
}


/**
 * @brief The structure of the Lattice to be used as the CMFD mesh.
 * @param The number of cells in the x direction.
 * @param The number of cells in the y direction.
 */
void ModularTrackGenerator::setLatticeStructure(int cx, int cy){
  setNumX(cx);
  setNumY(cy);
}



std::vector< std::vector< std::vector<Track*> > > ModularTrackGenerator::getModularTracks(){
  return _modular_tracks;
}


/**
 * @brief Generates tracks for some number of azimuthal angles and track spacing
 * @details Computes the effective angles and track spacing. Computes the
 *          number of Tracks for each azimuthal angle, allocates memory for
 *          all Tracks at each angle and sets each Track's starting and ending
 *          Points, azimuthal angle, and azimuthal angle quadrature weight.
 */
void ModularTrackGenerator::generateTracks() {

  if (_geometry == NULL)
    log_printf(ERROR, "Unable to generate Tracks since no Geometry "
               "has been set for the TrackGenerator");

  /* Deletes Tracks arrays if Tracks have been generated */
  if (_contains_tracks) {
    delete [] _num_tracks;
    delete [] _num_segments;
    delete [] _num_x;
    delete [] _num_y;
    delete [] _azim_weights;

    for (int i = 0; i < _num_azim; i++)
      delete [] _tracks[i];

    delete [] _tracks;
  }

  initializeTrackFileDirectory();

  /* If not Tracks input file exists, generate Tracks */
  if (_use_input_file == false) {

    /* Allocate memory for the Tracks */
    try {
      _num_tracks = new int[_num_azim];
      _num_x = new int[_num_azim];
      _num_y = new int[_num_azim];
      _azim_weights = new FP_PRECISION[_num_azim];
      _tracks = new Track*[_num_azim];
    }
    catch (std::exception &e) {
      log_printf(ERROR, "Unable to allocate memory for TrackGenerator. "
                 "Backtrace:\n%s", e.what());
    }

    /* Check to make sure that height, width of the Geometry are nonzero */
    if (_geometry->getHeight() <= 0 || _geometry->getHeight() <= 0)
      log_printf(ERROR, "The total height and width of the Geometry must be"
                 "must be nonzero for Track generation. Create a CellFill which"
                 "is filled by the entire geometry and bounded by XPlanes"
                 "and YPlanes to enable the Geometry to determine the total"
                 "width and height of the model.");

    /* Generate Tracks, perform ray tracing across the geometry, and store
     * the data to a Track file */
    try {
      initializeTracks();
      recalibrateTracksToOrigin();
      decomposeTracks();
      segmentize();
      //dumpTracksToFile();
    }
    catch (std::exception &e) {
      log_printf(ERROR, "Unable to allocate memory needed to generate "
                 "Tracks. Backtrace:\n%s", e.what());
    }
  }

  initializeBoundaryConditions();
  return;
}


/**
 * @brief This method creates a directory to store Track files, and reads
 *        in ray tracing data for Tracks and segments from a Track file
 *        if one exists.
 * @details This method is called by the TrackGenerator::generateTracks()
 *          class method. If a Track file exists for this Geometry, number
 *          of azimuthal angles, and track spacing, then this method will
 *          import the ray tracing Track and segment data to fill the
 *          appropriate data structures.
 */
void ModularTrackGenerator::initializeTrackFileDirectory() {

  std::stringstream directory;
  struct stat buffer;
  std::stringstream test_filename;

  /** Create directory to store Track files with pre-generated ray tracing data
   *  if the directory does not yet exist */

  directory << get_output_directory() << "/tracks";
  struct stat st;
  if (!stat(directory.str().c_str(), &st) == 0)
    mkdir(directory.str().c_str(), S_IRWXU);

  if (_geometry->getCmfd() != NULL){
    test_filename << directory.str() << "/"
                  <<  _num_azim*2.0 << "_angles_"
                  << _spacing << "_cm_spacing_cmfd_"
                  << _geometry->getCmfd()->getNumX() 
                  << "x" << _geometry->getCmfd()->getNumY()
                  << "_decomp_" << _cx << "x" << _cy
                  << ".data";
    }
  else{
    test_filename << directory.str() << "/"
                  <<  _num_azim*2.0 << "_angles_"
                  << _spacing << "_cm_spacing_"
                  << "_decomp_" << _cx << "x" << _cy << ".data";
  }

  _tracks_filename = test_filename.str();

  /* Check to see if a Track file exists for this geometry, number of azimuthal
   * angles, and track spacing, and if so, import the ray tracing data */
  if (!stat(_tracks_filename.c_str(), &buffer)) {
    if (readTracksFromFile()) {
      _use_input_file = true;
      _contains_tracks = true;
    }
  }
}


/**
 * @brief Initializes Track azimuthal angles, start and end Points.
 * @details This method computes the azimuthal angles and effective track
 *          spacing to use to guarantee cyclic Track wrapping. Based on the
 *          angles and spacing, the number of Tracks per angle and the start
 *          and end Points for each Track are computed.
 */
void ModularTrackGenerator::decomposeTracks() {

  log_printf(INFO, "Decomposing tracks across coarse mesh...");

  Track* old_track;
  Track* new_track;
  double height = _geometry->getHeight();
  double width = _geometry->getWidth();
  double cell_width = width / _cx;
  double cell_height = height / _cy;
  int cell;
  double phi;
  double dist;

  /* initialize 3D vector of tracks */
  for (int cell = 0; cell < _cx*_cy; cell++){

    std::vector< std::vector<Track*> > *tracks1 = new std::vector< std::vector<Track*> >;
    _modular_tracks.push_back(*tracks1);
    
    for (int i = 0; i < _num_azim; i++) {
      std::vector<Track*> *tracks2 = new std::vector<Track*>;
      _modular_tracks.at(cell).push_back(*tracks2);          
    }
  }

  /* inialize domain lattice */
  _lattice = new Lattice(0, cell_width, cell_height);
  _lattice->setNumX(_cx);
  _lattice->setNumY(_cy);

  log_printf(INFO, "decomposing tracks");

  /* Loop over tracks generated for the geometry and decompose them into 
   * subtracks for each previous track. Add them to _modular_tracks */
  for (int i = 0; i < _num_azim; i++) {

    phi = _tracks[i][0].getPhi();

    /* Loop the tracks for this azimuthal angle */
    for (int j = 0; j < _num_tracks[i]; j++) {

      /* add first track to first cell */
      Track* first_track = new Track();
      first_track->setPhi(phi);
      first_track->getStart()->setX(_tracks[i][j].getStart()->getX());
      first_track->getStart()->setY(_tracks[i][j].getStart()->getY());

      /* find cell that contains track */
      cell = findNextCell(first_track->getStart(), phi);

      /* get distance to next cell */
      dist = nextCellDist(first_track->getStart(), phi);

      /* set track end point */
      setEndPoint(first_track, dist);

      /* add track to modular tracks vector */
      _modular_tracks.at(cell).at(i).push_back(first_track);

      /* indicate the tracks starting point is on a geometry boundary */
      first_track->setOnBoundaryIn(1);
      
      /* get next cell */
      cell = findNextCell(first_track->getEnd(), phi);
      Track* new_track = first_track;
      
      while(cell != -1){

        /* indicate that the end point of the previous track is not on 
         * a geometry boundary */
        new_track->setReflOut(0);
        new_track->setBCOut(1);
        new_track->setOnBoundaryOut(0);

        /* initialize new track */
        Track* old_track = new_track;
        new_track = new Track();
        new_track->setPhi(phi);
        new_track->getStart()->setX(old_track->getEnd()->getX());
        new_track->getStart()->setY(old_track->getEnd()->getY());

        /* connect old and new tracks */
        old_track->setTrackOut(new_track);
        new_track->setTrackIn(old_track);

        /* indicate that the start point of the new track is not on 
         * a geometry boundary */
        new_track->setReflIn(1);
        new_track->setBCIn(1);
        new_track->setOnBoundaryIn(0);
        
        /* get distance to next cell */
        dist = nextCellDist(new_track->getStart(), phi);
        
        /* set track end point */
        setEndPoint(new_track, dist);

        /* add track to modular tracks vector */
        _modular_tracks.at(cell).at(i).push_back(new_track);
        
        /* increment tracks counter */
        _tot_num_tracks++;
        
        /* get next cell */
        cell = findNextCell(new_track->getEnd(), phi);
      }

      /* indicate the the end point of the last track is on 
       * a geometry boundary */
      new_track->setOnBoundaryOut(1);

      /* add macro track linking first and last subtracts to original macro track */
      macro_track* new_macro_track = new macro_track;
      new_macro_track->_start = first_track;
      new_macro_track->_end = new_track;
      _modular_track_map[&(_tracks[i][j])] = *new_macro_track;
    }
  }

  int uid = 0;
  std::vector<Track*>::iterator iter;

  /* Set the uid for each track. This is used by the Solver to map the 
   * track angular flux to a 1D angular flux array */
  for (int cell = 0; cell < _cx*_cy; cell++){
    for (int i=0; i < _num_azim; i++){
      for (iter = _modular_tracks.at(cell).at(i).begin(); iter != _modular_tracks.at(cell).at(i).end(); ++iter){
        (*iter)->setUid(uid);
        uid++;
      }
    }
  }
}


void ModularTrackGenerator::setEndPoint(Track* track, double dist) {
  double angle = track->getPhi();
  double new_x = track->getStart()->getX() + cos(angle) * dist;
  double new_y = track->getStart()->getY() + sin(angle) * dist;
  track->getEnd()->setCoords(new_x, new_y);
}


int ModularTrackGenerator::findNextCell(Point* point, double angle) {
  double delta_x = cos(angle) * TINY_MOVE;
  double delta_y = sin(angle) * TINY_MOVE;
  point->setCoords(point->getX() + delta_x, point->getY() + delta_y);
  int cell = _lattice->getLatticeCell(point);
  point->setCoords(point->getX() - delta_x, point->getY() - delta_y);
  return cell;
}


double ModularTrackGenerator::nextCellDist(Point* point, double angle) {
  double delta_x = cos(angle) * TINY_MOVE;
  double delta_y = sin(angle) * TINY_MOVE;
  point->setCoords(point->getX() + delta_x, point->getY() + delta_y);
  double dist = _lattice->minSurfaceDist(point, angle);
  dist += pow(pow(delta_x, 2) + pow(delta_y, 2), 0.5);
  point->setCoords(point->getX() - delta_x, point->getY() - delta_y);
  return dist;
}


/**
 * @brief Initializes boundary conditions for each Track.
 * @details Sets boundary conditions by setting the incoming and outgoing Tracks
 *          for each Track using a special indexing scheme into the 2D jagged
 *          array of Tracks.
 */
void ModularTrackGenerator::initializeBoundaryConditions() {

  log_printf(INFO, "Initializing Track boundary conditions...");

  /* nxi = number of tracks starting on y-axis for angle i
   * nyi = number of tracks starting on y-axis for angle i
   * nti = total number of tracks for angle i */
  int nxi, nyi, nti;

  Track *curr;
  Track *refl;

  /* Loop over only half the angles since we will set the pointers for
   * connecting Tracks at the same time */
  for (int i = 0; i < floor(_num_azim / 2); i++) {
    nxi = _num_x[i];
    nyi = _num_y[i];
    nti = _num_tracks[i];
    curr = _tracks[i];
    refl = _tracks[_num_azim - i - 1];

    /* Loop over all of the Tracks for this angle */
    for (int j = 0; j < nti; j++) {

      /* More Tracks starting along x-axis than y-axis */
      if (nxi <= nyi) {

        /* Bottom to right hand side */
        if (j < nxi) {
          curr_track(i,j)._start->setTrackIn(refl_track(i,j)._start);
          refl_track(i,j)._start->setTrackIn(curr_track(i,j)._start);

          curr_track(i,j)._start->setReflIn(false);
          refl_track(i,j)._start->setReflIn(false);

          if (_geometry->getBCBottom() == REFLECTIVE) {
            curr_track(i,j)._start->setBCIn(1);
            refl_track(i,j)._start->setBCIn(1);
          }
          else {
            curr_track(i,j)._start->setBCIn(0);
            refl_track(i,j)._start->setBCIn(0);
          }

          curr_track(i,j)._end->setTrackOut(refl_track(i,2 * nxi - 1 - j)._start);
          refl_track(i,2 * nxi - 1 - j)._start->setTrackIn(curr_track(i,j)._end);

          curr_track(i,j)._end->setReflOut(false);
          refl_track(i,2 * nxi - 1 - j)._start->setReflIn(true);

          if (_geometry->getBCRight() == REFLECTIVE) {
            curr_track(i,j)._end->setBCOut(1);
            refl_track(i,2 * nxi - 1 - j)._start->setBCIn(1);
          }
          else {
            curr_track(i,j)._end->setBCOut(0);
            refl_track(i,2 * nxi - 1 - j)._start->setBCIn(0);
          }
        }

        /* Left hand side to right hand side */
        else if (j < nyi) {
          curr_track(i,j)._start->setTrackIn(refl_track(i,j - nxi)._end);
          refl_track(i,j - nxi)._end->setTrackOut(curr_track(i,j)._start);

          curr_track(i,j)._start->setReflIn(true);
          refl_track(i,j - nxi)._end->setReflOut(false);

          if (_geometry->getBCLeft() == REFLECTIVE) {
            curr_track(i,j)._start->setBCIn(1);
            refl_track(i,j - nxi)._end->setBCOut(1);
          }
          else {
            curr_track(i,j)._start->setBCIn(0);
            refl_track(i,j - nxi)._end->setBCOut(0);
          }

          curr_track(i,j)._end->setTrackOut(refl_track(i,j + nxi)._start);
          refl_track(i,j + nxi)._start->setTrackIn(curr_track(i,j)._end);

          curr_track(i,j)._end->setReflOut(false);
          refl_track(i,j + nxi)._start->setReflIn(true);

          if (_geometry->getBCRight() == REFLECTIVE) {
            curr_track(i,j)._end->setBCOut(1);
            refl_track(i,j + nxi)._start->setBCIn(1);
          }
          else {
            curr_track(i,j)._end->setBCOut(0);
            refl_track(i,j + nxi)._start->setBCIn(0);
          }
        }

        /* Left hand side to top (j > ny) */
        else {
          curr_track(i,j)._start->setTrackIn(refl_track(i,j - nxi)._end);
          refl_track(i,j - nxi)._end->setTrackOut(curr_track(i,j)._start);

          curr_track(i,j)._start->setReflIn(true);
          refl_track(i,j - nxi)._end->setReflOut(false);

          if (_geometry->getBCLeft() == REFLECTIVE) {
            curr_track(i,j)._start->setBCIn(1);
            refl_track(i,j - nxi)._end->setBCOut(1);
          }
          else {
            curr_track(i,j)._start->setBCIn(0);
            refl_track(i,j - nxi)._end->setBCOut(0);
          }

          curr_track(i,j)._end->setTrackOut(refl_track(i,2 * nti - nxi - j - 1)._end);
          refl_track(i,2 * nti - nxi - j - 1)._end->setTrackOut(curr_track(i,j)._end);

          curr_track(i,j)._end->setReflOut(true);
          refl_track(i,2 * nti - nxi - j - 1)._end->setReflOut(true);

          if (_geometry->getBCTop() == REFLECTIVE) {
            curr_track(i,j)._end->setBCOut(1);
            refl_track(i,2 * nti - nxi - j - 1)._end->setBCOut(1);
          }
          else {
            curr_track(i,j)._end->setBCOut(0);
            refl_track(i,2 * nti - nxi - j - 1)._end->setBCOut(0);
          }
        }
      }

      /* More Tracks starting on y-axis than on x-axis */
      else {

        /* Bottom to top */
        if (j < nxi - nyi) {
          curr_track(i,j)._start->setTrackIn(refl_track(i,j)._start);
          refl_track(i,j)._start->setTrackIn(curr_track(i,j)._start);

          curr_track(i,j)._start->setReflIn(false);
          refl_track(i,j)._start->setReflIn(false);

          if (_geometry->getBCBottom() == REFLECTIVE) {
            curr_track(i,j)._start->setBCIn(1);
            refl_track(i,j)._start->setBCIn(1);
         }
          else {
            curr_track(i,j)._start->setBCIn(0);
            refl_track(i,j)._start->setBCIn(0);
          }

          curr_track(i,j)._end->setTrackOut(refl_track(i,nti - (nxi - nyi) + j)._end);
          refl_track(i,nti - (nxi - nyi) + j)._end->setTrackOut(curr_track(i,j)._end);

          curr_track(i,j)._end->setReflOut(true);
          refl_track(i,nti - (nxi - nyi) + j)._end->setReflOut(true);

          if (_geometry->getBCTop() == REFLECTIVE) {
            curr_track(i,j)._end->setBCOut(1);
            refl_track(i,nti - (nxi - nyi) + j)._end->setBCOut(1);
          }
          else {
            curr_track(i,j)._end->setBCOut(0);
            refl_track(i,nti - (nxi - nyi) + j)._end->setBCOut(0);
          }
        }

        /* Bottom to right hand side */
        else if (j < nxi) {
          curr_track(i,j)._start->setTrackIn(refl_track(i,j)._start);
          refl_track(i,j)._start->setTrackIn(curr_track(i,j)._start);

          curr_track(i,j)._start->setReflIn(false);
          refl_track(i,j)._start->setReflIn(false);

          if (_geometry->getBCBottom() == REFLECTIVE) {
            curr_track(i,j)._start->setBCIn(1);
            refl_track(i,j)._start->setBCIn(1);
          }
          else {
            curr_track(i,j)._start->setBCIn(0);
            refl_track(i,j)._start->setBCIn(0);
          }

          curr_track(i,j)._end->setTrackOut(refl_track(i,nxi + (nxi - j) - 1)._start);
          refl_track(i,nxi + (nxi - j) - 1)._start->setTrackIn(curr_track(i,j)._end);

          curr_track(i,j)._end->setReflOut(false);
          refl_track(i,nxi + (nxi - j) - 1)._start->setReflIn(true);

          if (_geometry->getBCRight() == REFLECTIVE) {
            curr_track(i,j)._end->setBCOut(1);
            refl_track(i,nxi + (nxi - j) - 1)._start->setBCIn(1);
          }
          else {
            curr_track(i,j)._end->setBCOut(0);
            refl_track(i,nxi + (nxi - j) - 1)._start->setBCIn(0);
          }
        }

        /* Left-hand side to top (j > nx) */
        else {
          curr_track(i,j)._start->setTrackIn(refl_track(i,j - nxi)._end);
          refl_track(i,j - nxi)._end->setTrackOut(curr_track(i,j)._start);

          curr_track(i,j)._start->setReflIn(true);
          refl_track(i,j - nxi)._end->setReflOut(false);

          if (_geometry->getBCLeft() == REFLECTIVE) {
            curr_track(i,j)._start->setBCIn(1);
            refl_track(i,j - nxi)._end->setBCOut(1);
          }
          else {
            curr_track(i,j)._start->setBCIn(0);
            refl_track(i,j - nxi)._end->setBCOut(0);
          }

          curr_track(i,j)._end->setTrackOut(refl_track(i,nyi + (nti - j) - 1)._end);
          refl_track(i,nyi + (nti - j) - 1)._end->setTrackOut(curr_track(i,j)._end);

          curr_track(i,j)._end->setReflOut(true);
          refl_track(i,nyi + (nti - j) - 1)._end->setReflOut(true);

          if (_geometry->getBCTop() == REFLECTIVE) {
            curr_track(i,j)._end->setBCOut(1);
            refl_track(i,nyi + (nti - j) - 1)._end->setBCOut(1);
          }
          else {
            curr_track(i,j)._end->setBCOut(0);
            refl_track(i,nyi + (nti - j) - 1)._end->setBCOut(0);
          }
        }
      }
    }
  }

  return;
}


/**
 * @brief Generate segments for each Track across the Geometry.
 */
void ModularTrackGenerator::segmentize() {

  log_printf(NORMAL, "Ray tracing for modular track segmentation...");

  Track* track;
  std::vector<Track*>::iterator iter;

  if (_num_segments != NULL)
    delete [] _num_segments;

  /* This section loops over all Track and segmentizes each one if the
   * Tracks were not read in from an input file */
  if (!_use_input_file) {

    log_printf(NORMAL, "segmenting tracks...");

    /* Loop over all Tracks */
    #pragma omp parallel for private(iter, track)
    for (int cell = 0; cell < _cx*_cy; cell++){
      for (int i=0; i < _num_azim; i++) {
        for (iter = _modular_tracks.at(cell).at(i).begin(); 
             iter != _modular_tracks.at(cell).at(i).end(); ++iter){
          track = *iter;
          log_printf(DEBUG, "Segmenting Track %d/%d with i = %d, cell = %d",
                     track->getUid(), _modular_tracks.at(cell).at(i).size(),
                     i, cell);
          _geometry->segmentize(track);
        }
      }
    }

    log_printf(NORMAL, "done segmenting tracks...");

    /* Compute the total number of segments in the simulation */
    _num_segments = new int[_tot_num_tracks];
    _tot_num_segments = 0;

    for (int cell = 0; cell < _cx*_cy; cell++){
      for (int i=0; i < _num_azim; i++) {
        for (iter = _modular_tracks.at(cell).at(i).begin(); 
             iter != _modular_tracks.at(cell).at(i).end(); ++iter){
          track = *iter;
          _num_segments[track->getUid()] = track->getNumSegments();
          _tot_num_segments += _num_segments[track->getUid()];
        }
      }
    }
  }
  
  _contains_tracks = true;
  
  return;
}


/**
 * @brief Fills an array with the x,y coordinates for each Track segment.
 * @details This class method is intended to be called by the OpenMOC
 *          Python "plotter" module as a utility to assist in plotting
 *          segments. Although this method appears to require two arguments,
 *          in reality it only requires one due to SWIG and would be called
 *          from within Python as follows:
 *
 * @code
 *          num_segments = track_generator.getNumSegments()
 *          coords = track_generator.retrieveSegmentCoords(num_segments*5)
 * @endcode
 *
 * @param coords an array of coords of length 5 times the number of segments
 * @param num_segments the total number of Track segments
 */
void ModularTrackGenerator::retrieveSegmentCoords(double* coords, int num_segments) {

  if (num_segments != 5*getNumSegments())
    log_printf(ERROR, "Unable to retrieve the Track segment coordinates since "
               "the TrackGenerator contains %d segments with %d coordinates "
               "but an array of length %d was input",
               getNumSegments(), 5*getNumSegments(), num_segments);

  segment* curr_segment = NULL;
  double x0, x1, y0, y1;
  double phi;
  segment* segments;
  Track* track;

  int counter = 0;
  std::vector<Track*>::iterator iter;

  /* Loop over Track segments and populate array with their FSR ID and *
   * start/end points */
  for (int cell=0; cell < _cx*_cy; cell++) {
    for (int i=0; i < _num_azim; i++) {
      for (iter = _modular_tracks.at(cell).at(i).begin(); 
          iter != _modular_tracks.at(cell).at(i).end(); ++iter){
      
        track = *iter;
        x0 = track->getStart()->getX();
        y0 = track->getStart()->getY();
        phi = track->getPhi();
        
        segments = track->getSegments();

        for (int s=0; s < track->getNumSegments(); s++) {
          curr_segment = &segments[s];

          coords[counter] = curr_segment->_region_id;
          //coords[counter] = cell;

          coords[counter+1] = x0;
          coords[counter+2] = y0;
          
          x1 = x0 + cos(phi) * curr_segment->_length;
          y1 = y0 + sin(phi) * curr_segment->_length;
          
          coords[counter+3] = x1;
          coords[counter+4] = y1;
          
          x0 = x1;
          y0 = y1;
          
          counter += 5;
        }
      }
    }
  }

  return;
}


/**
 * @brief Fills an array with the x,y coordinates for each Track.
 * @details This class method is intended to be called by the OpenMOC
 *          Python "plotter" module as a utility to assist in plotting
 *          tracks. Although this method appears to require two arguments,
 *          in reality it only requires on due to SWIG and would be called
 *          from within Python as follows:
 *
 * @code
 *          num_tracks = track_generator.getNumTracks()
 *          coords = track_generator.retrieveTrackCoords(num_tracks*4)
 * @endcode
 *
 * @param coords an array of coords of length 4 times the number of Tracks
 * @param num_tracks the total number of Tracks
 */
void ModularTrackGenerator::retrieveTrackCoords(double* coords, int num_tracks) {

  if (num_tracks != 4*getNumTracks())
    log_printf(ERROR, "Unable to retrieve the Track coordinates since the "
               "TrackGenerator contains %d Tracks with %d coordinates but an "
               "array of length %d was input",
               getNumTracks(), 4*getNumTracks(), num_tracks);

  int counter = 0;
  std::vector<Track*>::iterator iter;
  Track* track;

  /* Loop over Track segments and populate array with their FSR ID and *
   * start/end points */
  for (int cell=0; cell < _cx*_cy; cell++) {
    for (int i=0; i < _num_azim; i++) {
      for (iter = _modular_tracks.at(cell).at(i).begin(); 
          iter != _modular_tracks.at(cell).at(i).end(); ++iter){
      
        track = *iter;
        coords[counter] = track->getStart()->getX();
        coords[counter+1] = track->getStart()->getY();
        coords[counter+2] = track->getEnd()->getX();
        coords[counter+3] = track->getEnd()->getY();
        counter += 4;
      }
    }
  }

  return;
}


