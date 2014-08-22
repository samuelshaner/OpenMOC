#include "ModularCPUSolver.h"


/**
 * @brief Constructor initializes array pointers for Tracks and Materials.
 * @details The constructor retrieves the number of energy groups and FSRs
 *          and azimuthal angles from the Geometry and TrackGenerator if
 *          passed in as parameters by the user. The constructor initalizes
 *          the number of OpenMP threads to a default of 1.
 * @param geometry an optional pointer to the Geometry
 * @param track_generator an optional pointer to the TrackGenerator
 * @param cmfd an optional pointer to a Cmfd object object
 */
ModularCPUSolver::ModularCPUSolver(Geometry* geometry, TrackGenerator* track_generator) 
    : CPUSolver(geometry, track_generator) {

}


/**
 * @brief Destructor deletes array for OpenMP mutual exclusion locks for
 *        FSR scalar flux updates, and calls Solver parent class destructor
 *        to deletes arrays for fluxes and sources.
 */
ModularCPUSolver::~ModularCPUSolver() {
}


/**
 * @brief Allocates memory for Track boundary angular flux and leakage
 *        and FSR scalar flux arrays.
 * @details Deletes memory for old flux arrays if they were allocated for a
 *          previous simulation.
 */
void ModularCPUSolver::initializeFluxArrays() {

  /* Delete old flux arrays if they exist */
  if (_boundary_flux != NULL)
    delete [] _boundary_flux;

  if (_boundary_leakage != NULL)
    delete [] _boundary_leakage;

  if (_scalar_flux != NULL)
    delete [] _scalar_flux;

  int size;

  /* Allocate memory for the Track boundary flux and leakage arrays */
  try{
    size = 2 * _tot_num_tracks * _polar_times_groups;
    _boundary_flux = new FP_PRECISION[size];
    _boundary_leakage = new FP_PRECISION[size];

    /* Allocate an array for the FSR scalar flux */
    size = _num_FSRs * _num_groups;
    _scalar_flux = new FP_PRECISION[size];
  }
  catch(std::exception &e) {
    log_printf(ERROR, "Could not allocate memory for the Solver's fluxes. "
               "Backtrace:%s", e.what());
  }
}


/**
 * @brief Initializes the FSR volumes and Materials array.
 * @details This method assigns each FSR a unique, monotonically increasing
 *          ID, sets the Material for each FSR, and assigns a volume based on
 *          the cumulative length of all of the segments inside the FSR.
 */
void ModularCPUSolver::initializeFSRs() {

  log_printf(INFO, "Initializing flat source regions...");

  /* Delete old FSR arrays if they exist */
  if (_FSR_volumes != NULL)
    delete [] _FSR_volumes;

  if (_FSR_materials != NULL)
    delete [] _FSR_materials;

  _FSR_volumes = (FP_PRECISION*)calloc(_num_FSRs, sizeof(FP_PRECISION));
  _FSR_materials = new Material*[_num_FSRs];
  _FSR_locks = new omp_lock_t[_num_FSRs];
  _modular_track_generator = static_cast<ModularTrackGenerator*>(_track_generator);
  _modular_tracks = _modular_track_generator->getModularTracks();

  int num_segments;
  segment* curr_segment;
  segment* segments;
  FP_PRECISION volume;
  CellBasic* cell;
  Material* material;
  Universe* univ_zero = _geometry->getUniverse(0);
  Track* track;

  std::vector<Track*>::iterator iter;

  /* Set each FSR's "volume" by accumulating the total length of all segments
   * inside the FSR. Loop over azimuthal angles, Tracks and Track segments. */
  for (int cell = 0; cell < _modular_track_generator->getNumCells(); cell++){
    for (int i=0; i < _num_azim; i++) {
      for (iter = _modular_tracks.at(cell).at(i).begin(); iter != _modular_tracks.at(cell).at(i).end(); ++iter){

        track = *iter;

        num_segments = track->getNumSegments();
        segments = track->getSegments();

        for (int s=0; s < num_segments; s++) {
          curr_segment = &segments[s];
          volume = curr_segment->_length * _azim_weights[i];
          _FSR_volumes[curr_segment->_region_id] += volume;
        }
      }
    }
  }

  /* Loop over all FSRs to extract FSR material pointers */
  #pragma omp parallel for private(cell, material) schedule(guided)
  for (int r=0; r < _num_FSRs; r++) {

    /* Get the Cell corresponding to this FSR from the geometry */
    cell = _geometry->findCellContainingFSR(r);

    /* Get the Cell's Material and assign it to the FSR */
    material = _geometry->getMaterial(cell->getMaterial());
    _FSR_materials[r] = material;

    log_printf(DEBUG, "FSR ID = %d has Cell ID = %d and Material ID = %d "
               "and volume = %f", r, cell->getId(),
               _FSR_materials[r]->getUid(), _FSR_volumes[r]);
  }

  /* Loop over all FSRs to initialize OpenMP locks */
  #pragma omp parallel for schedule(guided)
  for (int r=0; r < _num_FSRs; r++)
     omp_init_lock(&_FSR_locks[r]);

  return;
}


/**
 * @brief This method performs one transport sweep of all azimuthal angles,
 *        Tracks, Track segments, polar angles and energy groups.
 * @details The method integrates the flux along each Track and updates the
 *          boundary fluxes for the corresponding output Track, while updating
 *          the scalar flux in each flat source region.
 */
void ModularCPUSolver::transportSweep() {

  Track* curr_track;
  int track_id;
  int num_segments;
  segment* curr_segment;
  segment* segments;
  FP_PRECISION* track_flux;
  int cx = _modular_track_generator->getNumX();
  int cy = _modular_track_generator->getNumY();
  int cell;
  log_printf(DEBUG, "Transport sweep with %d OpenMP threads", _num_threads);

  /* Initialize flux in each FSr to zero */
  flattenFSRFluxes(0.0);

  if (_cmfd->getOverlayMesh() && _cmfd->getUpdateFlux())
    zeroSurfaceCurrents();

  /* loop over black domain cells */
  #pragma omp parallel for private(curr_track, track_id, num_segments, \
    curr_segment, segments, track_flux, cell) schedule(guided)
  for (int y = 0; y < cy; y++){
    for (int x = y % 2; x < cx; x += 2){
      cell = y*cx+x;

      /* Loop over azimuthal angles */
      for (int i=0; i < _num_azim; i++) {

        std::vector<Track*>::iterator iter;

        /* iterate over tracks in this cell and azimuthal angle */
        for (iter = _modular_tracks.at(cell).at(i).begin(); iter != _modular_tracks.at(cell).at(i).end(); ++iter){

          /* Use local array accumulator to prevent false sharing*/
          FP_PRECISION* thread_fsr_flux;
          thread_fsr_flux = new FP_PRECISION[_num_groups];

          /* Initialize local pointers to important data structures */
          curr_track = *iter;
          num_segments = curr_track->getNumSegments();
          segments = curr_track->getSegments();
          track_id = curr_track->getUid();
          track_flux = &_boundary_flux(track_id,0,0,0);

          /* Loop over each Track segment in forward direction */
          for (int s=0; s < num_segments; s++) {
            curr_segment = &segments[s];
            scalarFluxTally(curr_segment, i, track_flux,
                            thread_fsr_flux,true);
          }

          /* Transfer boundary angular flux to outgoing Track */
          transferBoundaryFluxModular(track_id, i, true, track_flux, curr_track);
          
          /* Loop over each Track segment in reverse direction */
          track_flux += _polar_times_groups;
          
          for (int s=num_segments-1; s > -1; s--) {
              curr_segment = &segments[s];
              scalarFluxTally(curr_segment, i, track_flux,
                              thread_fsr_flux,false);
          }
          delete thread_fsr_flux;

          /* Transfer boundary angular flux to outgoing Track */
          transferBoundaryFluxModular(track_id, i, false, track_flux, curr_track);

        }
      }
    }
  }

  /* loop over red domain cells */
  #pragma omp parallel for private(curr_track, track_id, num_segments, \
    curr_segment, segments, track_flux, cell) schedule(guided)
  for (int y = 0; y < cy; y++){
    for (int x = 1 - y % 2; x < cx; x += 2){
      cell = y*cx+x;
      
      /* Loop over azimuthal angles */
      for (int i=0; i < _num_azim; i++) {

        std::vector<Track*>::iterator iter;

        /* iterate over tracks in this cell and azimuthal angle */
        for (iter = _modular_tracks.at(cell).at(i).begin(); iter != _modular_tracks.at(cell).at(i).end(); ++iter){

          /* Use local array accumulator to prevent false sharing*/
          FP_PRECISION* thread_fsr_flux;
          thread_fsr_flux = new FP_PRECISION[_num_groups];

          /* Initialize local pointers to important data structures */
          curr_track = *iter;
          num_segments = curr_track->getNumSegments();
          segments = curr_track->getSegments();
          track_id = curr_track->getUid();
          track_flux = &_boundary_flux(track_id,0,0,0);

          /* Loop over each Track segment in forward direction */
          for (int s=0; s < num_segments; s++) {
            curr_segment = &segments[s];
            scalarFluxTally(curr_segment, i, track_flux,
                            thread_fsr_flux,true);
          }

          /* Transfer boundary angular flux to outgoing Track */
          transferBoundaryFluxModular(track_id, i, true, track_flux, curr_track);
          
          /* Loop over each Track segment in reverse direction */
          track_flux += _polar_times_groups;
          
          for (int s=num_segments-1; s > -1; s--) {
              curr_segment = &segments[s];
              scalarFluxTally(curr_segment, i, track_flux,
                              thread_fsr_flux,false);
          }
          delete thread_fsr_flux;
          
          /* Transfer boundary angular flux to outgoing Track */
          transferBoundaryFluxModular(track_id, i, false, track_flux, curr_track);
        }
      }
    }
  }

  return;
}


/**
 * @brief Updates the boundary flux for a Track given boundary conditions.
 * @details For reflective boundary conditions, the outgoing boundary flux
 *          for the Track is given to the reflecting Track. For vacuum
 *          boundary conditions, the outgoing flux tallied as leakage.
 * @param track_id the ID number for the Track of interest
 * @param azim_index a pointer to the azimuthal angle index for this segment
 * @param direction the Track direction (forward - true, reverse - false)
 * @param track_flux a pointer to the Track's outgoing angular flux
 */
void ModularCPUSolver::transferBoundaryFluxModular(int track_id,
                                                   int azim_index,
                                                   bool direction,
                                                   FP_PRECISION* track_flux,
                                                   Track* track) {
  int start;
  int bc;
  FP_PRECISION* track_leakage;
  int track_out_id;
  int on_boundary;

  /* Extract boundary conditions for this Track and the pointer to the
   * outgoing reflective Track, and index into the leakage array */

  /* For the "forward" direction */
  if (direction) {
    start = track->isReflOut() * _polar_times_groups;
    bc = (int)track->getBCOut();
    track_leakage = &_boundary_leakage(track_id,0);
    track_out_id = track->getTrackOut()->getUid();
    on_boundary = (int)track->getOnBoundaryOut();
  }

  /* For the "reverse" direction */
  else {
    start = track->isReflIn() * _polar_times_groups;
    bc = (int)track->getBCIn();
    track_leakage = &_boundary_leakage(track_id,_polar_times_groups);
    track_out_id = track->getTrackIn()->getUid();
    on_boundary = (int)track->getOnBoundaryIn();
  }

  FP_PRECISION* track_out_flux = &_boundary_flux(track_out_id,0,0,start);

  /* Loop over polar angles and energy groups */
  for (int e=0; e < _num_groups; e++) {
    for (int p=0; p < _num_polar; p++) {
      track_out_flux(p,e) = track_flux(p,e) * bc;
      track_leakage(p,e) = track_flux(p,e) * on_boundary *
                            _polar_weights(azim_index,p) * (!bc);
    }
  }
}


/**
 * @brief Checks that each FSR has at least one Track segment crossing it
 *        and if not, throws an exception and prints an error message.
 * @details This method is for internal use only and is called by the
 *          Solver::convergeSource() method and should not be called
 *          directly by the user.
 */
void ModularCPUSolver::checkTrackSpacing() {

  int* FSR_segment_tallies = new int[_num_FSRs];
  int num_segments;
  segment* curr_segment;
  segment* segments;
  Cell* cell;
  Track* track;

  /* Set each tally to zero to begin with */
  #pragma omp parallel for
  for (int r=0; r < _num_FSRs; r++)
    FSR_segment_tallies[r] = 0;

  std::vector<Track*>::iterator iter;

  for (int cell = 0; cell < _modular_track_generator->getNumCells(); cell++){
    for (int i=0; i < _num_azim; i++) {
      for (iter = _modular_tracks.at(cell).at(i).begin(); iter != _modular_tracks.at(cell).at(i).end(); ++iter){

        track = *iter;
        num_segments = track->getNumSegments();
        segments = track->getSegments();

        for (int s=0; s < num_segments; s++) {
          curr_segment = &segments[s];
          FSR_segment_tallies[curr_segment->_region_id]++;
        }
      }
    }
  }

  /* Loop over all FSRs and if one FSR does not have tracks in it, print
   * error message to the screen and exit program */
  #pragma omp parallel for private (cell)
  for (int r=0; r < _num_FSRs; r++) {

    if (FSR_segment_tallies[r] == 0) {
      log_printf(ERROR, "No tracks were tallied inside FSR id = %d. Please "
                 "reduce your track spacing, increase the number of azimuthal"
                 "angles, or increase the size of the FSRs", r);
    }
  }

  delete [] FSR_segment_tallies;
}
