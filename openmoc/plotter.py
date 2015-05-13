#
# @file plotter.py
# @package openmoc.plotter
# @brief The plotter module provides utility functions to plot data from
#        OpenMOCs C++ classes, in particular, the geomery, including Material,
#        Cells and flat source regions, and fluxes and pin powers.
# @author William Boyd (wboyd@mit.edu)
# @date March 10, 2013

import sys

## @var openmoc
#  @brief The openmoc module in use in the Python script using the
#         openmoc.plotter module.
openmoc = ''

print('whoa!!!')

# Determine which OpenMOC module is being used
if 'openmoc.gnu.double' in sys.modules:
  openmoc = sys.modules['openmoc.gnu.double']
elif 'openmoc.gnu.single' in sys.modules:
  openmoc = sys.modules['openmoc.gnu.single']
elif 'openmoc.intel.double' in sys.modules:
  openmoc = sys.modules['openmoc.intel.double']
elif 'openmoc.intel.single' in sys.modules:
  openmoc = sys.modules['openmoc.intel.single']
elif 'openmoc.bgq.double' in sys.modules:
  openmoc = sys.modules['openmoc.bgq.double']
elif 'openmoc.bgq.single' in sys.modules:
  openmoc = sys.modules['openmoc.bgq.single']
else:
  import openmoc


import matplotlib

# force headless backend, or set 'backend' to 'Agg'
# in your ~/.matplotlib/matplotlibrc
matplotlib.use('Agg')

import matplotlib.pyplot as plt

# Force non-interactive mode, or set 'interactive' to False
# in your ~/.matplotlib/matplotlibrc
plt.ioff()

import matplotlib.colors as colors
import matplotlib.cm as cmx
import numpy as np
import numpy.random
import os, sys

# For Python 2.X.X
if (sys.version_info[0] == 2):
  from log import *
  from process import *
# For Python 3.X.X
else:
  from openmoc.log import *
  from openmoc.process import *


## A static variable for the output directory in which to save plots
subdirectory = "/plots/"

TINY_MOVE = openmoc.TINY_MOVE


##
# @brief Plots the characteristic tracks from an OpenMOC simulation.
# @details This method requires that Tracks have been generated by a
#          TrackGenerator object. A user may invoke this function from
#          an OpenMOC Python file as follows:
#
# @code
#         openmoc.plotter.plot_tracks(track_generator)
# @endcode
#
# @param track_generator the TrackGenerator which has generated Tracks
def plot_tracks(track_generator):

  global subdirectory

  directory = openmoc.get_output_directory() + subdirectory

  # Make directory if it does not exist
  if not os.path.exists(directory):
    os.makedirs(directory)

  # Error checking
  if not 'TrackGenerator' in str(type(track_generator)):
    py_printf('ERROR', 'Unable to plot Tracks since %s was input rather ' + \
              'than a TrackGenerator', str(type(track_generator)))

  if not track_generator.containsTracks():
    py_printf('ERROR', 'Unable to plot Tracks since the track ' + \
              'generator has not yet generated tracks')

  py_printf('NORMAL', 'Plotting the tracks...')

  # Retrieve data from TrackGenerator
  num_azim = track_generator.getNumAzim()
  spacing = track_generator.getTrackSpacing()
  num_tracks = track_generator.getNumTracks()
  coords = track_generator.retrieveTrackCoords(num_tracks*4)

  # Convert data to NumPy arrays
  coords = np.array(coords)
  x = coords[0::2]
  y = coords[1::2]

  # Make figure of line segments for each Track
  fig = plt.figure()
  for i in range(num_tracks):
    plt.plot([x[i*2], x[i*2+1]], [y[i*2], y[i*2+1]], 'b-')

  plt.xlim([x.min(), x.max()])
  plt.ylim([y.min(), y.max()])

  title = 'Tracks for ' + str(num_azim) + ' angles and ' + str(spacing) + \
            ' cm spacing'

  plt.title(title)

  filename = directory + 'tracks-' + str(num_azim) + '-angles-' + \
      str(spacing) + '-spacing.png'

  fig.savefig(filename, bbox_inches='tight')
  plt.close(fig)


##
# @brief Plots the characteristic Track segments from an OpenMOC simulation.
# @details This method requires that tracks have been generated by a
#          TrackGenerator object. Each segment is colored by the ID of the
#          unique flat flat source region it is within. A user may invoke
#          this function from an OpenMOC Python file as follows:
#
# @code
#         openmoc.plotter.plot_segments(track_generator)
# @endcode
#
# @param track_generator the TrackGenerator which has generated Tracks
def plot_segments(track_generator):

  global subdirectory

  directory = openmoc.get_output_directory() + subdirectory

  # Make directory if it does not exist
  if not os.path.exists(directory):
    os.makedirs(directory)

  # Error checking
  if not 'TrackGenerator' in str(type(track_generator)):
    py_printf('ERROR', 'Unable to plot Track segments since %s was input ' + \
              'rather than a TrackGenerator', str(type(track_generator)))

  if not track_generator.containsTracks():
    py_printf('ERROR', 'Unable to plot Track segments since the ' + \
              'TrackGenerator has not yet generated Tracks.')

  py_printf('NORMAL', 'Plotting the track segments...')

  # Retrieve data from TrackGenerator
  num_azim = track_generator.getNumAzim()
  spacing = track_generator.getTrackSpacing()
  num_segments = track_generator.getNumSegments()
  num_fsrs = track_generator.getGeometry().getNumFSRs()
  coords = track_generator.retrieveSegmentCoords(num_segments*5)

  # Convert data to NumPy arrays
  coords = np.array(coords)
  x = numpy.zeros(num_segments*2)
  y = numpy.zeros(num_segments*2)
  fsrs = numpy.zeros(num_segments)

  for i in range(num_segments):
    fsrs[i] = coords[i*5]
    x[i*2] = coords[i*5+1]
    y[i*2] = coords[i*5+2]
    x[i*2+1] = coords[i*5+3]
    y[i*2+1] = coords[i*5+4]

  # Create array of equally spaced randomized floats as a color map for plots
  # Seed the NumPy random number generator to ensure reproducible color maps
  numpy.random.seed(1)
  color_map = np.linspace(0., 1., num_fsrs, endpoint=False)
  numpy.random.shuffle(color_map)

  # Make figure of line segments for each track
  fig = plt.figure()

  for i in range(num_segments):

    # Create a color map corresponding to FSR IDs
    jet = cm = plt.get_cmap('jet')
    cNorm  = colors.Normalize(vmin=0, vmax=max(color_map))
    scalarMap = cmx.ScalarMappable(norm=cNorm)
    color = scalarMap.to_rgba(color_map[fsrs[i] % num_fsrs])
    plt.plot([x[i*2], x[i*2+1]], [y[i*2], y[i*2+1]], c=color)

  plt.xlim([x.min(), x.max()])
  plt.ylim([y.min(), y.max()])

  title = 'Segments for ' + str(num_azim) + ' angles and ' + str(spacing) + \
        ' cm spacing'

  plt.title(title)

  filename = directory + 'segments-' + str(num_azim) + '-angles-' + \
      str(spacing) + '-spacing.png'

  fig.savefig(filename, bbox_inches='tight')
  plt.close(fig)



##
# @brief This method takes in a Geometry object and plots a color-coded 2D
#        surface plot representing the Materials in the Geometry.
# @details The Geometry object must be initialized with Materials, Cells,
#          Universes and lattices before being passed into this method. A user
#          may invoke this function from an OpenMOC Python file as follows:
#
# @code
#         openmoc.plotter.plot_materials(geometry)
# @endcode
#
# @param geometry a geometry object which has been initialized with Materials,
#        Cells, Universes and Lattices
# @param gridsize an optional number of grid cells for the plot
# @param xlim optional list/tuple of the minimim/maximum x-coordinates
# @param ylim optional list/tuple of the minimim/maximum y-coordinates
def plot_materials(geometry, gridsize=250, xlim=None, ylim=None):

  global subdirectory

  directory = openmoc.get_output_directory() + subdirectory

  # Make directory if it does not exist
  if not os.path.exists(directory):
    os.makedirs(directory)

  # Error checking
  if not 'Geometry' in str(type(geometry)):
    py_printf('ERROR', 'Unable to plot the Materials since ' + \
                    'input was not a geometry class object')

  if not is_integer(gridsize):
    py_printf('ERROR', 'Unable to plot the Materials since ' + \
              'since the gridsize %d is not an integer', gridsize)

  if gridsize <= 0:
    py_printf('ERROR', 'Unable to plot the Materials ' + \
              'with a negative gridsize (%d)', gridsize)

  py_printf('NORMAL', 'Plotting the materials...')

  # Initialize a NumPy array for the surface colors
  surface = numpy.zeros((gridsize, gridsize), numpy.int64)

  # Retrieve the pixel coordinates
  coords = get_pixel_coords(geometry, gridsize, xlim, ylim)

  # Find the <aterial IDs for each grid point
  for i in range(gridsize):
    for j in range(gridsize):

      x = coords['x'][i]
      y = coords['y'][j]

      point = openmoc.LocalCoords(x, y)
      point.setUniverse(geometry.getRootUniverse())
      cell = geometry.findCellContainingCoords(point)

      # If we did not find a Cell for this region, use a -1 "bad" number color
      if cell is None:
        surface[j][i] = -1
      else:
        surface[j][i] = cell.getMaterial().getId()

  # Get the number of Materials in the Geometry
  materials = geometry.getAllMaterials()
  num_materials = len(materials)

  # Create array of all Material IDs and randomly (but reproducibly) permute it
  material_ids = [material_id for material_id in materials]
  numpy.random.seed(1)
  numpy.random.shuffle(material_ids)

  # Create an array of the colors (array indices) for each value in the surface
  colors = np.zeros((gridsize, gridsize))

  for material_id in np.unique(surface):
    index = material_ids.index(material_id)
    indices = np.where(surface == material_id)
    colors[indices] = index

  # Make Matplotlib color "bad" numbers (ie, NaN, INF) with transparent pixels
  cmap = plt.get_cmap('spectral')
  cmap.set_bad(alpha=0.0)

  # Plot a 2D color map of the Materials
  fig = plt.figure()
  colors = np.flipud(colors)
  plt.imshow(colors, extent=coords['bounds'],
             interpolation='nearest', cmap=cmap, vmin=0, vmax=num_materials)
  plt.title('Materials')
  filename = directory + 'materials.png'
  fig.savefig(filename, bbox_inches='tight')
  plt.close(fig)


##
# @brief This method takes in a Geometry object and plots a color-coded 2D
#        surface plot representing the Cells in the Geometry.
# @details The geometry object must be initialized with Materials, Cells,
#          Universes and Lattices before being passed into this method. A user
#          may invoke this function from an OpenMOC Python file as follows:
#
# @code
#         openmoc.plotter.plot_cells(geometry)
# @endcode
#
# @param geometry a Geometry object which has been initialized with Materials,
#        Cells, Universes and Lattices
# @param gridsize an optional number of grid cells for the plot
# @param xlim optional list/tuple of the minimim/maximum x-coordinates
# @param ylim optional list/tuple of the minimim/maximum y-coordinates
def plot_cells(geometry, gridsize=250, xlim=None, ylim=None):

  global subdirectory

  directory = openmoc.get_output_directory() + subdirectory

  # Make directory if it does not exist
  if not os.path.exists(directory):
    os.makedirs(directory)

  # Error checking
  if not 'Geometry' in str(type(geometry)):
    py_printf('ERROR', 'Unable to plot the Cells since ' + \
              'input was not a Geometry class object')

  if not is_integer(gridsize):
    py_printf('ERROR', 'Unable to plot the Cells since ' + \
                'since the gridsize %d is not an integer', gridsize)

  if gridsize <= 0:
    py_printf('ERROR', 'Unable to plot the Cells ' + \
              'with a negative gridsize (%d)', gridsize)

  py_printf('NORMAL', 'Plotting the cells...')

  # Initialize a NumPy array for the surface colors
  surface = np.zeros((gridsize, gridsize), numpy.int64)

  # Retrieve the pixel coordinates
  coords = get_pixel_coords(geometry, gridsize, xlim, ylim)

  # Find the Cell IDs for each grid point
  for i in range(gridsize):
    for j in range(gridsize):

      x = coords['x'][i]
      y = coords['y'][j]

      point = openmoc.LocalCoords(x, y)
      point.setUniverse(geometry.getRootUniverse())
      cell = geometry.findCellContainingCoords(point)

      # If we did not find a Cell for this region, use a -1 "bad" number color
      if cell is None:
        surface[j][i] = -1
      else:
        surface[j][i] = cell.getId()

  # Get the number of Material Cells in the Geometry
  material_cells = geometry.getAllMaterialCells()
  num_cells = len(material_cells)

  # Create array of all Cell IDs and randomly (but reproducibly) permute it
  cell_ids = [cell_id for cell_id in material_cells]
  numpy.random.seed(1)
  numpy.random.shuffle(cell_ids)

  # Create an array of the colors (array indices) for each value in the surface
  colors = np.zeros((gridsize, gridsize))

  for cell_id in np.unique(surface):
    index = cell_ids.index(cell_id)
    indices = np.where(surface == cell_id)
    colors[indices] = index

  # Make Matplotlib color "bad" numbers (ie, NaN, INF) with transparent pixels
  cmap = plt.get_cmap('spectral')
  cmap.set_bad(alpha=0.0)

  # Plot a 2D color map of the Cells
  fig = plt.figure()
  colors = np.flipud(colors)
  plt.imshow(colors, extent=coords['bounds'],
             interpolation='nearest', cmap=cmap, vmin=0, vmax=num_cells)
  plt.title('Cells')
  filename = directory + 'cells.png'

  fig.savefig(filename, bbox_inches='tight')
  plt.close(fig)



##
# @brief This method takes in a Geometry object and plots a color-coded 2D
#        surface plot representing the flat source regions in the Geometry.
# @details The Geometry object must be initialized with Materials, Cells,
#          Universes and Lattices before being passed into this method. A user
#          may invoke this function from an OpenMOC Python file as follows:
#
# @code
#         openmoc.plotter.plot_flat_source_regions(geometry)
# @endcode
#
# @param geometry a geometry object which has been initialized with Materials,
#        Cells, Universes and Lattices
# @param gridsize an optional number of grid cells for the plot
# @param xlim optional list/tuple of the minimim/maximum x-coordinates
# @param ylim optional list/tuple of the minimim/maximum y-coordinates
def plot_flat_source_regions(geometry, gridsize=250, xlim=None, ylim=None):

  global subdirectory

  directory = openmoc.get_output_directory() + subdirectory

  # Make directory if it does not exist
  if not os.path.exists(directory):
    os.makedirs(directory)

  # Error checking
  if not 'Geometry' in str(type(geometry)):
    py_printf('ERROR', 'Unable to plot the flat source regions since ' + \
              'input was not a geometry class object')

  if not is_integer(gridsize):
    py_printf('ERROR', 'Unable to plot the flat source regions since ' + \
              'since the gridsize %d is not an integer', gridsize)

  if gridsize <= 0:
    py_printf('ERROR', 'Unable to plot the flat source regions ' + \
              'with a negative gridsize (%d)', gridsize)

  py_printf('NORMAL', 'Plotting the flat source regions...')

  # Get the number of flat source regions
  num_fsrs = geometry.getNumFSRs()

  if num_fsrs == 0:
    py_printf('ERROR', 'Unable to plot the flat source regions ' + \
              'since no tracks have been generated.')

  # Initialize a NumPy array for the surface colors
  surface = numpy.zeros((gridsize, gridsize), dtype=np.int64)

  # Retrieve the pixel coordinates
  coords = get_pixel_coords(geometry, gridsize, xlim, ylim)

  # Find the flat source region IDs for each grid point
  for i in range(gridsize):
    for j in range(gridsize):

      x = coords['x'][i]
      y = coords['y'][j]

      local_coords = openmoc.LocalCoords(x, y)
      local_coords.setUniverse(geometry.getRootUniverse())
      geometry.findCellContainingCoords(local_coords)
      fsr_id = geometry.getFSRId(local_coords)

      # If we did not find a region for this region, use a -1 "bad" number color
      if fsr_id is None:
        surface[j][i] = -1
      else:
       surface[j][i] = fsr_id

      del local_coords

  # Replace each Cell ID with a random (but reproducible) color ID
  # NOTE: This color coding scheme only works for FSRs and CMFD cells and not
  # for Materials and Cells. The reason is that FSRs and CMFD cells are by
  # definition a sequence of consecutive, monotonically increasing integers.
  # Material and Cell IDs however may be any sequence of positive integers.
  all_ids = np.arange(num_fsrs, dtype=np.int64)

  id_colors = np.arange(num_fsrs, dtype=np.int64)
  numpy.random.seed(1)
  np.random.shuffle(id_colors)

  ids_to_colors = np.arange(num_fsrs, dtype=np.int64)
  ids_to_colors[all_ids] = id_colors

  colors = ids_to_colors.take(surface)

  # Make Matplotlib color "bad" numbers (ie, NaN, INF) with transparent pixels
  cmap = plt.get_cmap('spectral')
  cmap.set_bad(alpha=0.0)

  # Plot a 2D color map of the flat source regions
  fig = plt.figure()
  colors = np.flipud(colors)
  plt.imshow(colors, extent=coords['bounds'],
             interpolation='nearest', cmap=cmap, vmin=0, vmax=num_fsrs)
  plt.title('Flat Source Regions')
  filename = directory + 'flat-source-regions.png'
  fig.savefig(filename, bbox_inches='tight')
  plt.close(fig)


##
# @brief This method takes in a Geometry and Cmfd object and plots a
#        color-coded 2D surface plot representing the CMFD cells in a geometry.
# @details The Geometry object must be initialized with Materials, Cells,
#          Universes and Lattices before being passed into this method.
#          Plotting the CMFD cells requires that segments must have been
#          created for the geometry and FSR IDs assigned to regions. A user
#          may invoke this function from an OpenMOC Python file as follows:
#
# @code
#         openmoc.plotter.plot_cmfd_cells(geometry, cmfd)
# @endcode
#
# @param geometry a geometry object which has been initialized with Materials,
#        Cells, Universes and Lattices. Segments must have been created or 
#        extracted from a file.
# @param cmfd a Cmfd object which has been used with the geometry in 
#        generating segments. The Cmfd object must have the _overlay_mesh
#        flag set to true; otherwise, the map linking FSR IDs to CMFD cells
#        would not have been created.
# @param gridsize an optional number of grid cells for the plot
# @param xlim optional list/tuple of the minimim/maximum x-coordinates
# @param ylim optional list/tuple of the minimim/maximum y-coordinates
def plot_cmfd_cells(geometry, cmfd, gridsize=250, xlim=None, ylim=None):

  global subdirectory

  directory = openmoc.get_output_directory() + subdirectory

  # Make directory if it does not exist
  if not os.path.exists(directory):
    os.makedirs(directory)

  # Error checking
  if not 'Geometry' in str(type(geometry)):
    py_printf('ERROR', 'Unable to plot the CMFD cells since ' + \
              'input was not a geometry class object')

  if not 'Cmfd' in str(type(cmfd)):
    py_printf('ERROR', 'Unable to plot the CMFD cells since ' + \
              'input was not a CMFD class object')

  if not is_integer(gridsize):
    py_printf('ERROR', 'Unable to plot the CMFD cells since ' + \
              'since the gridsize %s is not an integer', str(gridsize))

  if gridsize <= 0:
    py_printf('ERROR', 'Unable to plot the CMFD cells ' + \
              'with a negative gridsize (%d)', gridsize)

  py_printf('NORMAL', 'Plotting the CMFD cells...')

  # Initialize a NumPy array for the surface colors
  surface = numpy.zeros((gridsize, gridsize), numpy.int64)

  # Retrieve the pixel coordinates
  coords = get_pixel_coords(geometry, gridsize, xlim, ylim)

  # Find the CMFD cell ID for each grid point
  for i in range(gridsize):
    for j in range(gridsize):

      x = coords['x'][i]
      y = coords['y'][j]

      local_coords = openmoc.LocalCoords(x, y)
      local_coords.setUniverse(geometry.getRootUniverse())
      geometry.findCellContainingCoords(local_coords)
      fsr_id = geometry.getFSRId(local_coords)
      cell_id = cmfd.convertFSRIdToCmfdCell(fsr_id)

      # If we did not find a cell for this point, use a -1 "bad" number color
      if np.isnan(cell_id):
        surface[j][i] = -1
      else:
       surface[j][i] = cell_id

  # Get the number of CMFD cells
  num_cmfd_cells = cmfd.getNumCells()

  # Replace each Cell ID with a random (but reproducible) color ID
  # NOTE: This color coding scheme only works for FSRs and CMFD cells and not
  # for Materials and Cells. The reason is that FSRs and CMFD cells are by
  # definition a sequence of consecutive, monotonically increasing integers.
  # Material and Cell IDs however may be any sequence of positive integers.
  all_ids = np.arange(num_cmfd_cells, dtype=np.int64)

  id_colors = np.arange(num_cmfd_cells, dtype=np.int64)
  numpy.random.seed(1)
  np.random.shuffle(id_colors)

  ids_to_colors = np.arange(num_cmfd_cells, dtype=np.int64)
  ids_to_colors[all_ids] = id_colors

  colors = ids_to_colors.take(surface)

  # Make Matplotlib color "bad" numbers (ie, NaN, INF) with transparent pixels
  cmap = plt.get_cmap('spectral')
  cmap.set_bad(alpha=0.0)

  # Plot a 2D color map of the CMFD cells
  fig = plt.figure()
  colors = np.flipud(colors)
  plt.imshow(colors, extent=coords['bounds'],
             interpolation='nearest', cmap=cmap)
  plt.title('CMFD cells')
  filename = directory + 'cmfd-cells.png'
  fig.savefig(filename, bbox_inches='tight')


##
# @brief This method takes in a Solver object and plots a color-coded 2D
#        surface plot representing the flat source region scalar fluxes.
# @details The Solver must have converged the flat source sources prior to
#          calling this routine. A user may invoke this function from an 
#          OpenMOC Python file as follows:
#
# @code
#         openmoc.plotter.plot_spatial_fluxes(solver, energy_groups=[1,7])
# @endcode
#
# @param solver a Solver object that has converged the source for the Geometry
# @param energy_groups a Python list of integer energy groups to plot
# @param gridsize an optional number of grid cells for the plot
# @param xlim optional list/tuple of the minimim/maximum x-coordinates
# @param ylim optional list/tuple of the minimim/maximum y-coordinates
def plot_spatial_fluxes(solver, energy_groups=[1],
                        gridsize=250, xlim=None, ylim=None):

  global subdirectory

  directory = openmoc.get_output_directory() + subdirectory

  # Make directory if it does not exist
  if not os.path.exists(directory):
    os.makedirs(directory)

  if not 'Solver' in str(type(solver)):
    py_printf('ERROR', 'Unable to plot the FSR flux since the ' + \
              'input did not contain a solver class object')

  geometry = solver.getGeometry()
  num_groups = geometry.getNumEnergyGroups()

  if isinstance(energy_groups, (list, tuple, np.ndarray)):
    for group in energy_groups:
      if not is_integer(group):
        py_printf('ERROR', 'Unable to plot the FSR flux since the ' + \
                  'energy_groups contains %s which is not a number', str(group))

      elif group <= 0:
        py_printf('ERROR', 'Unable to plot the FSR flux since the ' + \
                  'energy_groups contains %d which is less than the ' + \
                  'index for all energy groups', group)

      elif group > num_groups:
        py_printf('ERROR', 'Unable to plot the FSR flux since the ' + \
                  'energy_groups contains %d which is greater than ' + \
                  'the index for all energy groups', group)

  else:
    py_printf('ERROR', 'Unable to plot the FSR flux since the ' + \
              'energy_groups is not a Python tuple/list or NumPy array')

  if not is_integer(gridsize):
    py_printf('ERROR', 'Unable to plot the FSR flux since the ' + \
              'gridsize %s is not an integer', str(gridsize))

  if gridsize <= 0:
    py_printf('ERROR', 'Unable to plot the FSR flux with a ' + \
              'negative gridsize (%d)', gridsize)

  py_printf('NORMAL', 'Plotting the FSR scalar fluxes...')

  # Initialize a numpy array for the groupwise scalar fluxes
  fluxes = numpy.zeros((len(energy_groups), gridsize, gridsize))

  # Retrieve the pixel coordinates
  coords = get_pixel_coords(geometry, gridsize, xlim, ylim)

  for i in range(gridsize):
    for j in range(gridsize):

      # Find the flat source region IDs for each grid point
      x = coords['x'][i]
      y = coords['y'][j]

      point = openmoc.LocalCoords(x, y)
      point.setUniverse(geometry.getRootUniverse())
      geometry.findCellContainingCoords(point)
      fsr_id = geometry.getFSRId(point)

      # If we did not find a region for this region, use a -1 "bad" number color
      if np.isnan(fsr_id):
        fluxes[:,j,i] = -1

      # Get the scalar flux for each energy group in this FSR
      else:
        for index, group in enumerate(energy_groups):
          fluxes[index,j,i] = solver.getFSRScalarFlux(fsr_id, group)

  # Loop over all energy group and create a plot
  for index, group in enumerate(energy_groups):

    # Plot a 2D color map of the flat source regions
    fig = plt.figure()
    plt.imshow(np.flipud(fluxes[index,:,:]), extent=coords['bounds'])
    plt.colorbar()
    plt.title('FSR Scalar Flux (Group {0})'.format(group))
    filename = directory + 'fsr-flux-group-' + str(group) + '.png'
    fig.savefig(filename, bbox_inches='tight')
    plt.close(fig)


##
# @brief This method takes in a Solver object and plots the scalar
#        flux vs. energy for one or more flat source regions.
# @details The Solver must have converged the flat source sources prior to
#          calling this routine. The routine will generate a step plot of the
#          flat flux across each energy group. 
#
#          An optional parameter for the energy group bounds may be input. 
#          The group bounds should be input in increasing order of energy.
#          If group bounds are not specified, the routine will use equal 
#          width steps for each energy group.
#
#          A user may invoke this function from an OpenMOC Python file 
#          as follows:
#
# @code
#         openmoc.plotter.plot_energy_fluxes(solver, fsrs=[1,5,20],
#                                            group_bounds=[0., 0.625, 2e7])
# @endcode
#
# @param solver a Solver object that has converged the source for the Geometry
# @param fsrs the flat source region IDs of interest
# @param group_bounds an optional Python list of the energy group bounds (eV)
# @param norm a boolean indicating whether to normalize the flux
# @param loglog boolean indicating whether to plot use a log-log scale
def plot_energy_fluxes(solver, fsrs, group_bounds=None, norm=True, loglog=True):

  global subdirectory

  directory = openmoc.get_output_directory() + subdirectory

  # Make directory if it does not exist
  if not os.path.exists(directory):
    os.makedirs(directory)

  if not 'Solver' in str(type(solver)):
    py_printf('ERROR', 'Unable to plot the flux vs. energy ' + \
              'since input did not contain a Solver class object')

  geometry = solver.getGeometry()
  num_fsrs = geometry.getNumFSRs()
  num_groups = geometry.getNumEnergyGroups()

  if isinstance(fsrs, (tuple, list, np.ndarray)):
    for fsr in fsrs:
      if not is_integer(fsr):
        py_printf('ERROR', 'Unable to plot the flux vs. energy since ' + \
                  'the fsrs contains %s which is not an int', str(fsr))

      elif fsr < 0:
        py_printf('ERROR', 'Unable to plot the flux vs. energy since ' + \
                  'the fsrs contains %d which is less than zero', fsr)

      elif fsr >= num_fsrs:
        py_printf('ERROR', 'Unable to plot the flux vs. energy since ' + \
                  'the fsrs contains %d which is greater than the ' + \
                  'total number of FSRs %d', fsr, num_fsrs)

  else:
    py_printf('ERROR', 'Unable to plot the flux vs. energy since ' + \
              'the fsrs is not a Python tuple, list or NumPy array')

  if isinstance(group_bounds, (tuple, list, np.ndarray)):

    if not all(low < up for low, up in zip(group_bounds, group_bounds[1:])):
      py_printf('ERROR', 'Unable to plot the flux vs. energy since the ' + \
                'energy group bounds are not monotonically increasing')

    elif len(group_bounds) != geometry.getNumEnergyGroups()+1:
      py_printf('ERROR', 'Unable to plot the flux vs. energy since the ' + \
                'group bounds does not correspond to %d groups', num_groups)

    for bound in group_bounds:
      if not is_integer(bound) and not is_float(bound):
        py_printf('ERROR', 'Unable to plot the flux vs. energy since the ' + \
                  'group bounds contains %s which is not a number', str(fsr))

      elif bound < 0:
        py_printf('ERROR', 'Unable to plot the flux vs. energy since the ' + \
                  'group bounds contains %f which is less than zero', bound)

  elif group_bounds is None:
    group_bounds = np.arange(num_groups+1, dtype=np.int)
    loglog = False

  else:
    py_printf('ERROR', 'Unable to plot the flux vs. energy since ' + \
              'the group bounds is not a Python tuple, list or NumPy array')

  py_printf('NORMAL', 'Plotting the scalar fluxes vs. energy...')

  # Compute difference in energy bounds for each group
  group_deltas = np.ediff1d(group_bounds)
  group_bounds = np.flipud(group_bounds)
  group_deltas = np.flipud(group_deltas)
 
  # Iterate over all flat source regions
  for fsr in fsrs:

    # Allocate memory for an array of this FSR's fluxes
    fluxes = np.zeros(num_groups, dtype=np.float)

    # Extract the flux in each energy group
    for group in range(num_groups):
        fluxes[group] = solver.getFSRScalarFlux(fsr, group+1)

    # Normalize fluxes to the total integrated flux
    if norm:
      fluxes /= np.sum(group_deltas * fluxes)

    # Initialize a separate plot for this FSR's fluxes
    fig = plt.figure()

    # Draw horizontal/vertical lines on the plot for each energy group
    for group in range(num_groups):
    
      # Horizontal line
      if loglog:
        plt.loglog(group_bounds[group:group+2], [fluxes[group]]*2, 
                   linewidth=3, c='b', label='openmoc', linestyle='-')
      else:
        plt.plot(group_bounds[group:group+2], [fluxes[group]]*2, 
                 linewidth=3, c='b', label='openmoc', linestyle='-')

      # Vertical lines
      if group < num_groups - 1:
        if loglog:
          plt.loglog([group_bounds[group+1]]*2, fluxes[group:group+2], 
                     c='b', linestyle='--')
        else:
          plt.plot([group_bounds[group+1]]*2, fluxes[group:group+2], 
                   c='b', linestyle='--')

    plt.xlabel('Energy')
    plt.ylabel('Flux')
    plt.xlim((min(group_bounds), max(group_bounds)))
    plt.grid()
    plt.title('FSR {0} Flux ({1} groups)'.format(fsr, num_groups))
    filename = directory + 'flux-fsr-' + str(fsr) + '.png'
    plt.savefig(filename, bbox_inches='tight')
    plt.close(fig)


##
# @brief This method plots a color-coded 2D surface plot representing the 
#        FSR fission rates in the Geometry.
# @details The Solver must have converged the flat source sources prior to
#          calling this routine. The routine will generate a step plot of the
#          flat flux across each energy group. 
#
#          A user may invoke this function from an OpenMOC Python file 
#          as follows:
#
# @code
#         openmoc.plotter.plot_fission_rates(solver)
# @endcode
#
# @param solver a Solver object that has converged the source for the Geometry
# @param gridsize an optional number of grid cells for the plot
# @param xlim optional list/tuple of the minimim/maximum x-coordinates
# @param ylim optional list/tuple of the minimim/maximum y-coordinates
def plot_fission_rates(solver, gridsize=250, xlim=None, ylim=None):

  global subdirectory

  directory = openmoc.get_output_directory() + subdirectory

  # Make directory if it does not exist
  if not os.path.exists(directory):
    os.makedirs(directory)

  if not 'Solver' in str(type(solver)):
    py_printf('ERROR', 'Unable to plot the fission rates ' + \
              'since input did not contain a solver class object')

  if not is_integer(gridsize):
    py_printf('ERROR', 'Unable to plot the fission rates ' + \
              'since since the gridsize %s is not an integer', str(gridsize))

  if gridsize <= 0:
    py_printf('ERROR', 'Unable to plot the fission rates ' + \
              'with a negative gridsize (%d)', gridsize)

  py_printf('NORMAL', 'Plotting the flat source region fission rates...')

  # Get geometry
  geometry = solver.getGeometry()

  # Compute the volume-weighted fission rates for each FSR
  fission_rates = solver.computeFSRFissionRates(geometry.getNumFSRs())

  # Initialize a numpy array of fission rates
  surface = numpy.zeros((gridsize, gridsize))

  # Retrieve the pixel coordinates
  coords = get_pixel_coords(geometry, gridsize, xlim, ylim)

  for i in range(gridsize):
    for j in range(gridsize):

      # Find the flat source region IDs for each grid point
      x = coords['y'][i]
      y = coords['x'][j]

      point = openmoc.LocalCoords(x, y)
      point.setUniverse(geometry.getRootUniverse())
      geometry.findCellContainingCoords(point)
      fsr_id = geometry.getFSRId(point)

      # If we did not find a region for this region, use a -1 "bad" number color
      if np.isnan(fsr_id):
        surface[j][i] = -1
      # Get the fission rate in this FSR
      else:
       surface[j][i] = fission_rates[fsr_id]

  # Plot a 2D color map of the flat source regions fission rates
  fig = plt.figure()
  plt.imshow(np.flipud(surface), extent=coords['bounds'])
  plt.colorbar()
  plt.title('Flat Source Region Fission Rates')
  filename = directory + 'fission-rates.png'
  fig.savefig(filename, bbox_inches='tight')


##
# @brief This is a helper method to define coordinates for a plotting window.
# @details This routine builds a coordinate surface map for the plotting
#          window defined for by the user. If no window was defined, then
#          this routine uses the outer bounding box around the geometry as
#          the plotting window.
# @param geometry a Geometry object which has been initialized with Materials,
#        Cells, Universes and Lattices
# @param gridsize an optional number of grid cells for the plot
# @param xlim optional list/tuple of the minimim/maximum x-coordinates
# @param ylim optional list/tuple of the minimim/maximum y-coordinates
# @return a dictionary with the plotting window map and bounding box
def get_pixel_coords(geometry, gridsize, xlim, ylim):

  # initialize variables to be returned
  bounds = [geometry.getMinX() + TINY_MOVE, geometry.getMaxX() - TINY_MOVE,
            geometry.getMinY() + TINY_MOVE, geometry.getMaxY() - TINY_MOVE]
  xcoords = None
  ycoords = None
  coords = dict()

  if not xlim is None:
    bounds[0] = xlim[0]
    bounds[1] = xlim[1]

  if not ylim is None:
    bounds[2] = ylim[0]
    bounds[3] = ylim[1]

  xcoords = np.linspace(bounds[0], bounds[1], gridsize)
  ycoords = np.linspace(bounds[2], bounds[3], gridsize)

  # add attributes to coords dictionary
  coords['x'] = xcoords
  coords['y'] = ycoords
  coords['bounds'] = bounds

  return coords
