import array
cimport cpython.array as array
from grid cimport Grid
from spatialDependence cimport SpatialDependence

cdef class Wind:
    """Defines wind advection.
    Most importantly, this computes the probability of a mosquito advecting from one active cell to another active cell.
    Either the method parseRawFile() or parseProcessedFile() must be called before this object is ready to provide advection information"""

    # file containing the raw wind data.  The units used should match those used in Grid (eg, km)
    cdef str raw_wind_fn

    # file containing the processed data, which is: given an active cell, what cells it will advect to
    cdef str processed_wind_fn

    # probability distribution of times an advecting mosquito will stay in the air and be advected by the wind.
    # This PDF is of the form [[time0, prob0], [time1, prob1], [time2, prob2], ...].
    # The times must be in the same time unit as the wind data (for instance, days).
    # The probabilities should sum to 1 (prob0 + prob1 + ... = 1).
    # For instance, if all mosquitoes that are being advected by the wind stay in the wind for 0.5 days, then pdf = [[0.5, 1.0]]
    # For instance, if 30% of mosquitoes stay in the wind for 0.1 days, and 70% stay in the wind for 0.4 days, then pdf = [[0.1, 0.3], [0.4, 0.7]]
    # For instance, if all mosquitoes stay in the wind for 0.5 days, BUT an internal timestep of 0.1 days is required to accurately track advective movements, then pdf = [[0.1, 0], [0.2, 0], [0.3, 0], [0.4, 0], [0.5, 1.0]]
    cdef list pdf
    # number of entries in pdf
    cdef unsigned num_pdf

    # whether the wind files are binary (=1) or ascii (=0)
    cdef unsigned binary_file_format

    # the grid (which defines cell_size, etc)
    cdef Grid grid
    
    # lower-left corner
    cdef float xmin, ymin

    # cell side-length
    cdef float cell_size

    # number of cells in x and y directions
    cdef unsigned nx, ny

    # total number of cells
    cdef unsigned num_cells

    # active_index[i] = active cell index of the global cell index i.  Here i ranges from 0 to num_cells - 1, and active_index values will range from 0 to num_active_cells - 1
    cdef array.array active_index

    # number of active cells
    cdef unsigned num_active_cells

    # whether advection_from, etc, has been build (raw_wind_fn or processed_wind_fn has been read) (0 = false, 1 = true)
    cdef int processed_data_computed

    # raw velocity in the x and y direction DIVIDED by cell_size defined in grid (ie, vel_x = number of grid cells/unit time)
    cdef array.array raw_velx
    cdef array.array raw_vely

    # processed advective transfer
    cdef array.array advection_from
    cdef array.array advection_to
    cdef array.array advection_p
    cdef unsigned num_advection

    # unsigned array template to make array creation faster
    cpdef array.array uint_template

    # float array template to make array creation faster
    cpdef array.array float_template

    # the parser that provides the wind information by parsing either raw_wind_fn or processed_wind_fn
    cpdef SpatialDependence windParser

    cpdef parseRawFile(self)
    """Parse the raw wind file specified in the constructor, then processes this data to produce advection_from, etc"""

    cpdef parseProcessedFile(self)
    """Parse the processed wind file specified in the constructor, putting the results in advection_from, etc"""

    cdef void processRawVelocities(self)
    "works out which active cells the advected mosquitoes will end up in"""

    cpdef outputProcessedCSV(self)
    "Outputs the advection information to a file, ready for later reading by parseProcessedFile()"""

    cpdef str getRawWindFilename(self)
    """Returns the raw wind filename"""

    cpdef str getProcessedWindFilename(self)
    """Returns the raw wind filename"""

    cpdef list getPDF(self)
    """Returns the probability distribution for advection"""

    cpdef Grid getGrid(self)
    """Returns the Grid object that this class depends on"""

    cpdef int getProcessedDataComputed(self)
    """Returns 1 if processed advection data has been computed"""

    cpdef array.array getAdvectionFrom(self)
    """Before calling this, check that getProcessedDataComputed() == 1.
    Returns an array containing active cell indices from which advection occurs.
    For each i, getAdvectionFrom()[i] is an active cell index.  Mosquitoes advect from this cell to
    getAdvectionTo()[i] with probability getAdvectionP()[i]"""
    
    cpdef array.array getAdvectionTo(self)
    """Before calling this, check that getProcessedDataComputed() == 1.
    Returns an array containing active cell indices to which advection occurs.
    For each i, getAdvectionFrom()[i] is an active cell index.  Mosquitoes advect from this cell to
    getAdvectionTo()[i] with probability getAdvectionP()[i]"""
    
    cpdef array.array getAdvectionP(self)
    """Before calling this, check that getProcessedDataComputed() == 1.
    Returns an array containing advection probability.active cell indices to which advection occurs.
    For each i, getAdvectionFrom()[i] is an active cell index.  Mosquitoes advect from this cell to
    getAdvectionTo()[i] with probability getAdvectionP()[i]"""

    cpdef unsigned getNumAdvection(self)
    """Returns the size of the advection_from array (= size of advection_to array = size of advection_p array)"""

    cpdef setBinaryFileFormat(self, unsigned value)
    """set binary_file_format to value, which must be 0 (indicating ascii plaintext CSV) or 1 (indicating binary)"""
    
    cpdef unsigned getBinaryFileFormat(self)
    """returns binary_file_format"""
    

