import array
cimport cpython.array as array

cdef class Grid:
    """Defines the grid, most imporantly the active cells and nearest-neighbour connections between them.
    See methods below, in particular internal_global_index and setActiveAndInactive for numbering schemes"""

    # lower-left corner
    cpdef float xmin, ymin

    # cell side-length
    cpdef float cell_size

    # number of cells in x and y directions
    cpdef unsigned nx, ny

    # number of cells
    cpdef unsigned num_cells

    # whether a cell is active (1) or not (0)
    cpdef array.array active

    # number of active cells
    cpdef unsigned num_active_cells

    # global_index[i] = global index of the active cell index i.  Here i ranges from 0 to num_active_cells - 1, and the global_index values will be between 0 and num_cells - 1
    cpdef array.array global_index

    # active[i] = active cell index of the global cell index i.  Here i ranges from 0 to num_cells - 1, and active_index values will range from 0 to num_active_cells - 1
    cpdef array.array active_index

    # connectivity information: connect_from[i] is a cell number that connects to cell number connect_to[i].  These arrays contain only active cells
    cpdef array.array connect_from
    cpdef array.array connect_to

    # total size of connect_from (which contains only the active set)
    cpdef unsigned num_connections

    cdef unsigned internal_global_index(self, unsigned x_ind, unsigned y_ind)
    """provides the index used in this class, given x and y cell number.
    Note that:
    (x_ind, y_ind) = 0 corresponds to the lower left-hand cell, which is internal_global_index = 0
    (x_ind, y_ind) = (nx - 1, 0) corresponds to the lower right-hand cell, which is internal_global_index = nx - 1
    (x_ind, y_ind) = (0, ny - 1) corresponds to the upper left-hand cell, which is internal_global_index = (ny - 1) * nx
    (x_ind, y_ind) = (nx - 1, ny - 1) corresponds to the upper right-hand cell, which is internal_global_index = (ny - 1) * nx + nx - 1
    """

    cpdef setActiveAndInactive(self, filename)
    """Parses filename to determine the active and inactive cells.
    Updates connectivity information, so getConnectionsFrom, etc, only returns active cells.
    filename must have the following format:
    - any empty lines are ignored
    - Data must be preceded by a line of the form '#xmin=a,ymin=b,cell_size=c,nx=d,ny=e' where a,b,c,d,e must match that given in the constructor
    This allows for error checking
    - Any line that starts with # is ignored (except for the '#xmin...' line)
    - Data is CSV form.
    The first row corresponds to the ymin set of cells
    The second row corresponds to the ymin+cell_size set of cells
    The third row corresponds to the ymin+2*cell_size set of cells, etc
    - Each row must have nx comma-separated entries.
    - The entries are either 0 (cell is inactive) or 1 (cell is active)"""
    

    cdef void computeNumActive(self)
    """Computes self.num_active_cells, self.global_index, self_active_index.
    Also, checks if there are entries that are not zero or one"""
    
    cdef void buildAdjacency(self)
    """Builds connect_from and connect_to, and num_connections, based on the self.active"""

    cpdef float getXmin(self)
    "Returns x coord of the lower-left corner"""

    cpdef float getYmin(self)
    "Returns y coord of the lower-left corner"""

    cpdef unsigned getNx(self)
    "Returns number of cells in x direction"""

    cpdef unsigned getNy(self)
    "Returns number of cells in y direction"""

    cpdef unsigned getNumCells(self)
    """Returns the total number of cells"""

    cpdef unsigned getNumActiveCells(self)
    """Returns the number of active cells"""

    cpdef unsigned getNumConnections(self)
    """Returns the number of connections within the active set of cells"""

    cpdef array.array getConnectionsFrom(self)
    """Returns the connections from the active cells"""

    cpdef array.array getConnectionsTo(self)
    """Returns the connections to the active cells"""

    cpdef array.array getGlobalIndex(self)
    """Returns the global index of the active cells"""

    cpdef array.array getActiveIndex(self)
    """Returns the active index of the the global cell array"""

    cpdef outputActiveCSV(self, filename)
    "Outputs the active information to a file"""
