#   #cython: boundscheck=False, wraparound=False, nonecheck=False
import array
cimport cpython.array as array

cdef class Grid:
    """Defines the grid"""

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

    # connectivity information: connect_from[i] is a cell number that connects to cell number connect_to[i].  These arrays contain only active cells
    cpdef array.array connect_from
    cpdef array.array connect_to

    # total size of connect_from (which contains only the active set)
    cpdef unsigned num_connections

    def __init__(self, float xmin, float ymin, float cell_size, unsigned nx, unsigned ny):

        if cell_size <= 0.0:
            raise ValueError("Cell size must be positive")

        self.xmin = xmin
        self.ymin = ymin
        self.cell_size = cell_size
        self.nx = nx
        self.ny = ny

        self.num_cells = nx * ny

        # initialise activity information to "all active"
        self.active = array.array('I', [1] * self.num_cells)
        self.num_active_cells = nx * ny

        # build connectivity information
        self.connect_from = array.array('I', [])
        self.connect_to = array.array('I', [])
        cdef unsigned x_ind
        cdef unsigned y_ind
        cdef unsigned ind
        # connecting to right neighbour
        for y_ind in range(self.ny):
            for x_ind in range(self.nx - 1):
                ind = self.internal_index(x_ind, y_ind)
                self.connect_from.extend(array.array('I', [ind]))
                self.connect_to.extend(array.array('I', [ind + 1]))
        # connecting to left neighbour
        for y_ind in range(self.ny):
            for x_ind in range(1, self.nx):
                ind = self.internal_index(x_ind, y_ind)
                self.connect_from.extend(array.array('I', [ind]))
                self.connect_to.extend(array.array('I', [ind - 1]))
        # connecting to bottom neighbour
        for y_ind in range(1, self.ny):
            for x_ind in range(self.nx):
                ind = self.internal_index(x_ind, y_ind)
                self.connect_from.extend(array.array('I', [ind]))
                self.connect_to.extend(array.array('I', [ind - self.nx]))
        # connecting to top neighbour
        for y_ind in range(self.ny - 1):
            for x_ind in range(self.nx):
                ind = self.internal_index(x_ind, y_ind)
                self.connect_from.extend(array.array('I', [ind]))
                self.connect_to.extend(array.array('I', [ind + self.nx]))
        self.num_connections = len(self.connect_from)
            

    cdef unsigned internal_index(self, unsigned x_ind, unsigned y_ind):
        """provides the index used in this class, given x and y cell number.
        Note that:
        (x_ind, y_ind) = 0 corresponds to the lower left-hand cell, which is internal_index = 0
        (x_ind, y_ind) = (nx - 1, 0) corresponds to the lower right-hand cell, which is internal_index = nx - 1
        (x_ind, y_ind) = (0, ny - 1) corresponds to the upper left-hand cell, which is internal_index = (ny - 1) * nx
        (x_ind, y_ind) = (nx - 1, ny - 1) corresponds to the upper right-hand cell, which is internal_index = (ny - 1) * nx + nx - 1
        """
        return y_ind * self.nx + x_ind

    cpdef float getXmin(self):
        "Returns x coord of the lower-left corner"""
        return self.xmin

    cpdef float getYmin(self):
        "Returns y coord of the lower-left corner"""
        return self.ymin

    cpdef unsigned getNx(self):
        "Returns number of cells in x direction"""
        return self.nx

    cpdef unsigned getNy(self):
        "Returns number of cells in y direction"""
        return self.ny

    cpdef unsigned getNumCells(self):
        """Returns the total number of cells"""
        return self.num_cells

    cpdef unsigned getNumActiveCells(self):
        """Returns the number of active cells"""
        return self.num_active_cells

    cpdef unsigned getNumConnections(self):
        """Returns the number of connections within the active set of cells"""
        return self.num_connections

    cpdef array.array getConnectionsFrom(self):
        """Returns the connections from the active cells"""
        return self.connect_from

    cpdef array.array getConnectionsTo(self):
        """Returns the connections to the active cells"""
        return self.connect_to


