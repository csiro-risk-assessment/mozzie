#   #cython: boundscheck=False, wraparound=False, nonecheck=False
import os
import time

# maximum number of nearest neighbours that are possible for each cell
cdef unsigned max_num_nn = 4

cdef class Grid:

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
        # This should never be resized, without eventually making its size = self.num_cells, for we use the unsafe access to the raw C pointer
        self.active = array.array('I', [1] * self.num_cells)
        self.active_filename = "None"

        # initialise active_index to dummy values "num_cells" which indicates active_index[i] is "not an active cell"
        # This should never be resized, for we use the unsafe access to the raw C pointer
        self.active_index = array.array('I', [self.num_cells] * self.num_cells) # use some dummy values

        # build active and adjacency information
        self.computeNumActive()
        self.buildAdjacency()


    cdef unsigned internal_global_index(self, unsigned x_ind, unsigned y_ind):
        """provides the index used in this class, given x and y cell number.
        Note that:
        (x_ind, y_ind) = 0 corresponds to the lower left-hand cell, which is internal_global_index = 0
        (x_ind, y_ind) = (nx - 1, 0) corresponds to the lower right-hand cell, which is internal_global_index = nx - 1
        (x_ind, y_ind) = (0, ny - 1) corresponds to the upper left-hand cell, which is internal_global_index = (ny - 1) * nx
        (x_ind, y_ind) = (nx - 1, ny - 1) corresponds to the upper right-hand cell, which is internal_global_index = (ny - 1) * nx + nx - 1
        """
        return y_ind * self.nx + x_ind

    cpdef setActiveAndInactive(self, filename):
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

        self.activeParser = SpatialDependence(self.xmin, self.ymin, self.cell_size, self.nx, self.ny)
        self.activeParser.parse(filename, "active_inactive", [])
        self.active = self.activeParser.getData0()
        self.active_filename = os.path.basename(filename)

        # build active and adjacency information
        self.computeNumActive()
        self.buildAdjacency()


    cdef void computeNumActive(self):
        """Computes self.num_active_cells, self.global_index.  Also, checks if there are entries that are not zero or one"""

        self.num_active_cells = 0
        self.global_index = array.array('I', [0] * self.num_cells) # maximum size it can be: below we resize it smaller
        cdef unsigned i
        for i in range(len(self.active)):
            if self.active.data.as_uints[i] == 0:
                self.active_index.data.as_uints[i] = self.num_cells
            elif self.active.data.as_uints[i] == 1:
                self.global_index.data.as_uints[self.num_active_cells] = i
                self.active_index.data.as_uints[i] = self.num_active_cells
                self.num_active_cells += 1
            else:
                raise ValueError("Data in active matrix that determines whether a cell is active or inactive must be only 0 or 1, but yours includes the number " + str(self.active[i]))
        array.resize(self.global_index, self.num_active_cells)

    cdef void buildAdjacency(self):
        """Builds connect_from and connect_to, and num_connections, based on the self.active"""

        # initialise connectivity information to a very large array.  This is for efficiency (can use data.as_uints below) and we resize appropriately at the very end
        self.connect_from = array.array('I', [0] * self.num_cells * max_num_nn)
        self.connect_to = array.array('I', [0] * self.num_cells * max_num_nn)
        self.num_connections = 0
        cdef unsigned x_ind
        cdef unsigned y_ind
        cdef unsigned ind
        # connecting to right neighbour
        for y_ind in range(self.ny):
            for x_ind in range(self.nx - 1):
                ind = self.internal_global_index(x_ind, y_ind)
                if (self.active.data.as_uints[ind] == 1 and self.active.data.as_uints[ind + 1] == 1):
                    self.connect_from.data.as_uints[self.num_connections] = self.active_index.data.as_uints[ind]
                    self.connect_to.data.as_uints[self.num_connections] = self.active_index.data.as_uints[ind + 1]
                    self.num_connections += 1
        # connecting to left neighbour
        for y_ind in range(self.ny):
            for x_ind in range(1, self.nx):
                ind = self.internal_global_index(x_ind, y_ind)
                if (self.active.data.as_uints[ind] == 1 and self.active.data.as_uints[ind - 1] == 1):
                    self.connect_from.data.as_uints[self.num_connections] = self.active_index.data.as_uints[ind]
                    self.connect_to.data.as_uints[self.num_connections] = self.active_index.data.as_uints[ind - 1]
                    self.num_connections += 1
        # connecting to bottom neighbour
        for y_ind in range(1, self.ny):
            for x_ind in range(self.nx):
                ind = self.internal_global_index(x_ind, y_ind)
                if (self.active.data.as_uints[ind] == 1 and self.active.data.as_uints[ind - self.nx] == 1):
                    self.connect_from.data.as_uints[self.num_connections] = self.active_index.data.as_uints[ind]
                    self.connect_to.data.as_uints[self.num_connections] = self.active_index.data.as_uints[ind - self.nx]
                    self.num_connections += 1
        # connecting to top neighbour
        for y_ind in range(self.ny - 1):
            for x_ind in range(self.nx):
                ind = self.internal_global_index(x_ind, y_ind)
                if (self.active.data.as_uints[ind] == 1 and self.active.data.as_uints[ind + self.nx] == 1):
                    self.connect_from.data.as_uints[self.num_connections] = self.active_index.data.as_uints[ind]
                    self.connect_to.data.as_uints[self.num_connections] = self.active_index.data.as_uints[ind + self.nx]
                    self.num_connections += 1
        array.resize(self.connect_from, self.num_connections)
        array.resize(self.connect_to, self.num_connections)



    cpdef str getActiveFilename(self):
        "Returns filename that defined the active/inactive cells"""
        return self.active_filename
    

    cpdef float getXmin(self):
        "Returns x coord of the lower-left corner"""
        return self.xmin

    cpdef float getYmin(self):
        "Returns y coord of the lower-left corner"""
        return self.ymin

    cpdef float getCellSize(self):
        "Returns cell side length"""
        return self.cell_size

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

    cpdef array.array getGlobalIndex(self):
        """Returns the global index of the active cells"""
        return self.global_index

    cpdef array.array getActiveIndex(self):
        """Returns the active index of the the global cell array"""
        return self.active_index

    cpdef outputActiveCSV(self, str filename):
        "Outputs the active information to a file"""
        f = open(filename, 'w')
        f.write("#Active cell information written at: " + time.asctime() + "\n")
        f.write("#Active cells defined by file " + self.active_filename + "\n")
        f.write("#xmin=" + str(self.xmin) + ",ymin=" + str(self.ymin) + ",cell_size=" + str(self.cell_size) + ",nx=" + str(self.nx) + ",ny=" + str(self.ny) + "\n")
        cdef unsigned x_ind
        cdef unsigned y_ind
        cdef unsigned ind
        for y_ind in range(self.ny):
            for x_ind in range(self.nx - 1):
                ind = self.internal_global_index(x_ind, y_ind)
                f.write(str(self.active[ind]) + ",")
            x_ind = self.nx - 1
            ind = self.internal_global_index(x_ind, y_ind)
            f.write(str(self.active[ind]) + "\n")
        f.close()
