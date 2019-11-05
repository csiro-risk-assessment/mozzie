#   #cython: boundscheck=False, wraparound=False, nonecheck=False
import time
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

    # global_index[i] = global index of the active cell index.  Here i ranges from 0 to num_active_cells - 1, and the global_index values will be between 0 and num_cells - 1
    cpdef array.array global_index

    # active[i] = active cell index of the global cell index.  Here i ranges from 0 to num_cells - 1, and active_index values will range from 0 to num_active_cells - 1
    cpdef array.array active_index

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

        try:
            with open(filename, 'r') as f:
                data = f.readlines()
        except:
            raise IOError('Cannot open or read ' + filename)

        self.active = array.array('I', [])
        checked_header = False
        header_string_error = "The header line in " + filename + " does not match #xmin=" + str(self.xmin) + ",ymin=" + str(self.ymin) + ",cell_size=" + str(self.cell_size) + ",nx=" + str(self.nx) + ",ny=" + str(self.ny)
        data_string_error = "Data in " + filename + " must be CSV formatted.  Each line in the file corresponds to a row (constant y) of cells, so must contain " + str(self.nx) + " entries (separated by commas).  Each entry must be either 0 (inactive) or 1 (active).  There must be " + str(self.ny) + " such rows.  The first row corresponds to cells at y=ymin, the next row at y=ymin+cell_size, etc"
        for line in data:
            if not line.strip():
                continue
            line = line.strip()
            if line.startswith("#xmin"):
                # catch the most common errors
                line = line.split(",")
                if len(line) != 5:
                    raise ValueError(header_string_error)
                if not (line[0].startswith("#xmin=") and line[1].startswith("ymin=") and line[2].startswith("cell_size=") and line[3].startswith("nx=") and line[4].startswith("ny=")):
                    raise ValueError(header_string_error)
                try:
                    xmin = float(line[0].split("=")[1])
                    ymin = float(line[1].split("=")[1])
                    cell_size = float(line[2].split("=")[1])
                    nx = int(line[3].split("=")[1])
                    ny = int(line[4].split("=")[1])
                except:
                    raise ValueError(header_string_error)
                if not (xmin == self.xmin and ymin == self.ymin and cell_size == self.cell_size and nx == int(self.nx) and ny == int(self.ny)):
                    raise ValueError(header_string_error)
                checked_header = True
                continue
            if line.startswith("#"):
                continue
            # can only get here if we should be reading data
            if not checked_header:
                raise ValueError(header_string_error)

            try:
                line = [int(i) for i in line.split(",")]
            except:
                raise ValueError(data_string_error)
            if len(line) != self.nx:
                raise ValueError(data_string_error)
            array.extend(self.active, array.array('I', line))

        if len(self.active) != self.num_cells:
            raise ValueError(data_string_error)

        # build active and adjacency information
        self.computeNumActive()
        self.buildAdjacency()


    cdef void computeNumActive(self):
        """Computes self.num_active_cells, self.global_index.  Also, checks if there are entries that are not zero or one"""

        self.num_active_cells = 0
        self.global_index = array.array('I', [])
        self.active_index = array.array('I', [self.num_cells] * self.num_cells) # use some dummy values
        cdef unsigned i
        for i in range(len(self.active)):
            if self.active[i] == 0:
                continue
            elif self.active[i] == 1:
                array.extend(self.global_index, array.array('I', [i]))
                self.active_index[i] = self.num_active_cells
                self.num_active_cells += 1
            else:
                raise ValueError("Data in active matrix that determines whether a cell is active or inactive must be only 0 or 1, but yours includes the number " + str(self.active[i]))

    cdef void buildAdjacency(self):
        """Builds connect_from and connect_to, and num_connections, based on the self.active"""

        # build connectivity information
        self.connect_from = array.array('I', [])
        self.connect_to = array.array('I', [])
        cdef unsigned x_ind
        cdef unsigned y_ind
        cdef unsigned ind
        # connecting to right neighbour
        for y_ind in range(self.ny):
            for x_ind in range(self.nx - 1):
                ind = self.internal_global_index(x_ind, y_ind)
                if (self.active[ind] == 1 and self.active[ind + 1] == 1):
                    self.connect_from.extend(array.array('I', [self.active_index[ind]]))
                    self.connect_to.extend(array.array('I', [self.active_index[ind + 1]]))
        # connecting to left neighbour
        for y_ind in range(self.ny):
            for x_ind in range(1, self.nx):
                ind = self.internal_global_index(x_ind, y_ind)
                if (self.active[ind] == 1 and self.active[ind - 1] == 1):
                    self.connect_from.extend(array.array('I', [self.active_index[ind]]))
                    self.connect_to.extend(array.array('I', [self.active_index[ind - 1]]))
        # connecting to bottom neighbour
        for y_ind in range(1, self.ny):
            for x_ind in range(self.nx):
                ind = self.internal_global_index(x_ind, y_ind)
                if (self.active[ind] == 1 and self.active[ind - self.nx] == 1):
                    self.connect_from.extend(array.array('I', [self.active_index[ind]]))
                    self.connect_to.extend(array.array('I', [self.active_index[ind - self.nx]]))
        # connecting to top neighbour
        for y_ind in range(self.ny - 1):
            for x_ind in range(self.nx):
                ind = self.internal_global_index(x_ind, y_ind)
                if (self.active[ind] == 1 and self.active[ind + self.nx] == 1):
                    self.connect_from.extend(array.array('I', [self.active_index[ind]]))
                    self.connect_to.extend(array.array('I', [self.active_index[ind + self.nx]]))
        self.num_connections = len(self.connect_from)



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

    cpdef array.array getGlobalIndex(self):
        """Returns the global index of the active cells"""
        return self.global_index

    cpdef array.array getActiveIndex(self):
        """Returns the active index of the the global cell array"""
        return self.active_index

    cpdef outputActiveCSV(self, filename):
        "Outputs the active information to a file"""
        f = open(filename, 'w')
        f.write("#Active cell information written at: " + time.asctime() + "\n")
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

