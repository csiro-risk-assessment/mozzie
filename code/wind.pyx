import time
from grid cimport Grid
import array
cimport cpython.array as array

cdef class Wind:

    # file containing the raw wind data
    cdef str raw_wind_fn

    # file containing the processed data, which is: given an active cell, what cells it will advect to
    cdef str processed_wind_fn

    # probability distribution: given a wind speed, this gives the probability of an advecting mozzie having a certain speed.
    # This PDF is of the form [[prob0, frac0], [prob1, frac1], [prob2, frac2], ...], where the probability
    # of a mosquito tagged as 'advecting', advecting by frac * wind_speed is prob.  The probs should sum to 1.0.
    # For instance, if advecting mosquitoes always advect at exactly the 0.5 * wind speed, pdf = [[1.0, 0.5]]
    # For instance, if 30% of advecting mosquitoes always advect 0.1 * wind speed, and the remainder at 0.8 * windspeed, pdf = [[0.3, 0.1], [0.7, 0.8]]
    cdef list pdf
    # number of entries in pdf
    cdef unsigned num_pdf

    # the grid (which defines cell_size, etc)
    cdef Grid grid
    
    # lower-left corner
    cpdef float xmin, ymin

    # cell side-length
    cpdef float cell_size

    # number of cells in x and y directions
    cpdef unsigned nx, ny

    # total number of cells
    cpdef unsigned num_cells

    # whether raw_wind_fn or processed_wind_fn has been read (0 = false, 1 = true)
    cpdef int fn_read

    # raw velocity in the x and y direction DIVIDED by cell_size defined in grid (ie, vel_x = number of grid cells/unit time)
    cpdef array.array raw_velx
    cpdef array.array raw_vely

    # processed advective transfer
    cpdef array.array advection_from
    cpdef array.array advection_to
    cpdef array.array advection_p
    cpdef unsigned num_advection

    cpdef array.array active_index

    def __init__(self, str raw_wind_fn, str processed_wind_fn, list pdf, Grid grid):

        self.raw_wind_fn = raw_wind_fn
        self.processed_wind_fn = processed_wind_fn

        try:
            self.pdf = [(float(e[0]), float(e[1])) for e in pdf]
            self.num_pdf = len(self.pdf)
        except:
            raise ValueError("PDF defining the advecting velocities in terms of the raw velocity must be of the form [[p0, f0], [p1, f1], [p2, f2]...], where p is the probability of a mosquito advecting by speed f*raw_speed.  The p's must sum to unity")
        if not (sum([p[0] for p in self.pdf]) > 0.999 and sum([p[0] for p in self.pdf]) < 1.001):
            raise ValueError("PDF defining the advecting velocities in terms of the raw velocity must be of the form [[p0, f0], [p1, f1], [p2, f2]...], where p is the probability of a mosquito advecting by speed f*raw_speed.  The p's must sum to unity")
            
        self.grid = grid
        self.xmin = self.grid.getXmin()
        self.ymin = self.grid.getYmin()
        self.cell_size = self.grid.getCellSize()
        self.nx = self.grid.getNx()
        self.ny = self.grid.getNy()
        self.num_cells = self.grid.getNumCells()
        self.active_index = self.grid.getActiveIndex()
        self.fn_read = 0

    cpdef void parseRawFile(self):
        """Parse the raw wind file specified in the constructor"""

        # utility indices
        cdef unsigned i, ind, active_ind
        # x and y indices
        cdef unsigned x_ind, y_ind
        
        try:
            with open(self.raw_wind_fn, 'r') as f:
                data = f.readlines()
        except:
            raise IOError('Cannot open or read ' + self.raw_wind_fn)

        self.raw_velx = array.array('f', [0.0] * self.grid.getNumCells())
        self.raw_vely = array.array('f', [0.0] * self.grid.getNumCells())
        checked_header = False
        data_string_error = "Data in " + self.raw_wind_fn + " must be CSV formatted.  Each line in the file corresponds to a row (constant y) of cells, so must contain " + str(2 * self.nx) + " entries (separated by commas).  Entries come (vel_x, vel_y) pairs, so a row is of the form vx0,vy0,vx1,vy1,vx2,vy2,....  There must be " + str(self.ny) + " such rows.  The first row corresponds to cells at y=ymin, the next row at y=ymin+cell_size, etc"
        y_ind = 0
        for line in data:
            if not line.strip():
                continue
            line = line.strip()
            if line.startswith("#xmin"):
                self.grid.checkHeaderLine(self.raw_wind_fn, line)
                checked_header = True
                continue
            if line.startswith("#"):
                continue
            # can only get here if we should be reading data
            if not checked_header:
                raise ValueError("Header line of the form #xmin=... not found in " + self.raw_wind_fn)

            try:
                line = [float(i) for i in line.split(",")]
            except:
                raise ValueError(data_string_error)
            if len(line) != 2 * self.nx:
                raise ValueError(data_string_error)

            if y_ind >= self.ny:
                raise ValueError(data_string_error)
            for i in range(self.nx):
                self.raw_velx[y_ind * self.nx + i] = line[2 * i] / self.cell_size
                self.raw_vely[y_ind * self.nx + i] = line[2 * i + 1] / self.cell_size
            y_ind += 1
                
        if y_ind != self.ny:
            raise ValueError(data_string_error)

        # now work out which active cells are the advected mosquitoes will end up in
        advection_from = array.array('I', [0] * self.num_cells * self.num_pdf)
        advection_to = array.array('I', [0] * self.num_cells * self.num_pdf)
        advection_p = array.array('f', [0.0] * self.num_cells * self.num_pdf)
        self.num_advection = 0
        cdef float velx, vely
        cdef int newx, newy # (x_ind, y_ind) of cell to which the mozzie will advect
        cdef unsigned to_ind # internal global index of cell to which mozzie will advect
        cdef unsigned active_to_ind # active index of cell to which mozzie will advect
        for y_ind in range(self.ny):
            for x_ind in range(self.nx):
                ind = self.grid.internal_global_index(x_ind, y_ind)
                active_ind = self.active_index[ind]
                if active_ind == self.num_cells:
                    # this is not an active cell: ignore it
                    continue
                velx = self.raw_velx[ind]
                vely = self.raw_vely[ind]
                for (p, f) in self.pdf:
                    # round to the nearest cell
                    newx = x_ind + int(velx * f + 0.5) 
                    newy = y_ind + int(vely * f + 0.5)
                    if (newx < 0 or newx >= self.nx or newy < 0 or newy >= self.ny):
                        # advects outside grid
                        continue
                    to_ind = self.grid.internal_global_index(newx, newy)
                    active_to_ind = self.active_index[to_ind]
                    if active_to_ind == self.num_cells:
                        # this is not an active cell: ignore it
                        continue
                    self.advection_from[self.num_advection] = active_ind
                    self.advection_to[self.num_advection] = active_to_ind
                    self.advection_p[self.num_advection] = p
                    self.num_advection += 1
                    
                    


            
        self.fn_read = 1
    

