import os
import time
from grid cimport Grid
import array
cimport cpython.array as array

cdef class Wind:

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

    cdef array.array active_index

    def __init__(self, str raw_wind_fn, str processed_wind_fn, list pdf, Grid grid):

        self.raw_wind_fn = raw_wind_fn
        self.processed_wind_fn = processed_wind_fn

        pdf_error_string = "PDF defining the advecting velocities in terms of the raw velocity must be of the form [[t0, p0], [t1, p1], [t2, p2]...], where p is the probability that an advecting mosquito will stay in the wind-stream for time t.  The t's must be positive and monotonically increasing.  The p's must sum to unity"
        try:
            self.pdf = []
            prev_time = 0.0
            for e in pdf:
                this_time = float(e[0])
                if this_time > prev_time:
                    self.pdf.append([this_time - prev_time, float(e[1])])
                    prev_time = this_time
                else:
                    raise ValueError(pdf_error_string)
            self.num_pdf = len(self.pdf)
        except:
            raise ValueError(pdf_error_string)
        if sum([p[1] for p in self.pdf]) > 1.0001 or sum([p[1] for p in self.pdf]) < 0.9999:
            raise ValueError(pdf_error_string)
            
        self.grid = grid
        self.xmin = self.grid.getXmin()
        self.ymin = self.grid.getYmin()
        self.cell_size = self.grid.getCellSize()
        self.nx = self.grid.getNx()
        self.ny = self.grid.getNy()
        self.num_cells = self.grid.getNumCells()
        self.active_index = self.grid.getActiveIndex()
        self.processed_data_computed = 0

        self.raw_velx = array.array('f', [])
        self.raw_vely = array.array('f', [])

        self.advection_from = array.array('I', [])
        self.advection_to = array.array('I', [])
        self.advection_p = array.array('f', [])
        self.num_advection = 0

    cpdef parseRawFile(self):
        """Parse the raw wind file specified in the constructor, then processes this data to produce advection_from, etc"""

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
                line = [float(vel_comp) for vel_comp in line.split(",")]
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

        self.processRawVelocities()
        self.processed_data_computed = 1

    cpdef parseProcessedFile(self):
        """Parse the processed wind file specified in the constructor, putting the results in advection_from, etc"""

        try:
            with open(self.processed_wind_fn, 'r') as f:
                data = f.readlines()
        except:
            raise IOError('Cannot open or read ' + self.processed_wind_fn)

        # size the arrays to the maximum possible length (there will be header lines in data)
        self.advection_from = array.array('I', [0] * len(data))
        self.advection_to = array.array('I', [0] * len(data))
        self.advection_p = array.array('f', [0.0] * len(data))
        self.num_advection = 0

        cdef unsigned afr, ato
        cdef float app

        checked_xmin_header = False
        checked_active_header = False
        checked_raw_fn = False
        checked_pdf = False
        data_string_error = "Data in " + self.processed_wind_fn + " must be CSV formatted.  Each line must contain integer1,integer2,float , where integer1 = active cell index from where mosquitoes are advecting; integer2 = active cell index to which mosquitoes are advecting; float = probability of this occuring"
        for line in data:
            if not line.strip():
                continue
            line = line.strip()
            if line.startswith("#xmin"):
                self.grid.checkHeaderLine(self.processed_wind_fn, line)
                checked_xmin_header = True
                continue
            if line.startswith("#Active"):
                if line != "#Active cells defined by file " + self.grid.getActiveFilename():
                    raise ValueError("Active cell filename in " + self.processed_wind_fn + " incompatible with that specified in the grid, which is " + self.grid.getActiveFilename())
                checked_active_header = True
                continue
            if line.startswith("#raw_vel_filename"):
                if line != "#raw_vel_filename=" + os.path.basename(self.raw_wind_fn):
                    raise ValueError("Raw velocity filename in " + self.processed_wind_fn + " incompatible with that specified in the wind constructor, which is " + os.path.basename(self.raw_wind_fn))
                checked_raw_fn = True
                continue
            if line.startswith("#processed_pdf"):
                if line != "#processed_pdf=" + str(self.pdf):
                    raise ValueError("PDF specified in " + self.processed_wind_fn + " incompatible with that specified in the wind constructor.  Remember the file's version is processed to give dt values rather than t values")
                checked_pdf = True
                continue
            if line.startswith("#"):
                continue

            # can only get here if we should be reading data
            if not (checked_xmin_header and checked_active_header and checked_raw_fn and checked_pdf):
                raise ValueError("Not all header lines found in " + self.processed_wind_fn)

            try:
                line = line.split(",")
                afr, ato, app = int(line[0]), int(line[1]), float(line[2])
            except:
                raise ValueError(data_string_error)
            if len(line) != 3:
                raise ValueError(data_string_error)

            self.advection_from[self.num_advection] = afr
            self.advection_to[self.num_advection] = ato
            self.advection_p[self.num_advection] = app
            self.num_advection += 1

        array.resize(self.advection_from, self.num_advection)
        array.resize(self.advection_to, self.num_advection)
        array.resize(self.advection_p, self.num_advection)

        self.processed_data_computed = 1

    cdef void processRawVelocities(self):
        "works out which active cells the advected mosquitoes will end up in"""
        # size the arrays to the maximum they're ever going to be
        self.advection_from = array.array('I', [0] * self.num_cells * self.num_pdf)
        self.advection_to = array.array('I', [0] * self.num_cells * self.num_pdf)
        self.advection_p = array.array('f', [0.0] * self.num_cells * self.num_pdf)
        self.num_advection = 0
        
        cdef unsigned x_ind, y_ind # x and y grid indices
        cdef float velx, vely # velocity
        cdef float xnew, ynew # position, as a float
        cdef float ti, pr # time and probability
        cdef unsigned from_ind, to_ind # global index from which mosquito is originally advecting, and index to which it has currently advected
        cdef unsigned active_from_ind, active_to_ind # active index from which mosquito is advecting, and index to which it has currently advected
        cdef int newx, newy # (x_ind, y_ind) of cell to which the mozzie will advect
        cdef unsigned maybe_new_ind # used to check whether a mosquito has actually entered a new cell
        cdef int has_advected = 0 # if ==1 then mosquito has advected somewhere (could be to original cell)
        cdef float xmin = -0.5 # minimum value of x allowable before mosquito has exited grid
        cdef float ymin = -0.5 # minimum value of y allowable before mosquito has exited grid
        cdef float xmax = float(self.nx - 1) + 0.5 # maximum value of x allowable before mosquito has exited grid
        cdef float ymax = float(self.ny - 1) + 0.5 # maximum value of y allowable before mosquito has exited grid
        for y_ind in range(self.ny):
            for x_ind in range(self.nx):
                has_advected = 0
                from_ind = self.grid.internal_global_index(x_ind, y_ind)
                active_from_ind = self.active_index[from_ind]
                if active_from_ind == self.num_cells:
                    # this is not an active cell: ignore it
                    continue
                # remember that we have divided raw velocities by cell_size, so
                xnew = float(x_ind)
                ynew = float(y_ind)
                # this makes the loop over self.pdf work
                to_ind = from_ind
                active_to_ind = active_from_ind
                for (ti, pr) in self.pdf:
                    # evaluate the velocity at the current position of the mosquito
                    velx = self.raw_velx[to_ind]
                    vely = self.raw_vely[to_ind]
                    # update its real-numbered position
                    xnew += velx * ti
                    ynew += vely * ti
                    if (xnew <= xmin or xnew >= xmax or ynew <= ymin or ynew >= ymax):
                        # mosquito has advected outside grid: assume it has completely gone forever
                        break
                    # round to the nearest cell
                    newx = int(xnew + 0.5)
                    newy = int(ynew + 0.5)
                    # find the index and the active index of the cell that it's now in
                    maybe_new_ind = self.grid.internal_global_index(newx, newy)
                    if maybe_new_ind == to_ind and has_advected == 1:
                        # mosquito has not actually entered a new cell, so if an active cell, just add to the previous probability
                        if active_to_ind != self.num_cells:
                            self.advection_p[self.num_advection - 1] += pr
                    else:
                         # mosquito has enetered new cell (or has_advected == 0, so we need to increment self.num_advected)
                        has_advected = 1
                        to_ind = maybe_new_ind
                        active_to_ind = self.active_index[to_ind]
                        # if active then record it
                        if active_to_ind != self.num_cells:
                            self.advection_from[self.num_advection] = active_from_ind
                            self.advection_to[self.num_advection] = active_to_ind
                            self.advection_p[self.num_advection] = pr
                            self.num_advection += 1
        array.resize(self.advection_from, self.num_advection)
        array.resize(self.advection_to, self.num_advection)
        array.resize(self.advection_p, self.num_advection)
            

    cpdef outputProcessedCSV(self):
        "Outputs the active information to a file"""
        f = open(self.processed_wind_fn, 'w')
        f.write("#Processed wind advection written at: " + time.asctime() + "\n")
        f.write("#Active cells defined by file " + self.grid.getActiveFilename() + "\n")
        f.write("#xmin=" + str(self.xmin) + ",ymin=" + str(self.ymin) + ",cell_size=" + str(self.cell_size) + ",nx=" + str(self.nx) + ",ny=" + str(self.ny) + "\n")
        f.write("#raw_vel_filename=" + os.path.basename(self.raw_wind_fn) + "\n")
        f.write("#processed_pdf=" + str(self.pdf) + "\n")
        f.write("#Data is of the form:\n")
        f.write("#active_cell_id_from,active_cell_id_to,probability\n")
        cdef unsigned ind
        for ind in range(self.num_advection):
            f.write(str(self.advection_from[ind]) + "," + str(self.advection_to[ind]) + "," + str(self.advection_p[ind]) + "\n")
        f.close()

    cpdef str getRawWindFilename(self):
        """Returns the raw wind filename"""
        return self.raw_wind_fn

    cpdef str getProcessedWindFilename(self):
        """Returns the raw wind filename"""
        return self.processed_wind_fn

    cpdef list getPDF(self):
        """Returns the probability distribution for advection"""
        return self.pdf

    cpdef Grid getGrid(self):
        """Returns the Grid object that this class depends on"""
        return self.Grid

    cpdef int getProcessedDataComputed(self):
        """Returns 1 if processed advection data has been computed"""
        return self.processed_data_computed

    cpdef array.array getAdvectionFrom(self):
        """Before calling this, check that getFileRead() == 1.
        Returns an array containing active cell indices from which advection occurs.
        For each i, getAdvectionFrom()[i] is an active cell index.  Mosquitoes advect from this cell to
        getAdvectionTo()[i] with probability getAdvectionP()[i]"""
        return self.advection_from
    
    cpdef array.array getAdvectionTo(self):
        """Before calling this, check that getFileRead() == 1.
        Returns an array containing active cell indices to which advection occurs.
        For each i, getAdvectionFrom()[i] is an active cell index.  Mosquitoes advect from this cell to
        getAdvectionTo()[i] with probability getAdvectionP()[i]"""
        return self.advection_to
    
    cpdef array.array getAdvectionP(self):
        """Before calling this, check that getFileRead() == 1.
        Returns an array containing advection probability.active cell indices to which advection occurs.
        For each i, getAdvectionFrom()[i] is an active cell index.  Mosquitoes advect from this cell to
        getAdvectionTo()[i] with probability getAdvectionP()[i]"""
        return self.advection_p

    cpdef unsigned getNumAdvection(self):
        """Returns the size of the advection_from array (= size of advection_to array = size of advection_p array)"""
        return self.num_advection
    

