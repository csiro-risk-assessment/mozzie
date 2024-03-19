import os
import time
from grid cimport Grid
cimport csvparser
import array
cimport cpython.array as array
from libc.stdlib cimport malloc, free

cdef class Wind:
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
            
        self.binary_file_format = 0
        self.grid = grid
        self.xmin = self.grid.getXmin()
        self.ymin = self.grid.getYmin()
        self.cell_size = self.grid.getCellSize()
        self.nx = self.grid.getNx()
        self.ny = self.grid.getNy()
        self.num_cells = self.grid.getNumCells()
        self.active_index = self.grid.getActiveIndex()
        self.num_active_cells = self.grid.getNumActiveCells()
        self.processed_data_computed = 0

        self.num_advection = 0

        self.uint_template = array.array('I', [])
        self.float_template = array.array('f', [])

    cpdef parseRawFile(self):
        """Parse the raw wind file specified in the constructor, then processes this data to produce advection_from, etc"""
        self.windParser = SpatialDependence(self.xmin, self.ymin, self.cell_size, self.nx, self.ny)
        if self.binary_file_format == 0:
            self.windParser.parse(self.raw_wind_fn, "wind_raw", [])
        else:
            self.windParser.parse(self.raw_wind_fn, "wind_raw_binary", [])
        self.raw_velx = self.windParser.getData0()
        self.raw_vely = self.windParser.getData1()
        cdef float scale = 1.0 / self.cell_size
        cdef unsigned ind
        for ind in range(self.num_cells):
            self.raw_velx.data.as_floats[ind] = scale * self.raw_velx.data.as_floats[ind]
            self.raw_vely.data.as_floats[ind] = scale * self.raw_vely.data.as_floats[ind]
        self.processRawVelocities()
        self.processed_data_computed = 1

    cpdef parseProcessedFile(self):
        """Parse the processed wind file specified in the constructor, putting the results in advection_from, etc"""
        self.windParser = SpatialDependence(self.xmin, self.ymin, self.cell_size, self.nx, self.ny)
        if self.binary_file_format == 0:
            self.windParser.parse(self.processed_wind_fn, "wind_processed", ["#Active cells defined by file " + self.grid.getActiveFilename(), "#raw_vel_filename=" + os.path.basename(self.raw_wind_fn), "#processed_pdf=" + str(self.pdf)])
        else:
            self.windParser.parse(self.processed_wind_fn, "wind_processed_binary", ["#Active cells defined by file " + self.grid.getActiveFilename(), "#raw_vel_filename=" + os.path.basename(self.raw_wind_fn), "#processed_pdf=" + str(self.pdf)])
        self.advection_from = self.windParser.getData0()
        self.advection_to = self.windParser.getData1()
        self.advection_p = self.windParser.getData2()
        self.num_advection = len(self.advection_from)

        cdef unsigned ind
        for ind in range(self.num_advection):
            if self.advection_from.data.as_uints[ind] >= self.num_active_cells or self.advection_to.data.as_uints[ind] >= self.num_active_cells or self.advection_p.data.as_floats[ind] > 1.0 or self.advection_p.data.as_floats[ind] < 0.0:
                raise ValueError("Data in " + self.processed_wind_fn + " is incorrectly bounded.  Bad data line number = " + str(ind))

        self.processed_data_computed = 1

    cdef void processRawVelocities(self):
        "works out which active cells the advected mosquitoes will end up in"""
        # size the arrays to the maximum they're ever going to be
        self.advection_from = array.clone(self.uint_template, self.num_cells * self.num_pdf, zero = False)
        self.advection_to = array.clone(self.uint_template, self.num_cells * self.num_pdf, zero = False)
        self.advection_p = array.clone(self.float_template, self.num_cells * self.num_pdf, zero = False)
        self.num_advection = 0

        cdef array.array dts = array.array('f', [e[0] for e in self.pdf])
        cdef array.array prs = array.array('f', [e[1] for e in self.pdf])
        cdef unsigned pdf_ind
        
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
                active_from_ind = self.active_index.data.as_uints[from_ind]
                if active_from_ind == self.num_cells:
                    # this is not an active cell: ignore it
                    continue
                # remember that we have divided raw velocities by cell_size, so
                xnew = float(x_ind)
                ynew = float(y_ind)
                # this makes the loop over self.pdf work
                to_ind = from_ind
                active_to_ind = active_from_ind
                for pdf_ind in range(self.num_pdf):
                    ti = dts.data.as_floats[pdf_ind]
                    pr = prs.data.as_floats[pdf_ind]
                    # evaluate the velocity at the current position of the mosquito
                    velx = self.raw_velx.data.as_floats[to_ind]
                    vely = self.raw_vely.data.as_floats[to_ind]
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
                            self.advection_p.data.as_floats[self.num_advection - 1] = self.advection_p.data.as_floats[self.num_advection - 1] + pr
                    else:
                         # mosquito has enetered new cell (or has_advected == 0, so we need to increment self.num_advected)
                        has_advected = 1
                        to_ind = maybe_new_ind
                        active_to_ind = self.active_index.data.as_uints[to_ind]
                        # if active then record it
                        if active_to_ind != self.num_cells:
                            self.advection_from.data.as_uints[self.num_advection] = active_from_ind
                            self.advection_to.data.as_uints[self.num_advection] = active_to_ind
                            self.advection_p.data.as_floats[self.num_advection] = pr
                            self.num_advection += 1
        array.resize(self.advection_from, self.num_advection)
        array.resize(self.advection_to, self.num_advection)
        array.resize(self.advection_p, self.num_advection)
            

    cpdef outputProcessedCSV(self):
        "Outputs the active information to a file"""
        cdef unsigned ind
        cdef float* float_data = NULL
        cdef unsigned* uint = NULL
        header_output = ""
        header_output += "#Processed wind advection written at: " + time.asctime() + "\n"
        header_output += "#Active cells defined by file " + self.grid.getActiveFilename() + "\n"
        header_output += "#xmin=" + str(self.xmin) + ",ymin=" + str(self.ymin) + ",cell_size=" + str(self.cell_size) + ",nx=" + str(self.nx) + ",ny=" + str(self.ny) + "\n"
        header_output += "#raw_vel_filename=" + os.path.basename(self.raw_wind_fn) + "\n"
        header_output += "#processed_pdf=" + str(self.pdf) + "\n"
        header_output += "#Data is of the form:\n"
        header_output += "#active_cell_id_from,active_cell_id_to,probability\n"
        if (self.binary_file_format == 0):
            f = open(self.processed_wind_fn, 'w')
            f.write(header_output)
            for ind in range(self.num_advection):
                f.write(str(self.advection_from.data.as_uints[ind]) + "," + str(self.advection_to.data.as_uints[ind]) + "," + str(self.advection_p.data.as_floats[ind]) + "\n")
            f.close()
        else:
            uint = <unsigned *> malloc(2 * self.num_advection * sizeof(unsigned))
            if not uint:
                raise MemoryError()
            float_data = <float *> malloc(self.num_advection * sizeof(float))
            if not float_data:
                raise MemoryError()
            for ind in range(self.num_advection):
                uint[2 * ind] = self.advection_from.data.as_uints[ind]
                uint[2 * ind + 1] = self.advection_to.data.as_uints[ind]
                float_data[ind] = self.advection_p.data.as_floats[ind]
            err = csvparser.writeProcessedWind(self.processed_wind_fn.encode(), header_output.encode(), uint, float_data, self.num_advection, 1)
            free(float_data)
            free(uint)

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
        return self.grid

    cpdef int getProcessedDataComputed(self):
        """Returns 1 if processed advection data has been computed"""
        return self.processed_data_computed

    cpdef array.array getAdvectionFrom(self):
        """Before calling this, check that getProcessedDataComputed() == 1.
        Returns an array containing active cell indices from which advection occurs.
        For each i, getAdvectionFrom()[i] is an active cell index.  Mosquitoes advect from this cell to
        getAdvectionTo()[i] with probability getAdvectionP()[i]"""
        return self.advection_from
    
    cpdef array.array getAdvectionTo(self):
        """Before calling this, check that getProcessedDataComputed() == 1.
        Returns an array containing active cell indices to which advection occurs.
        For each i, getAdvectionFrom()[i] is an active cell index.  Mosquitoes advect from this cell to
        getAdvectionTo()[i] with probability getAdvectionP()[i]"""
        return self.advection_to
    
    cpdef array.array getAdvectionP(self):
        """Before calling this, check that getProcessedDataComputed() == 1.
        Returns an array containing advection probability.active cell indices to which advection occurs.
        For each i, getAdvectionFrom()[i] is an active cell index.  Mosquitoes advect from this cell to
        getAdvectionTo()[i] with probability getAdvectionP()[i]"""
        return self.advection_p

    cpdef unsigned getNumAdvection(self):
        """Returns the size of the advection_from array (= size of advection_to array = size of advection_p array)"""
        return self.num_advection
    

    cpdef setBinaryFileFormat(self, unsigned value):
        if not (value == 0 or value == 1):
            raise ValueError("FileFormat must be 0 or 1")
        self.binary_file_format = value

    cpdef unsigned getBinaryFileFormat(self):
        return self.binary_file_format

        
