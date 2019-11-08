import array
cimport cpython.array as array
from wind cimport Wind
from grid cimport Grid
from populationsAndParameters cimport PopulationsAndParameters

cdef class SpatialDynamics:
    """Performs the diffusion and advection over the grid"""

    # the grid
    cdef Grid grid

    # all the quantities (populations and parameters) of all active cells
    cpdef array.array all_quantities

    # number of active cells
    cdef unsigned num_active_cells

    # Number of quantities (populations and parameters) at a single Cell
    cpdef unsigned num_quantities_at_cell

    # size-length of the cells
    cdef float cell_size

    # Number of diffusing populations at a single Cell
    cpdef unsigned num_diffusing_populations_at_cell

    # The diffusing population indices
    cpdef array.array diffusing_indices

    # Total number of diffusing populations over all cells
    cpdef unsigned num_diffusing_populations_total

    # change of the diffusing populations.  This is a class variable to avoid repeated allocations every time diffuse() is called
    cdef array.array change_diff

    # all the diffusing populations.  This is a class variable to avoid repeated allocations every time diffuse() is called
    cdef array.array all_diffusing_populations

    # connections_from
    cdef array.array connections_from

    # connections_to
    cdef array.array connections_to

    # number of connections
    cdef unsigned num_connections

    # Number of advecting populations at a single Cell
    cpdef unsigned num_advecting_populations_at_cell

    # The advecting population indices
    cpdef array.array advecting_indices

    # Total number of advecting populations over all cells
    cpdef unsigned num_advecting_populations_total

    # change of the advecting populations.  This is a class variable to avoid repeated allocations every time advect() is called
    cdef array.array change_adv

    # all the advecting populations.  This is a class variable to avoid repeated allocations every time advect() is called
    cdef array.array all_advecting_populations

    def __init__(self, Grid grid, PopulationsAndParameters pap):

        self.grid = grid
        self.all_quantities = pap.quantities

        cdef unsigned i

        self.num_active_cells = self.grid.getNumActiveCells()
        self.num_quantities_at_cell = pap.getCell().getNumberOfPopulations() + pap.getCell().getNumberOfParameters()
        if len(self.all_quantities) != self.num_active_cells * self.num_quantities_at_cell:
            raise ValueError("Incorrect size of all_quantities: " + str(len(self.all_quantities)) + " which should be the product of " + str(self.num_active_cells) + " and " + str(self.num_quantities_at_cell))
        
        self.cell_size = self.grid.getCellSize()
        
        cdef array.array float_template = array.array('f', [])

        self.num_diffusing_populations_at_cell = pap.getCell().getNumberOfDiffusingPopulations()
        self.diffusing_indices = pap.getCell().getDiffusingIndices()
        self.num_diffusing_populations_total = self.num_active_cells * self.num_diffusing_populations_at_cell
        # initialize change_diff
        self.change_diff = array.clone(float_template, self.num_diffusing_populations_total, zero = False)
        # initialize all_diffusing_populations
        self.all_diffusing_populations = array.clone(float_template, self.num_diffusing_populations_total, zero = False)

        self.connections_from = self.grid.getConnectionsFrom()
        self.connections_to = self.grid.getConnectionsTo()
        self.num_connections = self.grid.getNumConnections()
        # Note(1): In diffuse() we need to multiply by self.num_diffusing_populations_at_cell
        # Note(1): Let's do that now instead
        for i in range(self.num_connections):
            self.connections_from.data.as_uints[i] = self.num_diffusing_populations_at_cell * self.connections_from.data.as_uints[i]
            self.connections_to.data.as_uints[i] = self.num_diffusing_populations_at_cell * self.connections_to.data.as_uints[i]
        
        self.num_advecting_populations_at_cell = pap.getCell().getNumberOfAdvectingPopulations()
        self.advecting_indices = pap.getCell().getAdvectingIndices()
        self.num_advecting_populations_total = self.num_active_cells * self.num_advecting_populations_at_cell
        # initialize change_adv
        self.change_adv = array.clone(float_template, self.num_advecting_populations_total, zero = False)
        # initialize all_advecting_populations
        self.all_advecting_populations = array.clone(float_template, self.num_advecting_populations_total, zero = False)


    cpdef void diffuse(self, float dt, float diffusion_coeff):
        """One timestep of diffusion"""

        # fraction of population that diffuses to nearest neighbour
        # in comparison to the diffusion equation,
        # diffusion_d = (diffusion_coefficient) * 4 * dt / (dx)^2
        # (the 4 comes from the number of nearest neighbours: evaluating the laplacian on a square grid)
        cdef int num_nearest_neighbours = 4
        cpdef float diffusion_d = diffusion_coeff * num_nearest_neighbours * dt / (self.cell_size * self.cell_size)
        # diffusion_d / 4
        cpdef float diff_d = diffusion_d / num_nearest_neighbours

        # active cell index
        cdef unsigned ind
        # population counter
        cdef unsigned p
        # utility indeces
        cdef unsigned i, j, k

        # grab all the diffusing populations
        for ind in range(self.num_active_cells):
            i = ind * self.num_diffusing_populations_at_cell
            j = ind * self.num_quantities_at_cell
            for p in range(self.num_diffusing_populations_at_cell):
                self.all_diffusing_populations.data.as_floats[i + p] = self.all_quantities.data.as_floats[j + self.diffusing_indices.data.as_uints[p]]

        # initialise the self.change_diff in populations, which is just the amount that comes out of the cells
        for i in range(self.num_diffusing_populations_total):
            self.change_diff.data.as_floats[i] = - diffusion_d * self.all_diffusing_populations.data.as_floats[i]

        # disperse diff_d * population to neighbours
        cdef unsigned from_index
        cdef unsigned to_index
        for i in range(self.num_connections):
            # Note(1): This is why the optimisation was performed above.  Now we don't have to do:
            # Note(1): from_index = self.num_diffusing_populations_at_cell * self.connections_from.data.as_uints[i]
            # Note(1): to_index = self.num_diffusing_populations_at_cell * self.connections_to.data.as_uints[i]
            from_index = self.connections_from.data.as_uints[i]
            to_index = self.connections_to.data.as_uints[i]
            for p in range(self.num_diffusing_populations_at_cell):
                self.change_diff.data.as_floats[to_index + p] = self.change_diff.data.as_floats[to_index + p] + diff_d * self.all_diffusing_populations.data.as_floats[from_index + p]

        # add the result to the populations
        for ind in range(self.num_active_cells):
            i = ind * self.num_diffusing_populations_at_cell
            j = ind * self.num_quantities_at_cell
            for p in range(self.num_diffusing_populations_at_cell):
                k = j + self.diffusing_indices.data.as_uints[p]
                self.all_quantities.data.as_floats[k] = self.all_quantities.data.as_floats[k] + self.change_diff.data.as_floats[i + p]
                
    cpdef advect(self, float advection_fraction, Wind wind):
        """One timestep of advection, using the given wind.
        advection_fraction of all populations that the Cell has labelled as 'advecting' will experience advection"""
        if wind.getProcessedDataComputed() != 1:
            raise ValueError("Wind must have been processed before being used")

        cdef array.array afr = wind.getAdvectionFrom()
        cdef array.array ato = wind.getAdvectionTo()
        cdef array.array apr = wind.getAdvectionP()
        cdef unsigned num_advections = wind.getNumAdvection()

        # active cell index
        cdef unsigned ind
        # population counter
        cdef unsigned p
        # utility indeces
        cdef unsigned i, j, k
        # the quantity dumped at the "to" cell
        cdef float qu

        # grab all the advecting populations
        for ind in range(self.num_active_cells):
            i = ind * self.num_advecting_populations_at_cell
            j = ind * self.num_quantities_at_cell
            for p in range(self.num_advecting_populations_at_cell):
                self.all_advecting_populations.data.as_floats[i + p] = self.all_quantities.data.as_floats[j + self.advecting_indices.data.as_uints[p]]

        # initialise the self.change_adv in populations, which is just the amount that comes out of the cells
        for i in range(self.num_advecting_populations_total):
            self.change_adv.data.as_floats[i] = - advection_fraction * self.all_advecting_populations.data.as_floats[i]

        # use wind to disperse to neighbours
        cdef unsigned from_index
        cdef unsigned to_index
        for i in range(num_advections):
            from_index = self.num_advecting_populations_at_cell * afr.data.as_uints[i]
            to_index = self.num_advecting_populations_at_cell * ato.data.as_uints[i]
            qu = advection_fraction * apr.data.as_floats[i]
            for p in range(self.num_advecting_populations_at_cell):
                self.change_adv.data.as_floats[to_index + p] = self.change_adv.data.as_floats[to_index + p] + qu * self.all_advecting_populations.data.as_floats[from_index + p]

        # add the result to the populations
        for ind in range(self.num_active_cells):
            i = ind * self.num_advecting_populations_at_cell
            j = ind * self.num_quantities_at_cell
            for p in range(self.num_advecting_populations_at_cell):
                k = j + self.advecting_indices.data.as_uints[p]
                self.all_quantities.data.as_floats[k] = self.all_quantities.data.as_floats[k] + self.change_adv.data.as_floats[i + p]

        
    cpdef outputCSV(self, str filename, unsigned pop_or_param):
        """Outputs cell information for given population or parameter number to filename in CSV format"""
        if pop_or_param >= self.num_quantities_at_cell:
            raise ValueError("You requested pop_or_param number " + str(pop_or_param) + " but there are only " + str(self.num_quantities_at_cell) + " at each cell")
        # x index, y index
        cdef unsigned x_ind, y_ind
        # max of x
        cdef unsigned x_max = self.grid.getNx()
        # max of y
        cdef unsigned y_max = self.grid.getNy()
        # global cell index, active cell index
        cdef unsigned ind, active_ind
        # total number of cells
        cdef unsigned num_cells = self.grid.getNumCells()
        # active cell index, given global cell index
        cdef array.array active = self.grid.getActiveIndex()
        f = open(filename, "w")
        for y_ind in range(y_max):
            for x_ind in range(x_max - 1):
                ind = self.grid.internal_global_index(x_ind, y_ind)
                active_ind = active[ind]
                if active_ind == num_cells:
                    # inactive cell
                    f.write("0,")
                else:
                    f.write(str(self.all_quantities[active_ind * self.num_quantities_at_cell + pop_or_param]) + ",")
            x_ind = x_max - 1
            ind = self.grid.internal_global_index(x_ind, y_ind)
            active_ind = active[ind]
            if active_ind == num_cells:
                # inactive cell
                f.write("0\n")
            else:
                f.write(str(self.all_quantities[active_ind * self.num_quantities_at_cell + pop_or_param]) + "\n")
        f.close()

        
        
