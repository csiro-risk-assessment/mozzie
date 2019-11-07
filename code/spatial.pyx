import array
cimport cpython.array as array
from cell cimport Cell
from wind cimport Wind
from grid cimport Grid
from populations cimport Populations

cdef class Spatial:
    """Holds and manipulates information at all cells"""

    # time step size
    cpdef float dt
    
    # diffusion coefficient for diffusing mosquitoes
    cpdef float diffusion_coeff

    # the grid
    cdef Grid grid

    # all the populations (of all active cells)
    cpdef array.array all_pops

    # number of active cells
    cdef unsigned num_active_cells

    # Number of populations at a single Cell
    cpdef unsigned num_populations_at_cell

    # size-length of the cells
    cdef float cell_size

    # fraction of population that diffuses to nearest neighbour
    # in comparison to the diffusion equation,
    # diffusion_d = (diffusion_coefficient) * 4 * dt / (dx)^2
    # (the 4 comes from the number of nearest neighbours: evaluating the laplacian on a square grid)
    cpdef float diffusion_d

    # diffusion_d / 4
    cpdef float diff_d

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

    def __init__(self, float dt, float diffusion_coeff, Grid grid, Populations all_pops):

        self.dt = dt
        self.diffusion_coeff = diffusion_coeff
        self.grid = grid
        self.all_pops = all_pops.pop

        self.num_active_cells = self.grid.getNumActiveCells()
        self.num_populations_at_cell = Cell().getNumberOfPopulations()
        if len(self.all_pops) != self.num_active_cells * self.num_populations_at_cell:
            raise ValueError("Incorrect size of all_pops: " + str(len(self.all_pops)) + " which should be the product of " + str(self.num_active_cells) + " and " + str(self.num_populations_at_cell))
        
        self.cell_size = self.grid.getCellSize()
        cdef int num_nearest_neighbours = 4
        self.diffusion_d = self.diffusion_coeff * num_nearest_neighbours * self.dt / (self.cell_size * self.cell_size)
        self.diff_d = self.diffusion_d / num_nearest_neighbours
        
        self.num_diffusing_populations_at_cell = Cell().getNumberOfDiffusingPopulations()
        self.diffusing_indices = Cell().getDiffusingIndices()
        self.num_diffusing_populations_total = self.num_active_cells * self.num_diffusing_populations_at_cell
        # initialize change_diff
        cdef array.array float_template = array.array('f', [])
        self.change_diff = array.clone(float_template, self.num_diffusing_populations_total, zero = False)
        # initialize all_diffusing_populations
        self.all_diffusing_populations = array.clone(float_template, self.num_diffusing_populations_total, zero = False)

        self.connections_from = self.grid.getConnectionsFrom()
        self.connections_to = self.grid.getConnectionsTo()
        self.num_connections = self.grid.getNumConnections()
        # Note(1): In diffuse() we need to multiply by self.num_diffusing_populations_at_cell
        # Note(1): Let's do that now instead
        cdef unsigned i
        for i in range(self.num_connections):
            self.connections_from.data.as_uints[i] = self.num_diffusing_populations_at_cell * self.connections_from.data.as_uints[i]
            self.connections_to.data.as_uints[i] = self.num_diffusing_populations_at_cell * self.connections_to.data.as_uints[i]
        

    cpdef void diffuse(self):
        """One timestep of diffusion"""

        # active cell index
        cdef unsigned ind
        # population counter
        cdef unsigned p
        # utility indeces
        cdef unsigned i, j, k

        # grab all the diffusing populations
        for ind in range(self.num_active_cells):
            i = ind * self.num_diffusing_populations_at_cell
            j = ind * self.num_populations_at_cell
            for p in range(self.num_diffusing_populations_at_cell):
                self.all_diffusing_populations.data.as_floats[i + p] = self.all_pops.data.as_floats[j + self.diffusing_indices.data.as_uints[p]]

        # initialise the self.change_diff in populations, which is just the amount that comes out of the cells
        for i in range(self.num_diffusing_populations_total):
            self.change_diff.data.as_floats[i] = - self.diffusion_d * self.all_diffusing_populations.data.as_floats[i]

        # do the Brownian motion
        cdef unsigned from_index
        cdef unsigned to_index
        for i in range(self.num_connections):
            # Note(1): This is why the optimisation was performed above.  Now we don't have to do:
            # Note(1): from_index = self.num_diffusing_populations_at_cell * self.connections_from.data.as_uints[i]
            # Note(1): to_index = self.num_diffusing_populations_at_cell * self.connections_to.data.as_uints[i]
            from_index = self.connections_from.data.as_uints[i]
            to_index = self.connections_to.data.as_uints[i]
            for p in range(self.num_diffusing_populations_at_cell):
                self.change_diff.data.as_floats[to_index + p] = self.change_diff.data.as_floats[to_index + p] + self.diff_d * self.all_diffusing_populations.data.as_floats[from_index + p]

        # add the result to the populations
        for ind in range(self.num_active_cells):
            i = ind * self.num_diffusing_populations_at_cell
            j = ind * self.num_populations_at_cell
            for p in range(self.num_diffusing_populations_at_cell):
                k = j + self.diffusing_indices.data.as_uints[p]
                self.all_pops.data.as_floats[k] = self.all_pops.data.as_floats[k] + self.change_diff.data.as_floats[i + p]
                
    cpdef outputCSV(self, str filename, unsigned population):
        """Outputs cell information for given population number to filename in CSV format"""
        if population >= self.num_populations_at_cell:
            raise ValueError("You requested population number " + str(population) + " but there are only " + str(self.num_populations_at_cell) + " at each cell")
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
                    f.write(str(self.all_pops[active_ind * self.num_populations_at_cell + population]) + ",")
            x_ind = x_max - 1
            ind = self.grid.internal_global_index(x_ind, y_ind)
            active_ind = active[ind]
            if active_ind == num_cells:
                # inactive cell
                f.write("0\n")
            else:
                f.write(str(self.all_pops[active_ind * self.num_populations_at_cell + population]) + "\n")
        f.close()

        
        
