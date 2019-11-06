import array
cimport cpython.array as array
from cell cimport Cell
from wind cimport Wind
from grid cimport Grid

cdef class Spatial:
    """Holds and manipulates information at all cells"""

    # time step size
    cpdef float dt
    
    # diffusion coefficient for diffusing mosquitoes
    cpdef float diffusion_coeff

    # fraction of population that diffuses to nearest neighbour
    # in comparison to the diffusion equation,
    # diffusion_d = (diffusion_coefficient) * 4 * dt / (dx)^2
    # (the 4 comes from the number of nearest neighbours: evaluating the laplacian on a square grid)
    cpdef float diffusion_d

    # diffusion_d / 4
    cpdef float diff_d

    # Number of populations at a single Cell
    cpdef unsigned num_populations_at_cell

    # Number of diffusing populations at a single Cell
    cpdef unsigned num_diffusing_populations_at_cell

    # Total number of diffusing populations over all cells
    cpdef unsigned num_diffusing_populations_total

    # a list of all the Cells
    cpdef list all_cells

    # change of the diffusing populations.  This is a class variable to avoid repeated allocations every time diffuse() is called
    cdef array.array change_diff

    # all the diffusing populations.  This is a class variable to avoid repeated allocations every time diffuse() is called
    cdef array.array all_diffusing_populations

    # the grid
    cdef Grid grid

    # size-length of the cells
    cdef float cell_size

    # number of active cells
    cdef unsigned num_active_cells

    # connections_from
    cdef array.array connections_from

    # connections_to
    cdef array.array connections_to

    # number of connections
    cdef unsigned num_connections

    def __init__(self, float dt, float diffusion_coeff, Grid grid, list cells):

        self.dt = dt
        self.diffusion_coeff = diffusion_coeff
        self.grid = grid
        self.all_cells = cells

        self.num_active_cells = self.grid.getNumActiveCells()
        if len(self.all_cells) != self.num_active_cells:
            raise ValueError("Incorrect number of cells: " + len(self.all_cells) + " which should be " + self.grid.getNumActiveCells())
        
        self.cell_size = self.grid.getCellSize()
        cdef int num_nearest_neighbours = 4
        self.diffusion_d = self.diffusion_coeff * num_nearest_neighbours * self.dt / (self.cell_size * self.cell_size)
        self.diff_d = self.diffusion_d / num_nearest_neighbours
        
        self.num_populations_at_cell = Cell().getNumberOfPopulations()
        self.num_diffusing_populations_at_cell = Cell().getNumberOfDiffusingPopulations()


        self.num_diffusing_populations_total = self.num_active_cells * self.num_diffusing_populations_at_cell
        # initialize change_diff
        self.change_diff = array.clone(array.array('f', []), self.num_diffusing_populations_total, zero = False)
        # initialize all_diffusing_populations
        self.all_diffusing_populations = array.clone(array.array('f', []), self.num_diffusing_populations_total, zero = False)

        self.num_active_cells = self.grid.getNumActiveCells()
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
        # utility index
        cdef unsigned i

        # grab all the diffusing populations
        for ind in range(self.num_active_cells):
            i = ind * self.num_diffusing_populations_at_cell
            for p in range(self.num_diffusing_populations_at_cell):
                self.all_diffusing_populations.data.as_floats[i + p] = self.all_cells[ind].getDiffusingSingletonNoCheck(p)

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

        # add the result to the cells
        for ind in range(self.num_active_cells):
            i = ind * self.num_diffusing_populations_at_cell
            for p in range(self.num_diffusing_populations_at_cell):
                self.all_cells[ind].addDiffusingSingletonNoCheck(p, self.change_diff.data.as_floats[i + p])
                
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
                    f.write(str(self.all_cells[active_ind].getPopulation()[population]) + ",")
            x_ind = x_max - 1
            ind = self.grid.internal_global_index(x_ind, y_ind)
            active_ind = active[ind]
            if active_ind == num_cells:
                # inactive cell
                f.write("0\n")
            else:
                f.write(str(self.all_cells[active_ind].getPopulation()[population]) + "\n")
        f.close()

        
        
