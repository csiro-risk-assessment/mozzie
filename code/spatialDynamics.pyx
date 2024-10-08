# Copyright (c) 2024 Commonwealth Scientific and Industrial Research Organisation (CSIRO) ABN 41 687 119 230.
import time
import array
import sys
cimport cpython.array as array
from wind cimport Wind
from grid cimport Grid
from populationsAndParameters cimport PopulationsAndParameters
from cellDynamics cimport CellDynamicsBase

cdef class SpatialDynamics:
    """Performs the diffusion and advection over the grid"""

    # the grid
    cdef Grid grid

    # a reference to the cell in the PopulationsAndParameters
    cdef CellDynamicsBase cell

    # all the quantities (populations and parameters) of all active cells
    cdef array.array all_quantities

    # all the birth terms (by sex, genotype and species) of all active cells
    cdef array.array birth_quantities

    # all the Y' terms (by sex, genotype and species) of all active cells
    cdef array.array yyp_quantities

    # number of active cells
    cdef unsigned num_active_cells

    # Number of quantities (populations and parameters) at a single Cell
    cdef unsigned num_quantities_at_cell
    
    # Number of quantities (populations and parameters) at all cells
    cdef unsigned num_quantities_total

    # size-length of the cells
    cdef float cell_size

    # Number of diffusing populations at a single Cell
    cdef unsigned num_diffusing_populations_at_cell

    # The diffusing population indices
    cdef array.array diffusing_indices

    # Total number of diffusing populations over all cells
    cdef unsigned num_diffusing_populations_total

    # change of the diffusing populations.  This is a class variable to avoid repeated allocations every time diffuse() is called
    cdef array.array change_diff

    # all the diffusing populations.  This is a class variable to avoid repeated allocations every time diffuse() is called
    cdef array.array all_diffusing_populations

    # diffusion_coefficient * num_nearest_neighbours * dt / cell_size^2.  This is a class variable to avoid repeated allocations every time diffuse() is called
    cdef array.array diffusion_d

    # diffusion_coefficient * dt / cell_size^2.  This is a class variable to avoid repeated allocations every time diffuse() is called
    cdef array.array diff_d

    # connections_from
    cdef array.array connections_from

    # connections_to
    cdef array.array connections_to

    # number of connections
    cdef unsigned num_connections

    # Number of advecting populations at a single Cell
    cdef unsigned num_advecting_populations_at_cell

    # The advecting population indices
    cdef array.array advecting_indices

    # The advection classes of the advecting populations
    cdef array.array advection_classes

    # Total number of advecting populations over all cells
    cdef unsigned num_advecting_populations_total

    # change of the advecting populations.  This is a class variable to avoid repeated allocations every time advect() is called
    cdef array.array change_adv

    # all the advecting populations.  This is a class variable to avoid repeated allocations every time advect() is called
    cdef array.array all_advecting_populations

    # a C array for passing to the cell evolve method
    cdef array.array cell_params_and_props
    # a C array for rapid passing to the cell evolve method
    cdef float[:] c_cell_params_and_props

    def __init__(self, Grid grid, PopulationsAndParameters pap):

        self.grid = grid
        self.cell = pap.getCell()
        self.all_quantities = pap.quantities
        cdef unsigned i

        self.num_active_cells = self.grid.getNumActiveCells()
        self.num_quantities_at_cell = self.cell.getNumberOfPopulations() + self.cell.getNumberOfParameters()
        self.num_quantities_total = self.num_active_cells * self.num_quantities_at_cell
        
        if len(self.all_quantities) != self.num_quantities_total:
            raise ValueError("Incorrect size of all_quantities: " + str(len(self.all_quantities)) + " which should be the product of " + str(self.num_active_cells) + " and " + str(self.num_quantities_at_cell))
        
        self.cell_size = self.grid.getCellSize()
        
        cdef array.array float_template = array.array('f', [])
        self.num_diffusing_populations_at_cell = self.cell.getNumberOfDiffusingPopulations()
        self.diffusing_indices = self.cell.getDiffusingIndices()
        self.num_diffusing_populations_total = self.num_active_cells * self.num_diffusing_populations_at_cell
        self.birth_quantities = array.clone(array.array('f', []), self.num_diffusing_populations_total, zero = False)
        self.yyp_quantities = array.clone(array.array('f', []), self.num_diffusing_populations_total, zero = False)

        # initialize change_diff
        self.change_diff = array.clone(float_template, self.num_diffusing_populations_total, zero = False)
        # initialize all_diffusing_populations
        self.all_diffusing_populations = array.clone(float_template, self.num_diffusing_populations_total, zero = False)
        # initialize diffusion_d and diff_d
        self.diffusion_d = array.clone(float_template, self.num_diffusing_populations_at_cell, zero = False)
        self.diff_d = array.clone(float_template, self.num_diffusing_populations_at_cell, zero = False)

        self.connections_from = array.copy(self.grid.getConnectionsFrom())
        self.connections_to = array.copy(self.grid.getConnectionsTo())
        self.num_connections = self.grid.getNumConnections()
        # Note(1): In diffuse() we need to multiply by self.num_diffusing_populations_at_cell
        # Note(1): Let's do that now instead
        for i in range(self.num_connections):
            self.connections_from.data.as_uints[i] = self.num_diffusing_populations_at_cell * self.connections_from.data.as_uints[i]
            self.connections_to.data.as_uints[i] = self.num_diffusing_populations_at_cell * self.connections_to.data.as_uints[i]
        
        self.num_advecting_populations_at_cell = self.cell.getNumberOfAdvectingPopulations()
        self.advecting_indices = self.cell.getAdvectingIndices()
        self.advection_classes = self.cell.getAdvectionClass()
        self.num_advecting_populations_total = self.num_active_cells * self.num_advecting_populations_at_cell
        # initialize change_adv
        self.change_adv = array.clone(float_template, self.num_advecting_populations_total, zero = False)
        # initialize all_advecting_populations
        self.all_advecting_populations = array.clone(float_template, self.num_advecting_populations_total, zero = False)

        self.cell_params_and_props = array.clone(float_template, self.num_quantities_at_cell, zero = False)
        self.c_cell_params_and_props = self.cell_params_and_props

    cpdef array.array getBirthQuantities(self):
        return self.birth_quantities

    cpdef array.array getYypQuantities(self):
        return self.yyp_quantities

    cpdef setBirthQuantities(self, list birth_quantities):
        if len(birth_quantities) != self.num_diffusing_populations_total:
            raise ValueError("Incorrect size of birth_quantities, " + str(len(birth_quantities)) + ", which should be the total number of diffusing populations, " + str(self.num_diffusing_populations_total))
        for i in range(self.num_diffusing_populations_total):
            self.birth_quantities.data.as_floats[i] = birth_quantities[i]

    cpdef setYypQuantities(self, list yyp_quantities):
        if len(yyp_quantities) != self.num_diffusing_populations_total:
            raise ValueError("Incorrect size of yyp_quantities, " + str(len(yyp_quantities)) + ", which should be the total number of diffusing populations, " + str(self.num_diffusing_populations_total))
        for i in range(self.num_diffusing_populations_total):
            self.yyp_quantities.data.as_floats[i] = yyp_quantities[i]    

    cpdef array.array getAllQuantities(self):
        return self.all_quantities

    cpdef setAllQuantities(self, list all_quantities):
        if len(all_quantities) != self.num_quantities_total:
            raise ValueError("Incorrect size of all_quantities, " + str(len(all_quantities)) + ", which should be the total number of quantities, " + str(self.num_quantities_total))
        for i in range(self.num_quantities_total):
            self.all_quantities.data.as_floats[i] = all_quantities[i]

    cpdef void diffuse(self, float dt, float diffusion_coeff):
        """One timestep of diffusion"""

        # fraction of population that diffuses to nearest neighbour
        # in comparison to the diffusion equation,
        # diffusion_d = (diffusion_coefficient) * 4 * dt / (dx)^2
        # (the 4 comes from the number of nearest neighbours: evaluating the laplacian on a square grid)
        cdef int num_nearest_neighbours = 4
        cdef unsigned p
        for p in range(self.num_diffusing_populations_at_cell):
            self.diffusion_d.data.as_floats[p] = diffusion_coeff * num_nearest_neighbours * dt / (self.cell_size * self.cell_size)
            self.diff_d.data.as_floats[p] = diffusion_coeff * dt / (self.cell_size * self.cell_size)

        self.diffuseAll()

    cpdef diffuseVarying(self, float dt, float [:] diffusion_coeffs):
        """One timestep of diffusion with diffusion_coeffs that depend on population type (eg, sex)
        NOTE: diffusion_coeffs"""

        if len(diffusion_coeffs) != self.num_diffusing_populations_at_cell:
            raise ValueError("diffusion_coeffs incorrectly sized")
        # fraction of population that diffuses to nearest neighbour
        # in comparison to the diffusion equation,
        # diffusion_d = (diffusion_coefficient) * 4 * dt / (dx)^2
        # (the 4 comes from the number of nearest neighbours: evaluating the laplacian on a square grid)
        cdef int num_nearest_neighbours = 4
        cdef unsigned p
        for p in range(self.num_diffusing_populations_at_cell):
            self.diffusion_d.data.as_floats[p] = diffusion_coeffs[p] * num_nearest_neighbours * dt / (self.cell_size * self.cell_size)
            self.diff_d.data.as_floats[p] = diffusion_coeffs[p] * dt / (self.cell_size * self.cell_size)
        self.diffuseAll()


    cpdef void diffuseAll(self):
        """One timestep of diffusion"""

        # active cell index
        cdef unsigned ind
        # population counter
        cdef unsigned p
        # utility indices
        cdef unsigned i, j, k

        # grab all the diffusing populations
        for ind in range(self.num_active_cells):
            i = ind * self.num_diffusing_populations_at_cell
            j = ind * self.num_quantities_at_cell
            for p in range(self.num_diffusing_populations_at_cell):
                self.all_diffusing_populations.data.as_floats[i + p] = self.all_quantities.data.as_floats[j + self.diffusing_indices.data.as_uints[p]]

        # initialise the self.change_diff in populations, which is just the amount that comes out of the cells
        for ind in range(self.num_active_cells):
            i = ind * self.num_diffusing_populations_at_cell
            for p in range(self.num_diffusing_populations_at_cell):
                self.change_diff.data.as_floats[i + p] = - self.diffusion_d.data.as_floats[p] * self.all_diffusing_populations.data.as_floats[i + p]

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
                self.change_diff.data.as_floats[to_index + p] = self.change_diff.data.as_floats[to_index + p] + self.diff_d.data.as_floats[p] * self.all_diffusing_populations.data.as_floats[from_index + p]

        # add the result to the populations
        for ind in range(self.num_active_cells):
            i = ind * self.num_diffusing_populations_at_cell
            j = ind * self.num_quantities_at_cell
            for p in range(self.num_diffusing_populations_at_cell):
                k = j + self.diffusing_indices.data.as_uints[p]
                self.all_quantities.data.as_floats[k] = self.all_quantities.data.as_floats[k] + self.change_diff.data.as_floats[i + p]
                
                
    cpdef advect(self, float [:] advection_fraction, Wind wind):
        """One timestep of advection, using the given wind.
        advection_fraction of all populations that the Cell has labelled as 'advecting' will experience advection.  advection_class = i will experience advection_fraction[i]"""
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
        # utility indices
        cdef unsigned i, j, k

        # grab all the advecting populations
        for ind in range(self.num_active_cells):
            i = ind * self.num_advecting_populations_at_cell
            j = ind * self.num_quantities_at_cell
            for p in range(self.num_advecting_populations_at_cell):
                self.all_advecting_populations.data.as_floats[i + p] = self.all_quantities.data.as_floats[j + self.advecting_indices.data.as_uints[p]]

        # initialise the self.change_adv in populations, which is just the amount that comes out of the cells
        for ind in range(self.num_active_cells):
            i = ind * self.num_advecting_populations_at_cell
            for p in range(self.num_advecting_populations_at_cell):
                self.change_adv.data.as_floats[i + p] = - advection_fraction[self.advection_classes.data.as_uints[p]] * self.all_advecting_populations.data.as_floats[i + p]

        # use wind to disperse to neighbours
        cdef unsigned from_index
        cdef unsigned to_index
        for i in range(num_advections):
            from_index = self.num_advecting_populations_at_cell * afr.data.as_uints[i]
            to_index = self.num_advecting_populations_at_cell * ato.data.as_uints[i]
            for p in range(self.num_advecting_populations_at_cell):
                self.change_adv.data.as_floats[to_index + p] = self.change_adv.data.as_floats[to_index + p] + advection_fraction[self.advection_classes.data.as_uints[p]] * apr.data.as_floats[i] * self.all_advecting_populations.data.as_floats[from_index + p]

        # add the result to the populations
        for ind in range(self.num_active_cells):
            i = ind * self.num_advecting_populations_at_cell
            j = ind * self.num_quantities_at_cell
            for p in range(self.num_advecting_populations_at_cell):
                k = j + self.advecting_indices.data.as_uints[p]
                self.all_quantities.data.as_floats[k] = self.all_quantities.data.as_floats[k] + self.change_adv.data.as_floats[i + p]


    cpdef evolveCells(self, float timestep):
        """Evolves all the cells in the grid using the method defined in self.cell"""

        # active cell index
        cdef unsigned ind
        # population index
        cdef unsigned p
        # utility index
        cdef unsigned i, j

        for ind in range(self.num_active_cells):
            i = ind * self.num_quantities_at_cell
            j = ind * self.num_diffusing_populations_at_cell
            # copy into the c_cell_params_and_props local array
            for p in range(self.num_quantities_at_cell):
                self.c_cell_params_and_props[p] = self.all_quantities.data.as_floats[i + p]
            # evolve
            self.cell.evolve(timestep, self.c_cell_params_and_props)
            # copy back
            for p in range(self.num_quantities_at_cell):
                self.all_quantities.data.as_floats[i + p] = self.c_cell_params_and_props[p]
            getBMethod = getattr(self.cell, "getB", None)
            getYYprimeMethod = getattr(self.cell, "getYYprime", None)
            if callable(getBMethod) and callable(getYYprimeMethod):
                b = self.cell.getB()
                yyp = self.cell.getYYprime()
                for p in range(self.num_diffusing_populations_at_cell):
                    self.birth_quantities.data.as_floats[j + p] = b[p]
                    self.yyp_quantities.data.as_floats[j + p] = yyp[p]

    cpdef calcQm(self):
        """Calculates qm at all the cells in the grid, assuming that the populations on the grid are at equilibrium"""

        # active cell index
        cdef unsigned ind
        # population index
        cdef unsigned p
        # utility index
        cdef unsigned i
        # number of populations
        cdef unsigned num_pops = self.cell.getNumberOfPopulations()
        # qm values are held here
        cdef array.array qm_vals = array.clone(array.array('f', []), self.cell.getNumSpecies(), zero = False)

        for ind in range(self.num_active_cells):
            i = ind * self.num_quantities_at_cell
            # copy into the c_cell_params_and_props local array
            for p in range(self.num_quantities_at_cell):
                self.c_cell_params_and_props[p] = self.all_quantities.data.as_floats[i + p]
            # calculate qm at this cell
            qm_vals = self.cell.calcQm(self.c_cell_params_and_props)
            # place qm into the correct slots
            for p in range(num_pops, self.num_quantities_at_cell):
                self.all_quantities.data.as_floats[i + p] = qm_vals.data.as_floats[p - num_pops]

    cpdef outputCSV(self, str filename, int pop_or_param, str inactive_value, str additional_header_lines):
        """Outputs cell information for given population or parameter number to filename in CSV format.
        If pop_or_param >= 0 then it indicates the population (or parameter) number that should be outputted
        If pop_or_param < 0 then it indicates the birth_quantity that should be outputted (TODO: change this API)
        the value inactive_value (as a string) is used in the CSV file for inactive cells.
        A header line containing the current time will be written to the file
        A header line of the form #xmin... will be written to the file
        additional_header_lines will be written verbatim (including any \n, etc) into the file"""

        self.outputCSVsubset(0, 0, self.grid.getNx(), self.grid.getNy(), filename, pop_or_param, inactive_value, additional_header_lines)

    cpdef outputCSVsubset(self, unsigned x1, unsigned y1, unsigned x2, unsigned y2, str filename, int pop_or_param, str inactive_value, str additional_header_lines):
        """Outputs cell information for given population or parameter number to filename in CSV format.
        Performs this for grid cell number (not coordinates) with x1 <= x < x2 and y1 <= y < y2.
        If pop_or_param >= 0 then it indicates the population (or parameter) number that should be outputted
        If self.num_diffusing_populations_at_cell <= pop_or_param < 0 then it indicates the birth_quantity that should be outputted (TODO: change this API)
        If 2*self.num_diffusing_populations_at_cell <= pop_or_param < self.num_diffusing_populations_at_cell then it indicates the yyp that should be outputted (TODO: change this API)
        the value inactive_value (as a string) is used in the CSV file for inactive cells.
        A header line containing the current time will be written to the file
        A header line of the form #xmin... will be written to the file
        additional_header_lines will be written verbatim (including any \n, etc) into the file"""

        if pop_or_param >= int(self.num_quantities_at_cell):
            raise ValueError("You requested pop_or_param number " + str(pop_or_param) + " but there are only " + str(self.num_quantities_at_cell) + " at each cell")
        elif pop_or_param < -2*int(self.num_diffusing_populations_at_cell):
            raise ValueError("You requested pop_or_param number " + str(pop_or_param) + " but the minimum accepted value is " + str(-2*int(self.num_diffusing_populations_at_cell)))
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
        
        if x2 > x_max or y2 > y_max:
            raise ValueError("Cell range (" + str(x1) + "," + str(y1) + ") to (" + str(x2) + "," + str(y2) + ") outside range of grid with x_max = " + str(x_max) + " and y_max = " + str(y_max))
        
        f = open(filename, "w")
        f.write("#File written at: " + time.asctime() + "\n")
        f.write("#xmin=" + str(self.grid.getXmin()) + ",ymin=" + str(self.grid.getYmin()) + ",cell_size=" + str(self.grid.getCellSize()) + ",nx=" + str(x_max) + ",ny=" + str(y_max) + "\n")
        f.write(additional_header_lines)
        for y_ind in range(y1, y2):
            for x_ind in range(x1, x2 - 1):
                ind = self.grid.internal_global_index(x_ind, y_ind)
                active_ind = active[ind]
                if active_ind == num_cells:
                    # inactive cell
                    f.write(inactive_value + ",")
                elif pop_or_param >= 0:
                    f.write(str(self.all_quantities[active_ind * self.num_quantities_at_cell + pop_or_param]) + ",")
                elif pop_or_param >= -1*int(self.num_diffusing_populations_at_cell):
                    f.write(str(self.birth_quantities[active_ind * self.num_diffusing_populations_at_cell - pop_or_param - 1]) + ",") # (B)
                else:
                    f.write(str(self.yyp_quantities[active_ind * self.num_diffusing_populations_at_cell - pop_or_param - self.num_diffusing_populations_at_cell - 1]) + ",") # (Y')
            x_ind = x2 - 1
            ind = self.grid.internal_global_index(x_ind, y_ind)
            active_ind = active[ind]
            if active_ind == num_cells:
                # inactive cell
                f.write(inactive_value + "\n")
            elif pop_or_param >= 0:
                f.write(str(self.all_quantities[active_ind * self.num_quantities_at_cell + pop_or_param]) + "\n")
            elif pop_or_param >= -1*int(self.num_diffusing_populations_at_cell):
                f.write(str(self.birth_quantities[active_ind * self.num_diffusing_populations_at_cell - pop_or_param - 1]) + "\n")
            else:
                f.write(str(self.yyp_quantities[active_ind * self.num_diffusing_populations_at_cell - pop_or_param - self.num_diffusing_populations_at_cell - 1]) + "\n")
        f.close()
