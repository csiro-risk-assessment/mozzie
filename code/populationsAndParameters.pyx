import array
cimport cpython.array as array
from grid cimport Grid
from cellDynamics cimport CellDynamicsBase

cdef class PopulationsAndParameters:

    def __init__(self, Grid grid, CellDynamicsBase cell):
        self.grid = grid
        self.cell = cell

        self.num_active_cells = self.grid.getNumActiveCells()

        self.num_pop_per_cell = self.cell.getNumberOfPopulations()

        self.num_param_per_cell = self.cell.getNumberOfParameters()

        self.num_quantities_per_cell = self.num_pop_per_cell + self.num_param_per_cell

        self.num_quantities = self.num_active_cells * self.num_quantities_per_cell
        
        self.quantities = array.clone(array.array('f', []), self.num_active_cells * self.num_quantities_per_cell, zero = True)
        
        self.quantity = array.clone(array.array('f', []), self.num_quantities_per_cell, zero = True)


    cpdef array.array getQuantities(self):
        return self.quantities

    cpdef CellDynamicsBase getCell(self):
        return self.cell

    cdef unsigned getIndexOfCell(self, unsigned active_cell_index):
        """returns the index in self.quantities that corresponds to the 0th population at the given active_cell_index"""
        return active_cell_index * self.num_quantities_per_cell

    cpdef setPopulationAndParameters(self, unsigned active_cell_index, list pop_and_params):
        """Sets the populations and parameters at active_cell_index to the numbers given in the list 'pop_and_params'
        This is a slower method than accessing self.quantities directly because there is a lot of bounds checking"""
        if len(pop_and_params) != self.num_quantities_per_cell:
            raise ValueError("Length of pop_and_params in setPopulationAndParameters is " + str(len(pop_and_params)) + " which should be " + str(self.num_quantities_per_cell))
        if active_cell_index >= self.num_active_cells:
            raise ValueError("Active cell index is " + str(active_cell_index) + " which should be less than " + str(self.num_active_cells))
        cdef unsigned start_index = self.getIndexOfCell(active_cell_index)
        cdef unsigned i
        for i in range(self.num_quantities_per_cell):
            self.quantities.data.as_floats[start_index + i] = pop_and_params[i]

    cpdef setPopulationAndParametersFromXY(self, float x, float y, list pop_and_params):
        """Sets the populations and parameters at given (x, y) to the numbers given in the list 'pop_and_params'
        This is a slower method than accessing self.quantities directly because there is a lot of bounds checking"""
        if x < self.grid.getXmin() or x > self.grid.getXmin() + self.grid.getCellSize() * self.grid.getNx() or y < self.grid.getYmin() or y > self.grid.getYmin() + self.grid.getCellSize() * self.grid.getNy():
            raise ValueError("x or y in setPopulationAndParametersFromXY is not inside the grid")
        cdef unsigned x_ind = int((x - self.grid.xmin) / self.grid.cell_size + 0.5)
        cdef unsigned y_ind = int((y - self.grid.ymin) / self.grid.cell_size + 0.5)
        cdef unsigned global_ind = self.grid.internal_global_index(x_ind, y_ind)
        cdef unsigned active_ind = self.grid.active_index[global_ind]
        self.setPopulationAndParameters(active_ind, pop_and_params)

    cpdef array.array getPopulationAndParameters(self, unsigned active_cell_index):
        """Gets the populations and parameters at active_cell_index and returns the list 'self.quantity'
        This is a slower method than accessing self.quantities directly because there is a lot of bounds checking"""
        if active_cell_index >= self.num_active_cells:
            raise ValueError("Active cell index is " + str(active_cell_index) + " which should be less than " + str(self.num_active_cells))
        cdef unsigned start_index = self.getIndexOfCell(active_cell_index)
        cdef unsigned i
        for i in range(self.num_quantities_per_cell):
            self.quantity[i] = self.quantities.data.as_floats[start_index + i]
        return self.quantity

    cpdef array.array getPopulationAndParametersFromXY(self, float x, float y):
        """Gets the populations and parameters at given (x, y) and returns the list 'self.quantity'
        This is a slower method than accessing self.quantities directly because there is a lot of bounds checking"""
        if x < self.grid.getXmin() or x > self.grid.getXmin() + self.grid.getCellSize() * self.grid.getNx() or y < self.grid.getYmin() or y > self.grid.getYmin() + self.grid.getCellSize() * self.grid.getNy():
            raise ValueError("x or y in setPopulationAndParametersFromXY is not inside the grid")
        cdef unsigned x_ind = int((x - self.grid.xmin) / self.grid.cell_size + 0.5)
        cdef unsigned y_ind = int((y - self.grid.ymin) / self.grid.cell_size + 0.5)
        cdef unsigned global_ind = self.grid.internal_global_index(x_ind, y_ind)
        cdef unsigned active_ind = self.grid.active_index[global_ind]
        return self.getPopulationAndParameters(active_ind)

    cpdef setOverActiveGrid(self, unsigned pop_and_param_number, array.array pop_or_param_array):
        """pop_or_param_array is a float array defined over the entire active grid.
        The population or parameter with number pop_and_param_number is set to this array"""
        if pop_and_param_number >= self.num_quantities_per_cell:
            raise ValueError("pop_and_param_number is " + str(pop_and_param_number) + " but this must be less than the number of quantities per cell, which is " + str(self.num_quantities_per_cell))
        if len(pop_or_param_array) != self.num_active_cells:
            raise ValueError("length of pop_or_param_array is " + str(len(pop_or_param_array)) + " which must be equal to the number of active cells, which is " + str(self.num_active_cells))
        cdef unsigned i
        for i in range(self.num_active_cells):
            self.quantities.data.as_floats[i * self.num_quantities_per_cell + pop_and_param_number] = pop_or_param_array.data.as_floats[i]

