import array
cimport cpython.array as array
from grid cimport Grid
from cellDynamics cimport CellDynamicsBase

cdef class PopulationsAndParameters:
    """Holds mosquito-population information at every active cell of the entire grid"""

    # the grid, which allows extraction of things such as the number of active cells
    cdef Grid grid

    # the CellDynamics object, which allows extraction of things like number of populations per cell
    cdef CellDynamicsBase cell

    # number of active cells
    cdef unsigned num_active_cells

    # number of populations (maleGG, femaleGw, larvae, etc) at each cell
    cdef unsigned num_pop_per_cell

    # number of parameters (carrying capacity, mortality rate, etc) at each cell
    cdef unsigned num_param_per_cell

    # number of quantities (populations and parameters) at each cell
    cdef unsigned num_quantities_per_cell

    # total number of quantites handled by this class ( = num_active_cells * num_quantities_per_cell)
    cdef unsigned num_quantities

    # the array of quantities, which is of length num_quantities
    cdef array.array quantities
    
    # the array of quantity (a temporary array to store quantities from a single cell), which is of length num_quantities_per_cell
    cdef array.array quantity

    cpdef array.array getQuantities(self)
    """returns the quantities array"""

    cpdef CellDynamicsBase getCell(self)
    """returns a reference to the cell used in the constructor"""

    cdef unsigned getIndexOfCell(self, unsigned active_cell_index)
    """returns the index in self.pop that corresponds to the 0th population at the given active_cell_index"""

    cpdef setPopulationAndParameters(self, unsigned active_cell_index, list pop)
    """Sets the populations and parameters at active_cell_index to the numbers given in the list 'pop_and_params'
    This is a slower method than accessing self.quantities directly because there is a lot of bounds checking"""

    cpdef setPopulationAndParametersFromXY(self, float x, float y, list pop_and_params)
    """Sets the populations and parameters at given (x, y) to the numbers given in the list 'pop_and_params'
    This is a slower method than accessing self.quantities directly because there is a lot of bounds checking"""

    cpdef array.array getPopulationAndParameters(self, unsigned active_cell_index)
    """Gets the populations and parameters at active_cell_index and returns the list 'self.quantity'
    This is a slower method than accessing self.quantities directly because there is a lot of bounds checking"""

    cpdef array.array getPopulationAndParametersFromXY(self, float x, float y)
    """Gets the populations and parameters at given (x, y) and returns the list 'self.quantity'
    This is a slower method than accessing self.quantities directly because there is a lot of bounds checking"""

    cpdef setOverActiveGrid(self, unsigned pop_and_param_number, array.array pop_or_param_array)
    """pop_or_param_array is a float array defined over the entire active grid.
    The population or parameter with number pop_and_param_number is set to this array"""
