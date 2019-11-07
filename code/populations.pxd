import array
cimport cpython.array as array
from grid cimport Grid

cdef class Populations:
    """Holds mosquito-population information at every active cell of the entire grid"""

    # the grid, which allows extraction of things such as the number of active cells
    cdef Grid grid

    # number of active cells
    cdef unsigned num_active_cells

    # number of populations (maleGG, femaleGw, larvae, etc) at each cell
    cdef unsigned num_pop_per_cell

    # total number of populations handled by this class ( = num_active_cells * num_pop_per_cell)
    cdef unsigned num_pops

    # the array of populations, which is of length num_pops
    cdef array.array pop

    cdef array.array getPopulations(self)
    """returns the pop array"""

    cdef unsigned getIndexOfCell(self, unsigned active_cell_index)
    """returns the index in self.pop that corresponds to the 0th population at the given active_cell_index"""

    cpdef setPopulation(self, unsigned active_cell_index, list pop)
    """Sets the populations at active_cell_index to the numbers given in the list 'pop'"""
