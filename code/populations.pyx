import array
cimport cpython.array as array
from grid cimport Grid
from cell cimport Cell

cdef class Populations:

    def __init__(self, Grid grid):
        self.grid = grid
        self.num_active_cells = self.grid.getNumActiveCells()

        self.num_pop_per_cell = Cell().getNumberOfPopulations()

        self.num_pops = self.num_active_cells * self.num_pop_per_cell
        
        self.pop = array.array('f', [0.0] * self.num_active_cells * self.num_pop_per_cell)


    cdef array.array getPopulations(self):
        return self.pop

    cdef unsigned getIndexOfCell(self, unsigned active_cell_index):
        """returns the index in self.pop that corresponds to the 0th population at the given active_cell_index"""
        return active_cell_index * self.num_pop_per_cell

    cpdef setPopulation(self, unsigned active_cell_index, list pop):
        """Sets the populations at active_cell_index to the numbers given in the list 'pop'"""
        if len(pop) != self.num_pop_per_cell:
            raise ValueError("Length of pop in setPopulation is " + str(len(pop)) + " which should be " + str(self.num_pop_per_cell))
        if active_cell_index >= self.num_active_cells:
            raise ValueError("Active cell index is " + str(active_cell_index) + " which should be less than " + str(self.num_active_cells))
        cdef unsigned start_index = self.getIndexOfCell(active_cell_index)
        cdef unsigned i
        for i in range(self.num_pop_per_cell):
            self.pop.data.as_floats[start_index + i] = pop[i]

