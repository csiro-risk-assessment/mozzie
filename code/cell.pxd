import array
cimport cpython.array as array
cdef class Cell:
    """Holds and manipulates information at a single cell"""

    # the population array at the node
    cpdef array.array population
    
    cpdef unsigned getNumberOfPopulations(self)

    cpdef array.array getPopulation(self)

    cpdef setPopulations(self, list pops)

    cpdef unsigned getNumberOfDiffusingPopulations(self)

    cpdef float getDiffusingSingletonNoCheck(self, unsigned diffusing_index)

    cpdef array.array getDiffusingPopulation(self)

    cpdef void addDiffusingSingletonNoCheck(self, unsigned diffusing_index, float add_this)

    cdef void addDiffusingPopulationNoCheck(self, array.array add_this)

    cpdef addDiffusingPopulation(self, list add_this)
    """Adds the floats specified in the list add_this to the diffusing populations"""

    cpdef unsigned getNumberOfAdvectingPopulations(self)

    cpdef float getAdvectingSingletonNoCheck(self, unsigned advecting_index)

    cpdef array.array getAdvectingPopulation(self)

    cpdef void addAdvectingSingletonNoCheck(self, unsigned advecting_index, float add_this)

    cdef void addAdvectingPopulationNoCheck(self, array.array add_this)

    cpdef addAdvectingPopulation(self, list add_this)
    """Adds the floats specified in the list add_this to the advecting populations"""
