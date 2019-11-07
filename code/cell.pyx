import array
cimport cpython.array as array

# number of populations (maleGw, femaleGG, larvaeMaleGw, whatever) that this Cell looks after
cpdef unsigned num_populations = 12

# number of populations that diffuse
cpdef unsigned num_diffusing = 9

# self.population[diffusing_indices[i]] is the i^th population that is diffusing
# In the following 'I' means unsigned ints
cpdef array.array diffusing_indices = array.array('I', [0, 2, 3, 4, 5, 6, 8, 10, 11])

# number of populations that advect
cpdef unsigned num_advecting = 3

# self.population[advecting_indices[i]] is the i^th population that is advecting
# In the following 'I' means unsigned ints
cpdef array.array advecting_indices = array.array('I', [0, 4, 10])

cdef class Cell:
    """Holds and manipulates information at a single cell"""

    def __init__(self):
        """Initialise the Cell with a zero population"""
        self.population = array.clone(array.array('f', []), num_populations, zero = True)

    cpdef unsigned getNumberOfPopulations(self):
        return num_populations

    cpdef setPopulations(self, list pops):
        """Set the populations in the cell"""
        if len(pops) != num_populations:
            raise ValueError("Number of population values must be " + str(num_populations))
        self.population = array.array('f', pops)
    
    cpdef array.array getPopulation(self):
        return self.population

    cpdef unsigned getNumberOfDiffusingPopulations(self):
        return num_diffusing

    cpdef array.array getDiffusingIndices(self):
        return diffusing_indices

    cpdef float getDiffusingSingletonNoCheck(self, unsigned diffusing_index):
        return self.population.data.as_floats[diffusing_indices.data.as_uints[diffusing_index]]

    cpdef array.array getDiffusingPopulation(self):
        return array.array('f', [self.population[ind] for ind in diffusing_indices])

    cpdef void addDiffusingSingletonNoCheck(self, unsigned diffusing_index, float add_this):
        self.population.data.as_floats[diffusing_indices.data.as_uints[diffusing_index]] += add_this

    cdef void addDiffusingPopulationNoCheck(self, array.array add_this):
        for ind in range(num_diffusing):
            self.addDiffusingSingletonNoCheck(ind, add_this.data.as_floats[ind])

    cpdef addDiffusingPopulation(self, list add_this):
        """Adds the floats specified in the list add_this to the diffusing populations"""
        if len(add_this) != num_diffusing:
            raise ValueError("Attempted to add a diffusing population of length " + str(len(add_this)) + " but should have length " + str(num_diffusing))
        self.addDiffusingPopulationNoCheck(array.array('f', add_this))

        

    cpdef unsigned getNumberOfAdvectingPopulations(self):
        return num_advecting

    cpdef float getAdvectingSingletonNoCheck(self, unsigned advecting_index):
        return self.population.data.as_floats[advecting_indices.data.as_uints[advecting_index]]

    cpdef array.array getAdvectingPopulation(self):
        return array.array('f', [self.population[ind] for ind in advecting_indices])

    cpdef void addAdvectingSingletonNoCheck(self, unsigned advecting_index, float add_this):
        self.population.data.as_floats[advecting_indices.data.as_uints[advecting_index]] += add_this

    cdef void addAdvectingPopulationNoCheck(self, array.array add_this):
        for ind in range(num_advecting):
            self.addAdvectingSingletonNoCheck(ind, add_this.data.as_floats[ind])

    cpdef addAdvectingPopulation(self, list add_this):
        """Adds the floats specified in the list add_this to the advecting populations"""
        if len(add_this) != num_advecting:
            raise ValueError("Attempted to add an advecting population of length " + str(len(add_this)) + " but should have length " + str(num_advecting))
        self.addAdvectingPopulationNoCheck(array.array('f', add_this))
