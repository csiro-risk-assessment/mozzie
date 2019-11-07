import array
cimport cpython.array as array

cdef class CellDynamics:
    """Manipulates information at a single cell, in particular this class solves lifecycle ODEs"""

    # number of populations (maleGw, femaleGG, larvaeMaleGw, whatever) in the Cell dynamics
    cpdef unsigned num_populations

    # number of populations that diffuse
    cpdef unsigned num_diffusing

    # diffusing_indices[i] is the i^th population that is diffusing
    # 0 <= diffusing_indices[i] < num_populations.  0 <= i < num_diffusing
    cpdef array.array diffusing_indices

    # number of populations that advect
    cpdef unsigned num_advecting

    # advecting_indices[i] is the i^th population that is advecting
    # 0 <= advecting_indices[i] < num_populations.  0 <= i < num_advecting
    cpdef array.array advecting_indices

    # number of parameters (carrying capacity, mortality rate, etc) in the Cell dynamics
    cpdef unsigned num_parameters
    
    def __init__(self):
        """Initialise the Cell with zeroes"""
        self.num_populations = 0
        self.num_diffusing = 0
        self.diffusing_indices = array.array('I', [])
        self.num_advecting = 0
        self.advecting_indices = array.array('I', [])
        self.num_parameters = 0

    cpdef unsigned getNumberOfPopulations(self):
        return self.num_populations

    cpdef unsigned getNumberOfDiffusingPopulations(self):
        return self.num_diffusing

    cpdef array.array getDiffusingIndices(self):
        return self.diffusing_indices

    cpdef unsigned getNumberOfAdvectingPopulations(self):
        return self.num_advecting

    cpdef array.array getAdvectingIndices(self):
        return self.advecting_indices

    cpdef unsigned getNumberOfParameters(self):
        return self.num_parameters

    

