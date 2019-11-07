import array
cimport cpython.array as array

cdef class CellDynamicsBase:
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
    
    cpdef unsigned getNumberOfPopulations(self)
    """Returns number of populations (maleGw, femaleGG, larvaeMaleGw, whatever) in the Cell dynamics"""

    cpdef unsigned getNumberOfDiffusingPopulations(self)
    """Returns the number of populations that diffuse"""

    cpdef array.array getDiffusingIndices(self)
    """Returns: diffusing_indices[i] is the i^th populations that is diffusing
    0 <= diffusing_indices[i] < num_populations.  0 <= i < num_diffusing"""

    cpdef unsigned getNumberOfAdvectingPopulations(self)
    """Returns the number of populations that advect"""

    cpdef array.array getAdvectingIndices(self)
    """Returns: advecting_indices[i] is the i^th populations that is advecting
    0 <= advecting_indices[i] < num_populations.  0 <= i < num_advecting"""

    cpdef unsigned getNumberOfParameters(self)
    """Returns the number of parameters (carrying capacity, mortality rate, etc) in the cell dynamics"""

    cpdef void evolve(self, float timestep, float[:] pops_and_params)
    """Performs one timestep of evolution
    Note: pops_and_params is a pointer to an array of floats, so any changes made to its values will be evident to the calling function
    Note: pops_and_params will be of size num_populations + num_parameters, and should not be resized
    Note: the first num_populations entries of pops_and_params will be the population values, while the remainder are the parameter values"""


cdef class CellDynamicsStatic15_9_3_2(CellDynamicsBase):
    """No dynamics within the cell (all populations are static as far as the cell is concerned)
    15 populations, 9 of them are diffusing and 3 advecting, with 2 parameters"""
    

cdef class CellDynamicsLogistic1_1(CellDynamicsBase):
    """Logistic growth of a single diffusing and advecting population
    1 parameter, which is the carrying capacity"""
    
    # growth rate
    cdef float growth_rate
