import array
cimport cpython.array as array
### Following doesn't work with python2.7, yet
###cimport numpy as np
###import numpy as np

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

cdef class CellDynamicsBeeton2_2(CellDynamicsBase):
    """Dynamics of a hybridising 2-subspecies model described in
    Beeton, Hosack, Wilkins, Forbes, Ickowicz and Hayes, Journal of Theoretical Biology 2019.
    There are 2 parameters, which are the carrying capacities Kx and Ky.
    The dynamics is defined in Eqns(4.1) and (4.2) in Beeton et al."""
    cpdef float mux # mu_x
    cpdef float muy # mu_y
    cpdef float gax # gamma_x
    cpdef float gay # gamma_y
    cpdef float axy # alpha_{xy}
    cpdef float ayx # alpha_{yx}
    cpdef float w   # w

    cpdef void setMuX(self, float mux)
    """Sets mu_x"""

    cpdef void setMuY(self, float muy)
    """Sets mu_y"""

    cpdef void setGaX(self, float gax)
    """Sets gamma_x"""

    cpdef void setGaY(self, float gay)
    """Sets gamma_y"""

    cpdef void setAxy(self, float axy)
    """Sets alpha_{xy}"""

    cpdef void setAyx(self, float ayx)
    """Sets alpha_{yx}"""

    cpdef void setW(self, float w)
    """Sets w"""

    cpdef float getMuX(self)
    """Gets mu_x"""

    cpdef float getMuY(self)
    """Gets mu_y"""

    cpdef float getGaX(self)
    """Gets gamma_x"""

    cpdef float getGaY(self)
    """Gets gamma_y"""

    cpdef float getAxy(self)
    """Gets alpha_{xy}"""

    cpdef float getAyx(self)
    """Gets alpha_{yx}"""

    cpdef float getW(self)
    """Gets w"""

cdef class CellDynamicsMosquito23(CellDynamicsBase):
    """Solves Mosquito ODE with 2 sexes and 3 genotypes"""

    # male, female, always
    cdef unsigned num_sexes

    # ww, Gw, GG, always
    cdef unsigned num_genotypes

    # iheritance_cube[x,y,z] = probability of mother genotype x, father genotype y producing offspring genotype z
    # where index 0, 1, 2 = ww, Gw, GG respectively
    ###cpdef np.ndarray inheritance_cube

    # doco
    cdef float mu_larvae

    # codo
    cdef float mu_adult

    # doco
    cdef float fecundity

    # doco
    cdef float aging_rate

    # age categories are: larvae0, larvae1, larvae2, ..., larvaeN, adults
    cdef unsigned num_ages

    # number of species
    cdef unsigned num_species

    # accuracy of PMB
    cdef float accuracy

    # doco
    ###cdef np.ndarray fecundity_proportion

    # doco
    ###cdef np.ndarray ipm

    # doco
    ###cdef np.ndarray ipf

    # the population numbers, as a numpy array
    ###cdef np.ndarray xx

    # the carrying capacity
    cdef float kk

    cdef void setInternalParameters(self, unsigned num_ages, unsigned num_species, float accuracy)
    """Given num_ages, num_species and accuracy, set all internal parameters (num_populations, diffusing_indices, fecundity_proportion, etc) that depend on these"""

    cpdef void setMuLarvae(self, float mu_larvae)
    """Set mu_larvae"""

    cpdef float getMuLarvae(self)
    """Get mu_larvae"""
    
    cpdef void setMuAdult(self, float mu_adult)
    """Set mu_adult"""

    cpdef float getMuAdult(self)
    """Get mu_adult"""
    
    cpdef void setFecundity(self, float fecundity)
    """Set fecundity"""

    cpdef float getFecundity(self)
    """Get fecundity"""
    
    cpdef void setAgingRate(self, float aging_rate)
    """Set aging_rate"""

    cpdef float getAgingRate(self)
    """Get aging_rate"""
    
    cpdef void setNumAges(self, unsigned num_ages)
    """Set number of ages (larvae0, larvae1, ...., adults.  This calls setInternalParameters to appropriately set other parameters, given num_ages"""

    cpdef unsigned getNumAges(self)
    """Get number of ages"""
    
    cpdef void setNumSpecies(self, unsigned num_species)
    """Set number of species.  This calls setInternalParameters to appropriately set other parameters, given num_species"""

    cpdef unsigned getNumSpecies(self)
    """Get number of species"""
    
    cpdef void setAccuracy(self, float accuracy)
    """Set accuracy of PMB.  This calls setInternalParameters to appropriately set other parameters, given accuracy"""

    cpdef float getAccuracy(self)
    """Get accuracy"""
    
