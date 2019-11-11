import array
cimport cpython.array as array
import numpy as np
from scipy.integrate import solve_ivp

cdef class CellDynamicsBase:

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

    cpdef void evolve(self, float timestep, float[:] pops_and_params):
        return


cdef class CellDynamicsStatic15_9_3_2(CellDynamicsBase):
    """No dynamics within the cell (all populations are static as far as the cell is concerned)
    15 populations, 9 of them are diffusing and 3 advecting, with 2 parameters"""
    def __init__(self):
        super().__init__()
        self.num_populations = 15
        self.num_diffusing = 9
        self.diffusing_indices = array.array('I', [0, 2, 3, 4, 5, 6, 8, 10, 11])
        self.num_advecting = 3
        self.advecting_indices = array.array('I', [0, 4, 10])
        self.num_parameters = 2

    cpdef void evolve(self, float timestep, float[:] pops_and_params):
        """No dynamics here"""
        return
    

cdef class CellDynamicsLogistic1_1(CellDynamicsBase):
    """Logistic growth of a single diffusing and advecting population
    1 parameter, which is the carrying capacity"""

    def __init__(self):
        super().__init__()
        self.num_populations = 1
        self.num_diffusing = 1
        self.diffusing_indices = array.array('I', [0])
        self.num_advecting = 1
        self.advecting_indices = array.array('I', [0])
        self.num_parameters = 1
        self.growth_rate = 0.01 # fixed for this example

    cpdef void evolve(self, float timestep, float[:] pops_and_params):
        # explicit time-stepping of logistic growth
        # as defined in the doco for "evolve",
        # pops_and_params[0] is the current population
        # pops_and_params[1] is the carrying capacity
        pops_and_params[0] = pops_and_params[0] + timestep * self.growth_rate * pops_and_params[0] * (1.0 - pops_and_params[0] / pops_and_params[1])


cdef class CellDynamicsBeeton2_2(CellDynamicsBase):
    def __init__(self):
        super().__init__()
        self.num_populations = 2
        self.num_diffusing = 2 # both subspecies diffuse
        self.diffusing_indices = array.array('I', [0, 1])
        self.num_advecting = 2 # both subspecies advect
        self.advecting_indices = array.array('I', [0, 1])
        self.num_parameters = 2 # carrying capacity Kx, and carrying capacity Ky
        # Default parameters in Equations (4.1) and (4.2) of Beeton et al, are as set in Fig5c
        self.mux = 0.7
        self.muy = 0.8
        self.gax = 1.0
        self.gay = 1.0
        self.axy = 0.4
        self.ayx = 0.4
        self.w = 0.05

    cpdef void setMuX(self, float mux):
        self.mux = mux

    cpdef void setMuY(self, float muy):
        self.muy = muy

    cpdef void setGaX(self, float gax):
        self.gax = gax

    cpdef void setGaY(self, float gay):
        self.gay = gay

    cpdef void setAxy(self, float axy):
        self.axy = axy

    cpdef void setAyx(self, float ayx):
        self.ayx = ayx

    cpdef void setW(self, float w):
        self.w = w

    cpdef float getMuX(self):
        return self.mux

    cpdef float getMuY(self):
        return self.muy

    cpdef float getGaX(self):
        return self.gax

    cpdef float getGaY(self):
        return self.gay

    cpdef float getAxy(self):
        return self.axy

    cpdef float getAyx(self):
        return self.ayx

    cpdef float getW(self):
        return self.w


    cpdef void evolve(self, float timestep, float[:] pops_and_params):
        # explicit time-stepping of Beeton et al.'s equations
        # as defined in the doco for "evolve",
        # pops_and_params[0] is the current x
        # pops_and_params[1] is the current y
        # pops_and_params[2] is the carrying capacity Kx
        # pops_and_params[3] is the carrying capacity Ky
        cpdef float x = pops_and_params[0]
        cpdef float y = pops_and_params[1]
        cpdef float kx = pops_and_params[2]
        cpdef float ky = pops_and_params[3]
        cpdef float f = - self.mux + (1.0 - (x + self.axy * y) / kx) * (x / (x + self.w * y)) * self.gax
        cpdef float g = - self.muy + (1.0 - (self.ayx * x + y) / ky) * (self.gay + self.w * x / (x + self.w * y) * self.gax)
        # max means (x,y) doesn't stray outside physical quadrant, even if timestep is too large
        pops_and_params[0] = max(0.0, x + timestep * f * x)
        pops_and_params[1] = max(0.0, y + timestep * g * y)


cdef class CellDynamicsMosquito23(CellDynamicsBase):
    def __init__(self):
        super().__init__()

        # Following quantities are always fixed
        self.num_sexes = 2 # male, female, always
        self.num_genotypes = 3 # ww, Gw, GG, always
        self.num_parameters = 1 # the only spatially-varying quantity is the carrying capacity, K
        # iheritance_cube[x,y,z] = probability of mother genotype x, father genotype y producing offspring genotype z
        # where index 0, 1, 2 = ww, Gw, GG respectively
        self.inheritance_cube = np.array([[[1., 0., 0.],# ww x ww
                                           [ .5, .5, 0. ],
                                           [ 0., 1., 0. ]],
                                          [[ .5, .5, 0. ],
                                           [ .25, .5, .25 ], # Gw x Gw
                                           [ 0., .25, .25 ]],
                                          [[ 0., 1., 0. ],
                                           [ 0., .5, .5 ],
                                           [ 0., 0., 1. ]]]) # GG x GG

        # Following parameters may be set by user.
        self.mu_larvae = 0.1
        self.mu_adult = 0.1
        self.fecundity = 0.9
        self.aging_rate = 0.1

        self.setInternalParameters(2, 1, 0.95) # default to num_ages = 2, num_species = 1, accuracy = 0.95

    cdef void setInternalParameters(self, unsigned num_ages, unsigned num_species, float accuracy):
        self.num_ages = num_ages # age categories are: larvae0, larvae1, larvae2, ..., larvaeN, adults
        self.num_species = num_species
        self.accuracy = accuracy
        
        self.num_populations = self.num_ages * self.num_sexes * self.num_genotypes * self.num_species

        self.num_diffusing = 6 # only age = num_ages - 1 (that is, adults) diffuses
        cdef unsigned offset = (self.num_ages - 1) * self.num_sexes * self.num_genotypes * self.num_species
        moving_indices = [offset + i for i in range(6)]
        self.diffusing_indices = array.array('I', moving_indices) # the age = num_ages = 1 mosquitoes
        self.num_advecting = 6 # only age = num_ages - 1 (that is, adults) diffuses
        self.advecting_indices = array.array('I', moving_indices) # the age = num_ages = 1 mosquitoes

        # calculate fecundity allocated to male and female larvae
        # fecundity_proportion[z,x,y] = proportion of total fecundity of mother genotype x, father genotype y to producing offspring sex z
        # where index of z: 0, 1 = male, female
        # note - proportions will not add to 1 for male bias
        cdef float fprop = 1.0 / self.accuracy - 1.0
        prow = np.array([0.5] * self.num_genotypes)
        self.fecundity_proportion = np.array([[prow, prow, prow], [prow, fprop * prow, fprop * prow]]) # has length self.num_sexes

        # pre-calculate inheritance_cube * fecundity_proportion for males and females
        self.ipm = np.zeros((self.num_genotypes, self.num_genotypes, self.num_genotypes))
        self.ipf = np.zeros((self.num_genotypes, self.num_genotypes, self.num_genotypes))

        for j in range(self.num_genotypes):
            self.ipm[:,:,j] = self.inheritance_cube[:,:,j] * self.fecundity_proportion[0,:,:]
            self.ipf[:,:,j] = self.inheritance_cube[:,:,j] * self.fecundity_proportion[1,:,:]

        # transpose these for use in later matrix multiplication
        # (useful form to have in case of implicit)
        for j in range(self.num_genotypes):
            self.ipm[j,:,:] = np.transpose(self.ipm[j,:,:])
            self.ipf[j,:,:] = np.transpose(self.ipf[j,:,:])

        # Following parameters are for use in evolve: most efficient to allocate them now
        self.xx = np.ones(self.num_populations)
        self.kk = 1.0
        

    cpdef void setMuLarvae(self, float mu_larvae):
        self.mu_larvae = mu_larvae

    cpdef float getMuLarvae(self):
        return self.mu_larvae
    
    cpdef void setMuAdult(self, float mu_adult):
        self.mu_adult = mu_adult

    cpdef float getMuAdult(self):
        return self.mu_adult
    
    cpdef void setFecundity(self, float fecundity):
        self.fecundity = fecundity

    cpdef float getFecundity(self):
        return self.fecundity
    
    cpdef void setAgingRate(self, float aging_rate):
        self.aging_rate = aging_rate

    cpdef float getAgingRate(self):
        return self.aging_rate
    
    cpdef void setNumAges(self, unsigned num_ages):
        self.setInternalParameters(num_ages, self.num_species, self.accuracy)

    cpdef unsigned getNumAges(self):
        return self.num_ages
    
    cpdef void setNumSpecies(self, unsigned num_species):
        self.setInternalParameters(self.num_ages, num_species, self.accuracy)

    cpdef unsigned getNumSpecies(self):
        return self.num_species
    
    cpdef void setAccuracy(self, float accuracy):
        self.setInternalParameters(self.num_ages, self.num_species, accuracy)

    cpdef float getAccuracy(self):
        return self.accuracy

    def fun(self, t, y):
        """Evaluates d(populations)/dt"""
        Y = np.reshape(y, (self.num_ages, self.num_sexes, self.num_genotypes, self.num_species))
        n = Y[:-1,:,:,:].sum() # total number of larvae (all but the last age class)
        ratio = Y[-1,0,:,:] # ratio of adult males (last age class) of given genotype and species
        ratio = ratio / ratio.sum()
        mat = np.zeros((self.num_populations, self.num_populations))
        # for each father's genotype and species, add the relevant proportion of fecundity for 
        # mothers of each genotype and that species, producing proportions of offspring 
        # of the relevant genotype (and the same species)
        for i in range(self.genotypes): # fathers' genotypes
            for j in range(self.num_species): # species
                mat[j:(self.num_genotypes * self.num_species):self.num_species, (self.num_genotypes * self.num_genotypes * (self.num_ages - 1) * self.num_species + j):self.num_populations:self.num_species] += self.ipm[i,:,:] * ratio[i,j]
                mat[(self.num_genotypes * self.num_species + j):(self.num_genotypes * self.num_sexes * self.nu_species):self.num_species, (self.genotypes * self.num_genotypes * (self.num_ages - 1) * self.num_species + j): self.num_populations: self.num_species] += self.ipf[i,:,:] * ratio[i,j]
        mat *= (1 - n / self.kk) * self.fecundity # scaling by fecundity and density dependence
        # mortality
        mat[range(self.num_genotypes * self.num_sexes * (self.num_ages - 1) * self.num_species), range(self.num_genotypes * self.num_sexes * (self.num_ages - 1) * self.num_species)] = - self.mu_larvae - self.aging_rate
        mat[range(self.num_genotypes * self.num_sexes * (self.num_ages - 1) * self.num_species, self.num_populations), range(self.num_genotypes * self.num_sexes * (self.num_ages - 1) * self.num_species, self.num_populations)] = - self.mu_adult
        # aging
        for i in range(self.num_ages - 1):
            mat[range(self.num_genotypes * self.num_sexes * (i + 1) * self.num_species, self.num_genotypes * self.num_sexes * (i + 2) * self.num_species), range(self.num_genotypes * self.num_sexes * i * self.num_species, self.num_genotypes * self.num_sexes * (i + 1) * self.num_species)] = self.aging_rate
        dXdt = mat.dot(y)
        return dXdt

    cpdef void evolve(self, float timestep, float[:] pops_and_params):
        # copy into self.xx for use in self.fun
        cdef unsigned ind
        for ind in range(self.num_populations):
            self.xx[ind] = pops_and_params[ind]
        self.kk = pops_and_params[self.num_populations]
        # solve
        sol = solve_ivp(self.fun, [0.0, timestep], self.xx)
        # copoy back
        for ind in range(self.num_populations):
            pops_and_params[ind] = sol.y[ind, -1]
        
        

    
