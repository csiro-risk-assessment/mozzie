import sys
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

        self.num_genotypes2 = self.num_genotypes * self.num_genotypes

        # Following parameters may be set by user.
        self.mu_larvae = 0.1
        self.mu_adult = 0.1
        self.fecundity = 0.9
        self.aging_rate = 0.1

        # Set default carrying capacity
        self.kk = 1.0

        self.inheritance_cube = [[[1., 0., 0.],# ww x ww
                                  [ .5, .5, 0. ], # ww x Gw
                                  [ 0., 1., 0. ]], # ww x GG
                                 [[ .5, .5, 0. ], # Gw x ww
                                  [ .25, .5, .25 ], # Gw x Gw
                                  [ 0., .5, .5 ]], # Gw x GG
                                 [[ 0., 1., 0. ], # GG x ww
                                  [ 0., .5, .5 ], # GG x Gw
                                  [ 0., 0., 1. ]]] # GG x GG
        
        # size ipm_array and ipf_array correctly
        self.ipm_array = array.clone(array.array('f', []), self.num_genotypes * self.num_genotypes * self.num_genotypes, zero = False)
        self.ipf_array = array.clone(array.array('f', []), self.num_genotypes * self.num_genotypes * self.num_genotypes, zero = False)
        
        self.setInternalParameters(2, 1, 0.95) # default to num_ages = 2, num_species = 1, accuracy = 0.95

    cdef void setInternalParameters(self, unsigned num_ages, unsigned num_species, float accuracy):
        self.num_ages = num_ages # age categories are: larvae0, larvae1, larvae2, ..., larvaeN, adults
        self.num_species = num_species
        self.accuracy = accuracy
        
        self.num_populations = self.num_ages * self.num_sexes * self.num_genotypes * self.num_species

        # loop counters
        cdef unsigned gt0, gt1, gt2

        # set diffusing and advecting information
        self.num_diffusing = self.num_sexes * self.num_genotypes * self.num_species # only adults diffuse
        self.num_advecting = self.num_sexes * self.num_genotypes * self.num_species # only adults advect
        self.diffusing_indices = array.clone(array.array('I', []), self.num_diffusing, zero = False)
        self.advecting_indices = array.clone(array.array('I', []), self.num_advecting, zero = False)
        cdef unsigned sex, genotype, species, offset
        cdef unsigned ind = 0
        cdef unsigned age = self.num_ages - 1
        for sex in range(self.num_sexes):
            for genotype in range(self.num_genotypes):
                for species in range(self.num_species):
                    offset = species + genotype * self.num_species + sex * self.num_species * self.num_genotypes + age * self.num_species * self.num_genotypes * self.num_sexes
                    self.diffusing_indices.data.as_uints[ind] = offset
                    self.advecting_indices.data.as_uints[ind] = offset
                    ind = ind + 1

        # pre-calculate inheritance_cube * fecundity_proportion for males and females
        ipm = [[[0.0] * self.num_genotypes] * self.num_genotypes] * self.num_genotypes
        ipf = [[[0.0] * self.num_genotypes] * self.num_genotypes] * self.num_genotypes

        for gt2 in range(self.num_genotypes):
            for gt0 in range(self.num_genotypes):
                for gt1 in range(self.num_genotypes):
                    ipm[gt0][gt1][gt2] = self.inheritance_cube[gt0][gt1][gt2] * self.fecundity_proportion(0, gt0, gt1)
                    ipf[gt0][gt1][gt2] = self.inheritance_cube[gt0][gt1][gt2] * self.fecundity_proportion(1, gt0, gt1)

        # transpose these for use in later matrix multiplication
        # (useful form to have in case of implicit)
        # and put into cython array form
        for gt0 in range(self.num_genotypes):
            for gt1 in range(self.num_genotypes):
                for gt2 in range(self.num_genotypes):
                    self.setIPM(gt0, gt1, gt2, ipm[gt0][gt2][gt1])
                    self.setIPF(gt0, gt1, gt2, ipf[gt0][gt2][gt1])


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

        # useful indices
        cdef unsigned i, j, sp, gt0, gt1, ind0, ind1, ind2

        # Following lines can probably be optimised: currently lots of big copy-constructors
        Y = np.reshape(y, (self.num_ages, self.num_sexes, self.num_genotypes, self.num_species))
        if self.num_ages > 1:
			n = Y[:-1,:,:,:].sum() # total number of larvae (all but the last age class)
		else:
			n = Y.sum() # total number of mosquitoes (all assumed to be "adults")
        ratio = Y[-1,0,:,:] # ratio of adult males (last age class) of given genotype and species
        ratio = ratio / ratio.sum()

        # Nick wrote the following code blocks in vectorised form.  Andy has taken out vectorisation with a view to making "mat" and "ratio" cython arrays instead of a numpy objects
        
        mat = np.zeros((self.num_populations, self.num_populations))
        # for each father's genotype and species, add the relevant proportion of fecundity for 
        # mothers of each genotype and that species, producing proportions of offspring 
        # of the relevant genotype (and the same species)
        cdef unsigned offset_to_adult = (self.num_ages - 1) * self.num_species * self.num_genotypes * self.num_sexes
        for i in range(self.num_genotypes): # fathers' genotypes
            for sp in range(self.num_species): # species
                for gt0 in range(self.num_genotypes):
                    ind0 = sp + gt0 * self.num_species  # newborn (age=0) male (sex=0) of genotype gt0 and species sp
                    ind1 = sp + gt0 * self.num_species + self.num_species * self.num_genotypes # newborn (age=0) female (sex=1) of genotype gt0 and species sp
                    for gt1 in range(self.num_genotypes):
                        ind2 = offset_to_adult + sp + gt1 * self.num_species + self.num_species * self.num_genotypes # adult (age=num_ages-1) female (sex=0) of genotype gt1 and species sp
                        mat[ind0, ind2] += self.IPM(i, gt0, gt1) * ratio[i, sp]
                        mat[ind1, ind2] += self.IPF(i, gt0, gt1) * ratio[i, sp]
        mat *= (1 - n / self.kk) * self.fecundity # scaling by fecundity and density dependence

        # mortality, and aging into next age bracket
        cdef unsigned last_larvae =  (self.num_ages - 1) * self.num_species * self.num_genotypes * self.num_sexes
        if self.num_ages > 1:
			for i in range(last_larvae):
				mat[i][i] -= (self.mu_larvae + self.aging_rate)
        for i in range(last_larvae, self.num_populations):
            mat[i][i] -= self.mu_adult

        # aging from previous age bracket
        if self.num_ages > 1:
			cdef unsigned num_per_age =  self.num_species * self.num_genotypes * self.num_sexes
			cdef unsigned from_population
			for i in range(self.num_ages - 1):
				for j in range(num_per_age):
					from_population = num_per_age * i + j
					mat[from_population + num_per_age][from_population] = self.aging_rate

        dXdt = mat.dot(y)
        return dXdt

    cpdef void evolve(self, float timestep, float[:] pops_and_params):
        cdef unsigned ind
        self.kk = pops_and_params[self.num_populations]

        if self.kk <= 0.0:
            # instantly kill all populations
            for ind in range(self.num_populations):
                pops_and_params[ind] = 0.0
            return

        # copy into numpy array xx for use in self.fun
        xx = np.ones(self.num_populations)
        for ind in range(self.num_populations):
            xx[ind] = pops_and_params[ind]
        # solve
        sol = solve_ivp(self.fun, [0.0, timestep], xx)
        # copy back
        for ind in range(self.num_populations):
            pops_and_params[ind] = sol.y[ind, -1]
