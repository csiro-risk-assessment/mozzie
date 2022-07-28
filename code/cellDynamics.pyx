import sys
import array
cimport cpython.array as array
import numpy as np
from scipy.integrate import solve_ivp
from math import ceil, log, cos, sqrt
from libc.stdlib cimport rand, RAND_MAX
from libc.math cimport exp

cdef int binomial(int N, float p):
    cdef int count, wait
    cdef float tmp
    
    if ((N == 0) | (p == 0.)): # trivial cases
        return 0
    if ((N*p > 9.) & (N*(1-p) > 9.)): # normal approximation
        tmp = sqrt(-2*log(rand() / RAND_MAX)) * cos(2.*3.1415926535*rand()/RAND_MAX)
        count = int(tmp*sqrt(N*p*(1.-p)) + N*p + 0.5)
    else: # exact solution
        count = -1
        wait = 0
        if (p < 0.0001): # if p small approximate log(1-p) with -p
            p = -p
        else:
            p = log(1. - p)

        while (wait <= N):
            count += 1
            N -= wait
            wait = ceil( log(rand() / RAND_MAX) / p );
    return count
    
cdef int poisson(float l):
    cdef int k
    cdef float p, L
    
    if (l == 0.): # trivial case
        return 0
    elif (l > 10.): # normal approximation
        p = sqrt(-2*log(rand() / RAND_MAX)) * cos(2.*3.1415926535*rand()/RAND_MAX)
        k = int(p*sqrt(l) + l + 0.5)
        return k
    else: # exact solution
        L = exp(-l)
        k = 0
        p = 1.
        while (p > L):
            k += 1
            p *= rand() / RAND_MAX
        return k - 1

cdef class CellDynamicsBase:

    def __init__(self):
        """Initialise the Cell with zeroes"""
        
        self.num_populations = 0
        self.num_diffusing = 0
        self.diffusing_indices = array.array('I', [])
        self.num_advecting = 0
        self.advecting_indices = array.array('I', [])
        self.advection_class = array.array('I', [])
        self.num_parameters = 0
        self.small_value = 0.0

    cdef unsigned getNumberOfPopulations(self):
        return self.num_populations

    cdef unsigned getNumberOfDiffusingPopulations(self):
        return self.num_diffusing

    cdef array.array getDiffusingIndices(self):
        return self.diffusing_indices

    cdef unsigned getNumberOfAdvectingPopulations(self):
        return self.num_advecting

    cdef array.array getAdvectingIndices(self):
        return self.advecting_indices

    cdef array.array getAdvectionClass(self):
        return self.advection_class

    cdef setAdvectionClass(self, unsigned population_index, unsigned new_class):
        cdef int found = 0
        for ai in range(self.num_advecting):
            if self.advecting_indices[ai] == population_index:
                self.advection_class[ai] = new_class
                return
        raise ValueError("setAdvectionClass: population index " + str(population_index) + " has not defined to be advecting")

    cdef unsigned getNumberOfParameters(self):
        return self.num_parameters

    cdef void setSmallValue(self, float small_value):
        self.small_value = small_value

    cdef float getSmallValue(self):
        return self.small_value

    cdef void evolve(self, float timestep, float[:] pops_and_params):
        return

    cdef array.array calcQm(self, float[:] eqm_pops_and_params):
        """Not defined for this class"""
        return

    cdef unsigned getNumSpecies(self):
        return 0
    

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
        self.advection_class = array.array('I', [0, 0, 0])
        self.num_parameters = 2

    cdef void evolve(self, float timestep, float[:] pops_and_params):
        """No dynamics here"""
        
        return
    
    cdef array.array calcQm(self, float[:] eqm_pops_and_params):
        """Not defined for this class"""
        return

    cdef unsigned getNumSpecies(self):
        return 1
    
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
        self.advection_class = array.array('I', [0])
        self.num_parameters = 1
        self.growth_rate = 0.01 # fixed for this example

    cdef void evolve(self, float timestep, float[:] pops_and_params):
        # explicit time-stepping of logistic growth
        # as defined in the doco for "evolve",
        # pops_and_params[0] is the current population
        # pops_and_params[1] is the carrying capacity
        pops_and_params[0] = pops_and_params[0] + timestep * self.growth_rate * pops_and_params[0] * (1.0 - pops_and_params[0] / pops_and_params[1])

    cdef array.array calcQm(self, float[:] eqm_pops_and_params):
        """Not defined for this class"""
        return

    cdef unsigned getNumSpecies(self):
        return 1
    

cdef class CellDynamicsBeeton2_2(CellDynamicsBase):
    def __init__(self):
        super().__init__()
        self.num_populations = 2
        self.num_diffusing = 2 # both subspecies diffuse
        self.diffusing_indices = array.array('I', [0, 1])
        self.num_advecting = 2 # both subspecies advect (and there is just 1 advection class)
        self.advecting_indices = array.array('I', [0, 1])
        self.advection_class = array.array('I', [0, 0])
        self.num_parameters = 2 # carrying capacity Kx, and carrying capacity Ky
        # Default parameters in Equations (4.1) and (4.2) of Beeton et al, are as set in Fig5c
        self.mux = 0.7
        self.muy = 0.8
        self.gax = 1.0
        self.gay = 1.0
        self.axy = 0.4
        self.ayx = 0.4
        self.w = 0.05
        self.small = 0.1

    cdef void setMuX(self, float mux):
        self.mux = mux

    cdef void setMuY(self, float muy):
        self.muy = muy

    cdef void setGaX(self, float gax):
        self.gax = gax

    cdef void setGaY(self, float gay):
        self.gay = gay

    cdef void setAxy(self, float axy):
        self.axy = axy

    cdef void setAyx(self, float ayx):
        self.ayx = ayx

    cdef void setW(self, float w):
        self.w = w

    cdef void setSmall(self, float small):
        self.small = small

    cdef float getMuX(self):
        return self.mux

    cdef float getMuY(self):
        return self.muy

    cdef float getGaX(self):
        return self.gax

    cdef float getGaY(self):
        return self.gay

    cdef float getAxy(self):
        return self.axy

    cdef float getAyx(self):
        return self.ayx

    cdef float getW(self):
        return self.w

    cdef float getSmall(self):
        return self.small


    cdef void evolve(self, float timestep, float[:] pops_and_params):
        # explicit time-stepping of Beeton et al.'s equations
        # as defined in the doco for "evolve",
        # pops_and_params[0] is the current x
        # pops_and_params[1] is the current y
        # pops_and_params[2] is the carrying capacity Kx
        # pops_and_params[3] is the carrying capacity Ky
        cdef float x = pops_and_params[0]
        cdef float y = pops_and_params[1]
        cdef float kx = pops_and_params[2]
        cdef float ky = pops_and_params[3]
        cdef float f = - self.mux + (1.0 - (x + self.axy * y) / (kx + self.small)) * (x / (x + self.w * y + self.small)) * self.gax
        cdef float g = - self.muy + (1.0 - (self.ayx * x + y) / (ky + self.small)) * (self.gay + self.w * x / (x + self.w * y + self.small) * self.gax)
        # max means (x,y) doesn't stray outside physical quadrant, even if timestep is too large
        pops_and_params[0] = max(0.0, x + timestep * f * x)
        pops_and_params[1] = max(0.0, y + timestep * g * y)

    cdef array.array calcQm(self, float[:] eqm_pops_and_params):
        """Not defined for this class"""
        return

    cdef unsigned getNumSpecies(self):
        return 2
    
cdef class CellDynamicsMosquito23(CellDynamicsBase):
    def __init__(self):
        super().__init__()

        # num_sexes is always fixed to 2: male, female
        # this needs to be set early in constructor, so arrays (in particular genotypeRapidAccess) are sized correctly
        self.num_sexes = 2

        # Following parameters may be set by user.
        self.mu_larvae = 0.05
        self.mu_adult = 0.125
        self.fecundity = 9
        self.aging_rate = 0.1
        self.num_ages = 2
        self.num_species = 1
        self.accuracy = 0.95 # important to set this early in constructor so fecundity_proportion can be used to populate genotypeRapidAccess array
        self.num_species2 = 1

        # number of genotypes defaults to 3, which are: ww, Gw, GG
        self.setNumGenotypes(self.num_sexes, 3)

        self.num_parameters = 1 # the only spatially-varying quantity is the carrying capacity, K

        # Set default carrying capacities
        self.one_over_kk = array.clone(array.array('f', []), self.num_parameters, zero = True)
        for i in range(self.num_parameters):
            self.one_over_kk[i] = 1.0
            
        ## size inheritance correctly
        #self.inheritance_cube = array.clone(array.array('f', []), self.num_genotypes2 * self.num_genotypes, zero = False)
        #self.setInheritance()

        # allocate alpha array correctly, and set to the identity
        self.alpha = array.clone(array.array('f', []), 1, zero = True)
        self.setAlphaComponent(0, 0, 1.0)

        # allocate hyb array correctly, and set it to the identity
        self.hyb = array.clone(array.array('f', []), 1, zero = True)
        self.setHybridisationRate(0, 0, 0, 1.0)

        # allocate the mating array correctly, and set it to the identity
        # allocate hyb array correctly, and set it to the identity
        self.mating = array.clone(array.array('f', []), 1, zero = True)
        self.setMatingComponent(0, 0, 1.0)

        # allocate comp correctly
        self.comp = array.clone(array.array('f', []), self.num_species, zero = False)

        # allocate denom correctly
        self.denom = array.clone(array.array('f', []), self.num_species, zero = False)

        # allocate mat correctly
        self.mat = array.clone(array.array('f', []), self.num_populations * self.num_populations, zero = False)

        # allocate Xarray correctly
        self.Xarray = array.clone(array.array('f', []), self.num_populations, zero = False)

        # allocate the Runge-Kutta arrays correctly
        self.rk1 = array.clone(array.array('f', []), self.num_populations, zero = False)
        self.rk2 = array.clone(array.array('f', []), self.num_populations, zero = False)
        self.rk3 = array.clone(array.array('f', []), self.num_populations, zero = False)
        self.rk4 = array.clone(array.array('f', []), self.num_populations, zero = False)
        self.rky = array.clone(array.array('f', []), self.num_populations + self.num_parameters, zero = False)
        self.crky = self.rky

        # allocate rhs correctly
        self.rhs = array.clone(array.array('f', []), self.num_populations, zero = False)

        # allocate the change arrays correctly
        self.change = array.clone(array.array('f', []), self.num_populations, zero = False)
        self.cchange = self.change

        # provide defaults for some dummy floats, just to make the constructor fully construct everything
        self.species_stuff = 0.0
        self.genotype_stuff = 0.0
        self.tmp_float = 0.0
        self.xcol = 0.0

        # allocate the boolean arrays correctly
        self.species_present = array.clone(array.array('B', []), self.num_species, zero = False) 
        self.genotype_present = array.clone(array.array('B', []), self.num_genotypes, zero = False)

        # default to explicit_euler
        self.time_integration_method = 0

        # default to 1E-12 being the minimum timestep allowed
        self.min_dt = 1E-12

        # default to doing adaptive timestepping
        self.adaptive = 1
        
        # default to zero_cutoff = 1E-6, which means that if population numbers are less than 1E-6 at the end of a timestep, they will be set to zero
        self.zero_cutoff = 1E-6

        # default to min_cc = 1E-6, which means that if the carrying capacity is less than 1E-6 then no new larvae are produced
        self.setMinCarryingCapacity(1.0E-6)

        self.setInternalParameters(self.num_ages, self.num_species, self.accuracy)

    cdef void setInheritance(self):
        inheritance_list = [[[1., 0., 0.],# ww x ww
                             [ .5, .5, 0. ], # ww x Gw
                             [ 0., 1., 0. ]], # ww x GG
                            [[ .5, .5, 0. ], # Gw x ww
                             [ .25, .5, .25 ], # Gw x Gw
                             [ 0., .5, .5 ]], # Gw x GG
                            [[ 0., 1., 0. ], # GG x ww
                             [ 0., .5, .5 ], # GG x Gw
                             [ 0., 0., 1. ]]] # GG x GG
        cdef unsigned gt_father, gt_mother, gt_offspring
        for gt_father in range(self.num_genotypes):
            for gt_mother in range(self.num_genotypes):
                for gt_offspring in range(self.num_genotypes):
                    self.inheritance_cube.data.as_floats[gt_father + gt_mother * self.num_genotypes + gt_offspring * self.num_genotypes2] = inheritance_list[gt_father][gt_mother][gt_offspring]

    cdef void setInternalParameters(self, unsigned num_ages, unsigned num_species, float accuracy):
        self.num_ages = num_ages # age categories are: larvae0, larvae1, larvae2, ..., larvaeN, adults
        self.num_species = num_species
        self.num_species2 = num_species * num_species
        self.accuracy = accuracy # this must be done before the call to setGenotypeRapidAccess, below
        
        self.num_populations = self.num_ages * self.num_sexes * self.num_genotypes * self.num_species
        self.num_parameters = self.num_species

        self.one_over_kk = array.clone(array.array('f', []), self.num_parameters, zero = True)
        for i in range(self.num_parameters):
            self.one_over_kk[i] = 1.0
            
        # allocate the mat array correctly
        self.mat = array.clone(array.array('f', []), self.num_populations * self.num_populations, zero = False)
        # allocate Xarray correctly
        self.Xarray = array.clone(array.array('f', []), self.num_populations, zero = False)
        # allocate the Runge-Kutta arrays correctly
        self.rk1 = array.clone(array.array('f', []), self.num_populations, zero = False)
        self.rk2 = array.clone(array.array('f', []), self.num_populations, zero = False)
        self.rk3 = array.clone(array.array('f', []), self.num_populations, zero = False)
        self.rk4 = array.clone(array.array('f', []), self.num_populations, zero = False)
        self.rky = array.clone(array.array('f', []), self.num_populations + self.num_parameters, zero = False)
        self.crky = self.rky
        # allocate the rhs array correctly
        self.rhs = array.clone(array.array('f', []), self.num_populations, zero = False)
        # allocate the change arrays correctly
        self.change = array.clone(array.array('f', []), self.num_populations, zero = False)
        self.cchange = self.change

        # allocate the boolean arrays correctly
        self.species_present = array.clone(array.array('B', []), self.num_species, zero = False) 

        # set diffusing and advecting information
        self.num_diffusing = self.num_sexes * self.num_genotypes * self.num_species # only adults diffuse
        self.num_advecting = self.num_sexes * self.num_genotypes * self.num_species # only adults advect
        self.diffusing_indices = array.clone(array.array('I', []), self.num_diffusing, zero = False)
        self.advecting_indices = array.clone(array.array('I', []), self.num_advecting, zero = False)
        self.advection_class = array.array('I', [0] * self.num_advecting)
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

        # because genotypeRapidAccess depends on accuracy (through fecundity_proportion), we need to set it:
        self.setGenotypeRapidAccess()

    cdef setMinimumDt(self, float value):
        self.min_dt = value

    cdef float getMinimumDt(self):
        return self.min_dt

    cdef setAdaptive(self, unsigned value):
        self.adaptive = value

    cdef unsigned getAdaptive(self):
        return self.adaptive

    cdef setZeroCutoff(self, float value):
        self.zero_cutoff = value

    cdef float getZeroCutoff(self):
        return self.zero_cutoff

    cdef setMinCarryingCapacity(self, float value):
        self.min_cc = value
        self.one_over_min_cc = 1.0 / self.min_cc

    cdef float getMinCarryingCapacity(self):
        return self.min_cc

    cdef setAlphaComponent(self, unsigned sp0, unsigned sp1, float value):
        if sp0 >= self.num_species or sp1 >= self.num_species:
            raise ValueError("sp0 " + str(sp0) + " and sp1 " + str(sp1) + " must be less than the number of species, " + str(self.num_species))
        self.alpha.data.as_floats[sp0 + sp1 * self.num_species] = value

    cdef setHybridisationRate(self, unsigned species_father, unsigned species_mother, unsigned species_offspring, float value):
        if species_father >= self.num_species or species_mother >= self.num_species or species_offspring >= self.num_species:
            raise ValueError("All species numbers, " + str(species_father) + ", " + str(species_mother) + ", " + str(species_offspring) + " must be less than the number of species, " + str(self.num_species))
        self.hyb.data.as_floats[species_father + species_mother * self.num_species + species_offspring * self.num_species2] = value

    def getHybridisationRateFromPython(self, unsigned species_father, unsigned species_mother, unsigned species_offspring):
        """Returns the hybridisation rate for given father, mother and offspring.  This is a slow interface: use getAlphaComponent from all cython code"""
        
        if species_father >= self.num_species or species_mother >= self.num_species or species_offspring >= self.num_species:
            raise ValueError("All species numbers, " + str(species_father) + ", " + str(species_mother) + ", " + str(species_offspring) + " must be less than the number of species, " + str(self.num_species))
        return self.getHybridisationRate(species_father, species_mother, species_offspring)

    cdef setMatingComponent(self, unsigned species_father, unsigned species_mother, float value):
        if species_father >= self.num_species or species_mother >= self.num_species:
            raise ValueError("species_father " + str(species_father) + " and species_mother " + str(species_mother) + " must be less than the number of species, " + str(self.num_species))
        self.mating.data.as_floats[species_father + species_mother * self.num_species] = value

    cdef setFitnessComponent(self, unsigned genotype, float value):
        if genotype >= self.num_genotypes:
            raise ValueError("genotype " + str(genotype) + " must be less than the number of genotypes " + str(self.num_genotypes))
        self.fitness.data.as_floats[genotype] = value
        self.setGenotypeRapidAccess()

    cdef void setGenotypeRapidAccess(self):
        if len(self.genotypeRapidAccess) != self.num_sexes * self.num_genotypes2:
            raise ValueError("self.genotypeRapidAccess.size() != self.num_sexes * self.num_genotypes2 in setFitnessComponent.  Perhaps the code does not set num_sexes correctly")
        cdef unsigned ind1, sex, gtf, gtm
        ind1 = -1
        for sex in range(self.num_sexes): # row in M
            for gtf in range(self.num_genotypes): # female genotype
                for gtm in range(self.num_genotypes): # male genotype
                    ind1 += 1 # = sex * self.num_genotypes2 + gtf * self.num_genotypes + gtm
                    self.genotypeRapidAccess.data.as_floats[ind1] = self.fecundity_proportion(sex, gtf, gtm) * self.getFitnessComponent(gtm)


    def getFitnessComponentFromPython(self, unsigned genotype):
        """Python interface for getting a component of the Fitness vector.  This is a slow interface: use getFitnessComponent from all cython code"""
        
        if genotype >= self.num_genotypes:
            raise ValueError("Genotype " + str(genotype) + " must be less than the number of genotypes, " + str(self.num_genotypes))
        return self.getFitnessComponent(genotype)

    cdef void setMuLarvae(self, float mu_larvae):
        self.mu_larvae = mu_larvae

    cdef float getMuLarvae(self):
        return self.mu_larvae
    
    cdef void setMuAdult(self, float mu_adult):
        self.mu_adult = mu_adult

    cdef float getMuAdult(self):
        return self.mu_adult
    
    cdef void setFecundity(self, float fecundity):
        self.fecundity = fecundity

    cdef float getFecundity(self):
        return self.fecundity
    
    cdef void setAgingRate(self, float aging_rate):
        self.aging_rate = aging_rate

    cdef float getAgingRate(self):
        return self.aging_rate
    
    cdef void setNumAges(self, unsigned num_ages):
        self.setInternalParameters(num_ages, self.num_species, self.accuracy)

    cdef unsigned getNumAges(self):
        return self.num_ages
    
    cdef void setNumSpecies(self, unsigned num_species):
        self.setInternalParameters(self.num_ages, num_species, self.accuracy)
        self.alpha = array.clone(array.array('f', []), num_species * num_species, zero = True)
        cdef unsigned species
        for species in range(num_species):
            self.setAlphaComponent(species, species, 1.0)
        self.hyb = array.clone(array.array('f', []), num_species * num_species * num_species, zero = True)
        for species in range(num_species):
            self.setHybridisationRate(species, species, species, 1.0)
        self.mating = array.clone(array.array('f', []), num_species * num_species, zero = True)
        for species in range(num_species):
            self.setMatingComponent(species, species, 1.0)
        self.comp = array.clone(array.array('f', []), num_species, zero = False)
        self.denom = array.clone(array.array('f', []), num_species, zero = False)

    cdef void setNumGenotypes(self, unsigned num_sexes, unsigned num_genotypes):
        if num_sexes != self.num_sexes:
            raise ValueError("setNumGenotypes: num_sexes != self.num_sexes.  " + str(num_sexes) + "!=" + str(self.num_sexes) + ".  Probably there is a bug in the code")

        self.num_genotypes = num_genotypes
        self.num_genotypes2 = self.num_genotypes * self.num_genotypes

        # correctly size arrays that depend on num_genotypes
        self.inheritance_cube = array.clone(array.array('f', []), self.num_genotypes2 * self.num_genotypes, zero = False)
        self.genotypeRapidAccess = array.clone(array.array('f', []), num_sexes * self.num_genotypes2, zero = False)
        self.genotype_present = array.clone(array.array('B', []), num_genotypes, zero = False)
        self.fitness = array.clone(array.array('f', []), num_genotypes, zero = True)

        # calculate things that depend on num_genotypes
        self.setInheritance()
        for genotype in range(num_genotypes):
            self.setFitnessComponent(genotype, 1.0)

    def getAlphaComponentFromPython(self, unsigned sp0, unsigned sp1):
        """Python interface for getting a component of the alpha matrix (inter-specific competition).  This is a slow interface: use getAlphaComponent from all cython code"""
        
        return self.getAlphaComponent(sp0, sp1)

    def getMatingComponentFromPython(self, unsigned species_father, unsigned species_mother):
        """Python interface for getting a component of the mating matrix (relative probability of male mating with female).  This is a slow interface: use getMatingComponent from all cython code"""
        
        return self.getMatingComponent(species_father, species_mother)

    def getInheritanceFromPython(self, unsigned gt_father, unsigned gt_mother, unsigned gt_offspring):
        """Python interface for getting a component of the inheritance cube.  This is a slow interface: use getInheritance from all cython code"""
        
        if gt_father >= self.num_genotypes or gt_mother >= self.num_genotypes or gt_offspring >= self.num_genotypes:
            raise ValueError("All genotypes, " + str(gt_father) + ", " + str(gt_mother) + ", " + str(gt_offspring) + " must be less than the number of genotypes, " + str(self.num_genotypes))
        return self.getInheritance(gt_father, gt_mother, gt_offspring)

    cdef unsigned getNumSpecies(self):
        return self.num_species
    
    cdef void setAccuracy(self, float accuracy):
        self.setInternalParameters(self.num_ages, self.num_species, accuracy)

    cdef float getAccuracy(self):
        return self.accuracy

    cdef void computeRHS(self, float[:] x):
        """Evaluates d(populations)/dt"""

        # these indices are dummy indices (that are looped over) hence the subscript _d
        cdef unsigned age_d, sex_d, gt_d, sp_d, ind_d, tmp_d
        # these indices are summed over in the age=0 contributions
        cdef unsigned gtf, spf, gtm, spm
        # these indices are typically left-hand-side indices
        cdef unsigned age, sex, gt, sp
        # utility indices
        cdef unsigned ind, ind_c, ind0, ind1, ind2, cidx
        # indices into matrix M
        cdef unsigned col, row
        # index into mat
        cdef unsigned ind_mat
        # end of for-loops in newborn larvae competition code
        cdef unsigned end_index_for_competition = self.num_ages - 1 if self.num_ages > 1 else 1
	
        #array.zero(self.mat) # no longer using matrix multiplication, just putting results straight into RHS
        array.zero(self.rhs)
        
        # Calculate indicators whether a given species or genotype is present
        array.zero(self.species_present)
        ind_d = 0
        for sp in range(self.num_species):
            for gt_d in range(self.num_genotypes):
                for sex_d in range(self.num_sexes):
                    for age_d in range(self.num_ages):
                        ind_d = self.getIndex(sp, gt_d, sex_d, age_d)
                        if x[ind_d] > 0:
                            self.species_present.data.as_uchars[sp] = 1
                            break
                    if x[ind_d] > 0:
                        break
                if x[ind_d] > 0:
                    break
                
        array.zero(self.genotype_present)
        for gt in range(self.num_genotypes):
            for sp_d in range(self.num_species):
                for sex_d in range(self.num_sexes):
                    for age_d in range(self.num_ages):
                        ind_d = self.getIndex(sp_d, gt, sex_d, age_d)
                        if x[ind_d] > 0:
                            self.genotype_present.data.as_uchars[gt] = 1
                            break
                    if x[ind_d] > 0:
                        break
                if x[ind_d] > 0:
                    break

        # newborn larvae
        cdef unsigned newborn_calcs_needed = 0
        for p in range(self.num_parameters):
            if self.one_over_kk[p] < self.one_over_min_cc:
                newborn_calcs_needed = 1
                break
        if newborn_calcs_needed == 1:
            array.zero(self.comp)
            array.zero(self.denom)

            # first calculate all comp.data and denom.data where needed
            for sp in range(self.num_species):
                if self.species_present.data.as_uchars[sp] == 1: # can never generate a species if it does not exist
                    if self.num_parameters == 1: # only using one CC
                        cidx = 0
                    else:
                        cidx = sp # one CC for each species
                    if self.one_over_kk[cidx] < self.one_over_min_cc: # if there are any newborns of species sp (or any species) born at all
                        # define competition
                        for age_d in range(end_index_for_competition):
                            for sex_d in range(self.num_sexes):
                                for gt_d in range(self.num_genotypes):
                                    if self.genotype_present.data.as_uchars[gt_d] == 1: # if not present then x[ind] = 0
                                        for sp_d in range(self.num_species):
                                            if self.species_present.data.as_uchars[sp_d] == 1:  # if not present then x[ind] = 0
                                                ind_d = self.getIndex(sp_d, gt_d, sex_d, age_d)
                                                self.comp.data.as_floats[sp] = self.comp.data.as_floats[sp] + self.getAlphaComponent(sp, sp_d) * x[ind_d]
                        self.comp.data.as_floats[sp] = max(0.0, 1.0 - self.comp.data.as_floats[sp] * self.one_over_kk[cidx])

                        # define the denominator term
                        age_d = self.num_ages - 1 # adult
                        sex_d = 0 # male
                        for gt_d in range(self.num_genotypes):
                            if self.genotype_present.data.as_uchars[gt_d] == 1: # if not present then x[ind] = 0
                                for sp_d in range(self.num_species):
                                    if self.species_present.data.as_uchars[sp_d] == 1:  # if not present then x[ind] = 0
                                        ind_d = self.getIndex(sp_d, gt_d, sex_d, age_d)
                                        # note: here sp = female species
                                        self.denom.data.as_floats[sp] = self.denom.data.as_floats[sp] + self.getMatingComponent(sp_d, sp) * self.getFitnessComponent(gt_d) * x[ind_d]
                        # form the reciprocal of denom.  This is so that C doesn't have to check for division-by-zero in the big loops below
                        if self.denom.data.as_floats[sp] <= self.zero_cutoff:
                            # there must be zero adult males.  I presume this means there will be zero eggs layed, so setting denom=0 achieves this
                            self.denom.data.as_floats[sp] = 0.0
                        else:
                            self.denom.data.as_floats[sp] = 1.0 / self.denom.data.as_floats[sp]

            # now define competition

            for sp in range(self.num_species):
                if self.species_present.data.as_uchars[sp] == 1: # self.comp[sp] will be zero
                    if self.num_parameters == 1: # only using one CC
                        cidx = 0
                    else:
                        cidx = sp # one CC for each species
                    if self.one_over_kk[cidx] < self.one_over_min_cc: # if there are any newborns of species sp (or any species) born at all
                        # define competition
                        # now work out the contributions to the newborn ODEs
                        # In the following, mat is a 1D array (for efficiency)
                        # If visualised as a matrix, M, where the ODE is dot{X} = MX (X is a column vector) then:
                        #  - given an age, species, genotype and species, the index in the X vector is given by the function getIndex (inlined for efficiency)
                        #  - given a component, M_ij, the index in mat is j + i * self.num_populations
                        age = 0 # only newborn row in M
                        # Note: the species loop above (over sp) is a row in M
                        
                        for sex in range(self.num_sexes): # row in M
                            for gt in range(self.num_genotypes): # row in M
                                row = self.getIndex(sp, gt, sex, age)
                                # now want to set the column in M corresonding to female adults of genotype gtf and species spf
                                for gtf in range(self.num_genotypes): # female genotype
                                    if self.genotype_present.data.as_uchars[gtf] == 1: # if not present then x[col] will be zero
                                        for spf in range(self.num_species): # female species
                                            if self.species_present.data.as_uchars[spf] == 1: # if not present then x[col] will be zero
                                                col = self.getIndex(spf, gtf, 1, self.num_ages - 1) # species=spf, genotype=gtf, sex=female, age=adult
                                                self.xcol = x[col]
                                                self.tmp_float = 0.0
                                                ind_mat = col + row * self.num_populations  # index into mat corresponding to the row, and the aforementioned adult female
                                                for gtm in range(self.num_genotypes): # male genotype
                                                    if self.genotype_present.data.as_uchars[gtm] == 1: # if not present then x[ind] will be zero
                                                        ind1 = sex * self.num_genotypes2 + gtf * self.num_genotypes + gtm
                                                        self.genotype_stuff = self.getInheritance(gtm, gtf, gt)
                                                        for spm in range(self.num_species): # male species
                                                            if self.species_present.data.as_uchars[spm] == 1:  # if not present then x[ind] will be zero
                                                                ind = self.getIndex(spm, gtm, 0, self.num_ages - 1) # species=spm, genotype=gtm, sex=male, age=adult
                                                                self.species_stuff = self.getHybridisationRate(spm, spf, sp) * self.getMatingComponent(spm, spf)
                                                                self.tmp_float = self.tmp_float + self.species_stuff * self.genotypeRapidAccess.data.as_floats[ind1] * self.genotype_stuff * x[ind] * self.xcol
                                                # multiply rhs by things that don't depend on gtm or spm
                                                self.tmp_float *= self.comp.data.as_floats[sp] * self.fecundity * self.denom.data.as_floats[spf]
                                                self.rhs.data.as_floats[row] = self.rhs.data.as_floats[row] + self.tmp_float
            
        # mortality, and aging into/from neighbouring age brackets
        for sex in range(self.num_sexes):
            for gt in range(self.num_genotypes):
                if self.genotype_present.data.as_uchars[gt] == 1: # if not present then cannot create the genotype from current x (next timestep there may be juveniles of this gt that can age, but not this timestep)
                    for sp in range(self.num_species):
                        if self.species_present.data.as_uchars[sp] == 1: # if not present then cannot create the species out of nothing
                            age = 0
                            row = self.getIndex(sp, gt, sex, age)
                            ind = row + self.num_populations * row # diagonal entry
                            if self.num_ages > 1:
                                self.rhs.data.as_floats[row] = self.rhs.data.as_floats[row] - (self.mu_larvae + self.aging_rate)*x[row] # mortality + aging to older bracket
                                
                            for age in range(1, self.num_ages - 1):
                                row = self.getIndex(sp, gt, sex, age)
                                ind = row + self.num_populations * row # diagonal entry
                                self.rhs.data.as_floats[row] = self.rhs.data.as_floats[row] - (self.mu_larvae + self.aging_rate)*x[row] # mortality + aging to older bracket
                                col = self.getIndex(sp, gt, sex, age - 1)
                                ind = col + self.num_populations * row # below diagonal
                                self.rhs.data.as_floats[row] = self.rhs.data.as_floats[row] + self.aging_rate*x[col] # mortality + aging to older bracket
                            age = self.num_ages - 1
                            row = self.getIndex(sp, gt, sex, age)
                            ind = row + self.num_populations * row # diagonal entry
                            self.rhs.data.as_floats[row] = self.rhs.data.as_floats[row] - self.mu_adult * x[row]# mortality
                            if self.num_ages > 1:
                                col = self.getIndex(sp, gt, sex, age - 1)
                                ind = col + self.num_populations * row # below diagonal
                                self.rhs.data.as_floats[row] = self.rhs.data.as_floats[row] + self.aging_rate*x[col]



    cdef void evolve(self, float timestep, float[:] pops_and_params):
        cdef unsigned ind
        cdef unsigned sp

        for sp in range(self.num_parameters): # for each defined carrying capacity
            if pops_and_params[self.num_populations + sp] <= self.min_cc:
                self.one_over_kk[sp] = 2 * self.one_over_min_cc # signal to computeRHS that the carrying capacity has fallen below min_cc, so no newborns will be produced
            else:
                self.one_over_kk[sp] = 1.0 / pops_and_params[self.num_populations + sp]        

        cdef float time_done = 0.0
        cdef float dt = timestep
        cdef float new_dt = timestep
        while time_done < timestep:
            # try timestep dt
            self.popChange(dt, pops_and_params, self.cchange)
            new_dt = timestep # biggest value that ever needs to be used
            for ind in range(self.num_populations):
                if pops_and_params[ind] + self.cchange[ind] < -self.zero_cutoff:
                    # the change will send the population negative, which is deemed to be erroneous
                    # and due to dt being too big.  Using explicit-Euler timestepping, if we choose
                    # new_dt = dt * pops_and_params[ind] / (-cchange[ind]), then running popChange again
                    # will give exactly cchange[ind] = -pops_and_params[ind].  So multiply this
                    # new_dt by 0.9 to give a factor of safety
                    new_dt = min(new_dt, - 0.9 * dt * pops_and_params[ind] / self.cchange[ind])
                elif pops_and_params[ind] + self.cchange[ind] < 0.0:
                    # the change will send the population negative, but only by a tiny amount
                    # since it is greater than -self.zero_cutoff.
                    # This is probably due to a precision loss problem in Runge-Kutta.
                    # So set self.cchange so that pops_and_params[ind] + self.cchange[ind] = 0
                    self.cchange[ind] = -pops_and_params[ind]
            if self.adaptive == 0 or new_dt == timestep:
                # not doing adaptivity, or none of the populations went negative (less than -self.zero_cutoff), so success:
                for ind in range(self.num_populations):
                    pops_and_params[ind] = pops_and_params[ind] + self.cchange[ind]
                time_done = time_done + dt
                # can also increase dt for next time around (if there is a next time around)
                dt = min(1.1 * dt, timestep - time_done)
            else:
                # at least one of the populations went negative, so must re-do
                dt = new_dt
                if dt < self.min_dt:
                    sys.stderr.write("Minimum dt reached.  Exiting\n")
                    sys.exit(1)
        
        for ind in range(self.num_populations):
            if pops_and_params[ind] < self.zero_cutoff:
                pops_and_params[ind] = 0.0

    cdef array.array calcQm(self, float[:] eqm_pops_and_params):
        """Not defined for this class"""
        return

    cdef void popChange(self, float timestep, float[:] current_pops_and_params, float[:] cchange):
        cdef unsigned ind

        if self.time_integration_method == 0:
            self.computeRHS(current_pops_and_params)
            for ind in range(self.num_populations):
                cchange[ind] = timestep * self.rhs.data.as_floats[ind]

        elif self.time_integration_method == 1:
            # copy into numpy array xx for use in self.fun
            xx = np.ones(self.num_populations)
            for ind in range(self.num_populations):
                xx[ind] = current_pops_and_params[ind]
            # solve
            sol = solve_ivp(self.fun_for_scipy, [0.0, timestep], xx)
            # copy back
            for ind in range(self.num_populations):
                cchange[ind] = sol.y[ind, -1] - current_pops_and_params[ind]

        elif self.time_integration_method == 2:
            # step 1
            self.computeRHS(current_pops_and_params)
            for ind in range(self.num_populations):
                self.rk1.data.as_floats[ind] = timestep * self.rhs.data.as_floats[ind]
            # step 2
            for ind in range(self.num_populations):
                self.crky[ind] = current_pops_and_params[ind] + 0.5 * self.rk1.data.as_floats[ind]
            self.crky[self.num_populations:] = current_pops_and_params[self.num_populations:] # the carrying capacities
            self.computeRHS(self.crky)
            for ind in range(self.num_populations):
                self.rk2.data.as_floats[ind] = timestep * self.rhs.data.as_floats[ind]
            # step 3
            for ind in range(self.num_populations):
                self.crky[ind] = current_pops_and_params[ind] + 0.5 * self.rk2.data.as_floats[ind]
            self.crky[self.num_populations:] = current_pops_and_params[self.num_populations:] # the carrying capacities
            self.computeRHS(self.crky)
            for ind in range(self.num_populations):
                self.rk3.data.as_floats[ind] = timestep * self.rhs.data.as_floats[ind]
            # step 4
            for ind in range(self.num_populations):
                self.crky[ind] = current_pops_and_params[ind] + self.rk3.data.as_floats[ind]
            self.crky[self.num_populations:] = current_pops_and_params[self.num_populations:] # the carrying capacities
            self.computeRHS(self.crky)
            for ind in range(self.num_populations):
                self.rk4.data.as_floats[ind] = timestep * self.rhs.data.as_floats[ind]
            # put it all together
            for ind in range(self.num_populations):
                cchange[ind] = (1.0 / 6.0) * (self.rk1.data.as_floats[ind] + 2 * self.rk2.data.as_floats[ind] + 2 * self.rk3.data.as_floats[ind] + self.rk4.data.as_floats[ind])


    cdef setTimeIntegrationMethod(self, str method):
        if method == "explicit_euler":
            self.time_integration_method = 0
        elif method == "solve_ivp":
            self.time_integration_method = 1
        elif method == "runge_kutta4":
            self.time_integration_method = 2
        elif method == "stochastic":
            if self.__class__.__name__ != "CellDynamicsMosquito23G":
                raise ValueError("Stochastic time integration can only be used in CellDynamicsMosquito23G")
            self.time_integration_method = 3
        else:
            raise ValueError("Time integration method " + method + " not supported")
    
    def fun_for_scipy(self, t, y):
        """Evaluates d(populations)/dt"""
        
        # size cXarray correctly
        self.cXarray = self.Xarray
        # copy the "y" data into cXarray and then call computeRHS
        cdef unsigned ind
        for ind in range(self.num_populations):
            self.cXarray[ind] = y[ind]
        self.computeRHS(self.cXarray)
        dXdt = np.zeros((self.num_populations))
        for ind in range(self.num_populations):
            dXdt[ind] = self.rhs.data.as_floats[ind]
        return dXdt

    cdef float fecundity_proportion(self, unsigned offspring_sex, unsigned mother_gt, unsigned father_gt):
        return 0.5 if (offspring_sex == 0 or father_gt == 0) else 0.5 * (1.0 / self.accuracy - 1.0) # reduced fecundity version


cdef class CellDynamicsMosquito23F(CellDynamicsMosquito23):
    def __init__(self):
        super().__init__()

    cdef float fecundity_proportion(self, unsigned offspring_sex, unsigned mother_gt, unsigned father_gt):
        if (father_gt == 0):
            if (mother_gt == 0):
                return 0.5
            else:
                if (offspring_sex == 0): # female sex bias if modified mother, wildtype father
                    return 0.45
                else:
                    return 0.55
        else:
            if (offspring_sex == 0): # male sex bias if modified father (any mother)
                return 0.95
            else:
                return 0.05

cdef class CellDynamicsMosquito23G(CellDynamicsMosquito23F):
    def __init__(self):
        super().__init__()
        # default to NOT doing adaptive timestepping
        self.adaptive = 0

    cdef void popChange(self, float timestep, float[:] current_pops_and_params, float[:] cchange):
        if self.time_integration_method == 0 or self.time_integration_method == 1 or self.time_integration_method == 2:
            super().popChange(timestep, current_pops_and_params, cchange)
                
        elif self.time_integration_method == 3:
            self.computeRHS_stoc(current_pops_and_params)
            for ind in range(self.num_populations):
                cchange[ind] = self.rhs.data.as_floats[ind]

    cdef void computeRHS_stoc(self, float[:] x):
        """Evaluates change in populations for a timestep (assuming dt = 1 day)"""

        # these indices are dummy indices (that are looped over) hence the subscript _d
        cdef unsigned age_d, sex_d, gt_d, sp_d
        # these indices are summed over in the age=0 contributions
        cdef unsigned gtf, spf, gtm, spm
        # these indices are typically left-hand-side indices
        cdef unsigned age, sex, gt, sp
        # utility indices
        cdef unsigned ind, ind_c, ind0, ind1
        # indices into matrix M
        cdef unsigned col, row
        # index into mat
        cdef unsigned ind_mat
        # end of for-loops in newborn larvae competition code
        cdef unsigned end_index_for_competition = self.num_ages - 1 if self.num_ages > 1 else 1

        array.zero(self.mat)
        array.zero(self.rhs)
        
        # newborn larvae (TODO: incorporate multiple species as in CellDynamicsMosquito.computeRHS)
        if self.one_over_kk < self.one_over_min_cc:

            # ANDY QUESTION: duplicated block?
            # define competition
            array.zero(self.comp)
            for age_d in range(end_index_for_competition):
                for sex_d in range(self.num_sexes):
                    for gt_d in range(self.num_genotypes):
                        for sp_d in range(self.num_species):
                            ind_d = self.getIndex(sp_d, gt_d, sex_d, age_d)
                            for sp in range(self.num_species):
                                self.comp.data.as_floats[sp] = self.comp.data.as_floats[sp] + self.getAlphaComponent(sp, sp_d) * x[ind_d]
            for sp in range(self.num_species):
                self.comp.data.as_floats[sp] = max(0.0, 1.0 - self.comp.data.as_floats[sp] * self.one_over_kk)

            # ANDY QUESTION: duplicated block?
            # define the denominator term
            array.zero(self.denom)
            age_d = self.num_ages - 1 # adult
            sex_d = 0 # male
            for gt_d in range(self.num_genotypes):
                for sp_d in range(self.num_species):
                    ind_d = self.getIndex(sp_d, gt_d, sex_d, age_d)
                    for sp in range(self.num_species): # sp = female species
                        self.denom.data.as_floats[sp] = self.denom.data.as_floats[sp] + self.getMatingComponent(sp_d, sp) * x[ind_d]
            # form the reciprocal of denom.  This is so that C doesn't have to check for division-by-zero in the big loops below
            for sp in range(self.num_species):
                if self.denom.data.as_floats[sp] <= self.zero_cutoff:
                    # there must be zero adult males.  I presume this means there will be zero eggs layed, so setting denom=0 achieves this
                    self.denom.data.as_floats[sp] = 0.0
                else:
                    self.denom.data.as_floats[sp] = 1.0 / self.denom.data.as_floats[sp]


            # now work out the contributions to the newborn ODEs
            # In the following, mat is a 1D array (for efficiency)
            # If visualised as a matrix, M, where the ODE is dot{X} = MX (X is a column vector) then:
            #  - given an age, species, genotype and species, the index in the X vector is given by the function getIndex (inlined for efficiency)
            #  - given a component, M_ij, the index in mat is j + i * self.num_populations
            age = 0 # only newborn row in M
            for sex in range(self.num_sexes): # row in M
                for gt in range(self.num_genotypes): # row in M
                    for sp in range(self.num_species): # row in M
                        row = self.getIndex(sp, gt, sex, age)
                        # now want to set the column in M corresonding to female adults of genotype gtf and species spf
                        for gtf in range(self.num_genotypes): # female genotype
                            for spf in range(self.num_species): # female species
                                col = self.getIndex(spf, gtf, 1, self.num_ages - 1) # species=spf, genotype=gtf, sex=female, age=adult
                                ind_mat = col + row * self.num_populations  # index into mat corresponding to the row, and the aforementioned adult female
                                for gtm in range(self.num_genotypes): # male genotype
                                    for spm in range(self.num_species): # male species
                                        ind = self.getIndex(spm, gtm, 0, self.num_ages - 1) # species=spm, genotype=gtm, sex=male, age=adult
                                        self.mat.data.as_floats[ind_mat] += self.getHybridisationRate(spm, spf, sp) * self.getInheritance(gtm, gtf, gt) * self.getMatingComponent(spm, spf) * x[ind] * self.fecundity_proportion(sex, gtf, gtm)
                                # multiply mat by things that don't depend on gtm or spm
                                self.mat.data.as_floats[ind_mat] *= self.comp.data.as_floats[sp] * self.fecundity * self.denom.data.as_floats[spf]
                                # ANDY QUESTION: Is this the only change from the mosquito23 version?
                                self.rhs.data.as_floats[row] += poisson(self.mat.data.as_floats[ind_mat] * x[col])

        # mortality, and aging into/from neighbouring age brackets
        for sex in range(self.num_sexes):
            for gt in range(self.num_genotypes):
                for sp in range(self.num_species):
                    age = 0
                    row = self.getIndex(sp, gt, sex, age)
                    ind = row + self.num_populations * row # diagonal entry
                    if self.num_ages > 1:
                        self.rhs.data.as_floats[row] -= binomial(int(x[row]), 1 - exp(-self.mu_larvae)) # mortality
                        
                        # ANDY QUESTION: correctly re-setting mat[0] here?
                        self.mat.data.as_floats[0] = binomial(int(x[row]), 1 - exp(-self.aging_rate))
                        self.rhs.data.as_floats[row] -= self.mat.data.as_floats[0] # aging to older bracket
                        
                    for age in range(1, self.num_ages - 1):
                        row = self.getIndex(sp, gt, sex, age)
                        self.rhs.data.as_floats[row] -= binomial(int(x[row]), 1 - exp(-self.mu_larvae))
                        # ANDY QUESTION: Should [0] be [age - 1] ?
                        self.rhs.data.as_floats[row] += self.mat.data.as_floats[0] # contribution from younger age bracket
                        # ANDY QUESTION: correctly re-setting mat[0] here?
                        self.mat.data.as_floats[0] = binomial(int(x[row]), 1 - exp(-self.aging_rate))
                        self.rhs.data.as_floats[row] -= self.mat.data.as_floats[0]
                        
                    age = self.num_ages - 1
                    row = self.getIndex(sp, gt, sex, age)
                    self.rhs.data.as_floats[row] -= binomial(int(x[row]), 1 - exp(-self.mu_adult)) # mortality
                    if self.num_ages > 1:
                        # ANDY QUESTION: Should [0] be [age - 1] ?
                        self.rhs.data.as_floats[row] += self.mat.data.as_floats[0] # contribution from younger age bracket

        #for row in range(self.num_populations):
            #for col in range(self.num_populations):
                #self.rhs.data.as_floats[row] += self.mat.data.as_floats[col + self.num_populations * row] * x[col]

cdef class CellDynamicsMosquito26(CellDynamicsMosquito23):
    def __init__(self):
        super().__init__()
        
        self.accuracy = 0.5 # no sex bias
        
        self.setNumSpecies(2)

	# setNumSpecies initialises alpha to the identity.  Modify it as follows:
        self.setAlphaComponent(0, 1, 0.4)
        self.setAlphaComponent(1, 0, 0.4)

	# setNumSpecies initialises HybridisationRate(i, i, i) = 1, and all other components zero.  Modify it:
        self.setHybridisationRate(0, 1, 0, 0.5) # cross matings produce 50/50 gambiae + coluzzii
        self.setHybridisationRate(0, 1, 1, 0.5) 
        self.setHybridisationRate(1, 0, 0, 0.5)
        self.setHybridisationRate(1, 0, 1, 0.5) 

	# setNumSpecies initialises MatingComponent(i, i) = 1, and all other components zero.  Modify it:
        self.setMatingComponent(0, 1, 0.01) # relative cross mating w = 0.01
        self.setMatingComponent(1, 0, 0.01)

        # Probabilities used in the Inheritance cube
        self.w_prob = 0.5 * (1 - 0.995)
        self.c_prob = 0.5 * (1 + 0.995 * (1 - 0.02) * (1 - 0.0001))
        self.r_prob = 1. - self.w_prob - self.c_prob
        
        self.setNumGenotypes(self.num_sexes, 6)  # ww, wc, wr, cc, cr, rr.  This also initialises the Inheritance cube with the above probabilities

	# setNumGenotypes initialises FitnessComponent to 1 for all genotypes.  Modify it:
        h_e = h_n = 0.5
        s_e = 0.1
        s_n = 0.05
        self.setFitnessComponents26(h_e, h_n, s_e, s_n)

	# Each species has its own carrying capacity:
        self.num_parameters = self.num_species

        self.setInternalParameters(self.num_ages, self.num_species, self.accuracy)
        
    cdef void setInheritance(self):
        inheritance_list = [[[1., 0., 0., 0., 0., 0.], # ww x ww
                             [self.w_prob, self.c_prob, self.r_prob, 0., 0., 0.], # ww x wc
                             [0.5, 0., 0.5, 0., 0., 0.], # ww x wr
                             [0., 1., 0., 0., 0., 0.], # ww x cc
                             [0., 0.5, 0.5, 0., 0., 0.], # ww x cr
                             [0., 0., 1., 0., 0., 0.]], # ww x rr
                            [[self.w_prob, self.c_prob, self.r_prob, 0., 0., 0.], # wc x ww
                             [self.w_prob*self.w_prob, 2.*self.w_prob*self.c_prob, 2.*self.w_prob*self.r_prob, self.c_prob*self.c_prob, 2.*self.c_prob*self.r_prob, self.r_prob*self.r_prob], # wc x wc
                             [0.5*self.w_prob, 0.5*self.c_prob, 0.5*(self.w_prob + self.r_prob), 0., 0.5*self.c_prob, 0.5*self.r_prob], # wc x wr
                             [0., self.w_prob, 0., self.c_prob, self.r_prob, 0.], # wc x cc
                             [0., 0.5*self.w_prob, 0.5*self.w_prob, 0.5*self.c_prob, 0.5*(self.c_prob + self.r_prob), 0.5*self.r_prob], # wc x cr
                             [0., 0., self.w_prob, 0., self.c_prob, self.r_prob]], # wc x rr
                            [[0.5, 0., 0.5, 0., 0., 0.], # wr x ww
                             [0.5*self.w_prob, 0.5*self.c_prob, 0.5*(self.w_prob + self.r_prob), 0., 0.5*self.c_prob, 0.5*self.r_prob], # wr x wc
                             [0.25, 0., 0.5, 0., 0., 0.25], # wr x wr
                             [0., 0.5, 0.5, 0., 0., 0.], # wr x cc
                             [0., 0.25, 0.25, 0., 0.25, 0.25], # wr x cr
                             [0., 0., 0.5, 0., 0., 0.5]],# wr x rr
                            [[0., 1., 0., 0., 0., 0.], # cc x ww
                             [0., self.w_prob, 0., self.c_prob, self.r_prob, 0.], # cc x wc
                             [0., 0.5, 0.5, 0., 0., 0.], # cc x wr
                             [0., 0., 0., 1., 0., 0.], # cc x cc
                             [0., 0., 0., 0.5, 0.5, 0.], # cc x cr
                             [0., 0., 0., 0., 1., 0.]],# cc x rr
                            [[0., 0.5, 0.5, 0., 0., 0.], # cr x ww
                             [0., 0.5*self.w_prob, 0.5*self.w_prob, 0.5*self.c_prob, 0.5*(self.c_prob + self.r_prob), 0.5*self.r_prob], # cr x wc
                             [0., 0.25, 0.25, 0., 0.25, 0.25], # cr x wr
                             [0., 0., 0., 0.5, 0.5, 0.], # cr x cc
                             [0., 0., 0., 0.25, 0.5, 0.25], # cr x cr
                             [0., 0., 0., 0., 0.5, 0.5]],# cr x rr
                            [[0., 0., 1., 0., 0., 0.], # rr x ww
                             [0., 0., self.w_prob, 0., self.c_prob, self.r_prob], # rr x wc
                             [0., 0., 0.5, 0., 0., 0.5], # rr x wr
                             [0., 0., 0., 0., 1., 0.], # rr x cc
                             [0., 0., 0., 0., 0.5, 0.5], # rr x cr
                             [0., 0., 0., 0., 0., 1.]]] # rr x rr
        cdef unsigned gt_father, gt_mother, gt_offspring
        for gt_father in range(self.num_genotypes):
            for gt_mother in range(self.num_genotypes):
                for gt_offspring in range(self.num_genotypes):
                    self.inheritance_cube.data.as_floats[gt_father + gt_mother * self.num_genotypes + gt_offspring * self.num_genotypes2] = inheritance_list[gt_father][gt_mother][gt_offspring]

    cdef setFitnessComponents26(self, float h_e, float h_n, float s_e, float s_n):
        if self.num_genotypes != 6:
            raise ValueError("setFitnessComponents26 can only be used if the number of genotypes is 6")
        self.setFitnessComponent(0, 1) # ww
        self.setFitnessComponent(1, (1 - h_e * s_e) * (1 - h_n * s_n)) # wc
        self.setFitnessComponent(2, 1) # wr
        self.setFitnessComponent(3, (1 - s_e) * (1 - s_n)) # cc
        self.setFitnessComponent(4, (1 - h_e * s_e) * (1 - h_n * s_n)) # cr
        self.setFitnessComponent(5, 1) # rr

    cdef setInheritance26(self, float k_c, float k_j, float k_ne):
        if self.num_genotypes != 6:
            raise ValueError("setInheritance26 can only be used if the number of genotypes is 6")
        self.w_prob = 0.5 * (1 - k_c)
        self.c_prob = 0.5 * (1 + k_c * (1 - k_j) * (1 - k_ne))
        self.r_prob = 1 - self.w_prob - self.c_prob
        self.setInheritance()

    cdef float fecundity_proportion(self, unsigned offspring_sex, unsigned mother_gt, unsigned father_gt):
        if (father_gt == 1 or father_gt == 3 or father_gt == 4): # ww, wc, wr, cc, cr, rr.
            if (offspring_sex == 0): # male sex bias if modified father (any mother)
                return self.accuracy
            else:
                return 1. - self.accuracy
        else:
            return 0.5

cdef class CellDynamicsDelayBase(CellDynamicsBase):
    """Base class for lifecycle dynamics as governed by a delay differential equation"""
    
    def __init__(self, delay = 1, current_index = 0):
        """Constructor

        Parameters
        ----------
        delay : unsigned
            number of timesteps involved in the delay, so the total lag = delay * dt (default = 1)
        current_index: unsigned
            defines the generation that have most recently emerged as adults.  0 <= current_index <= delay.  (default = 0)
        """
        
        super().__init__()

        self.delay = delay
        self.current_index = current_index % (delay + 1)

    cdef unsigned getDelay(self):
        return self.delay

    cdef void incrementCurrentIndex(self):
        self.current_index = (self.current_index + 1) % (self.delay + 1)

    cdef unsigned getCurrentIndex(self):
        return self.current_index

cdef class CellDynamics26DelayBase(CellDynamicsDelayBase):
    """Base class for anyODE with 2 sexes, 6 genotypes, using a delay equation"""
    
    def __init__(self, num_species = 3, delay = 1, current_index = 0, death_rate = [[[1.0] * 3] * 6] * 2, competition = [[1, 0, 0], [0, 1, 0], [0, 0, 1]], emergence_rate = [1.0] * 3, activity = [[1, 0, 0], [0, 1, 0], [0, 0, 1]], reduction = [[1.0] * 6] * 6, hybridisation = [[[1, 0, 0], [0, 1, 0], [0, 0, 1]]] * 3, offspring_modifier = [[[1.0] * 3] * 3] * 2, sex_ratio = 0.5, female_bias = 0.5, m_w = 1E-6, m_c = 1E-6):
        """Constructor
        Note that num_sexes = 2 and num_genotypes = 6.  These two parameters could be arguments in the constructor, since all methods use self.num_sexes and self.num_genotypes (ie, no methods hardcode 2 and 6) but no tests exist for different num_sexes and num_genotypes.

        Parameters
        ----------
        num_species : unsigned
            number of mosquito subspecies (default = 3)
        delay : unsigned
            number of timesteps involved in the delay, so the total lag = delay * dt (default = 1)
        current_index: unsigned
            defines the generation that have most recently emerged as adults.  0 <= current_index <= delay.  (default = 0)
        death_rate: list
            death_rate[sex][genotype][mosquito_species].  All elements must be positive (default = 1.0)
        competition: list
            competition[species1][species2].  This is called alpha in the documentation (default = identity)
        emergence_rate: list
            emergence_rate[species].  Larvae per wildtype female at the end of the larval duration period, per day, per adult female, in the absence of density dependence (default = 0.0)
        activity: list
            activity[female_of_species1][male_of_species2] is activity level in the proportionate mixing (default = identity)
        reduction: list
            reduction[gM][gF] is the reduction in adults due to the construct (default = 1.0)
        hybridisation: list
            hybridisation[mM][mF][m] = probability that offspring of species m results from male of species mM and female of species mF (default = delta_{m, mF})
        offspring_modifier: list
            offspring_modifier[s][mM][mF] = suppression (if <1, or increased vigour if >1) of offspring of sex s that arises from male of species mM and female of species mF.
        sex_ratio : float
            probability that offspring of wc or cc or cr fathers are male (paternal male bias has sex_ratio > 0.5) (default = 0.5)
        female_bias : float
            probability that offspring of (mother wc or cc or cr + father ww) is female (usually female_bias > 0.5) (default = 0.5)
        m_w : float
            description.  Default value in report based on Beighton assuming spontaneous resistance (default = 1E-6)
        m_c : float
            description.  Default value in report based on Beighton assuming spontaneous resistance (default = 1E-6)
        """
        
        super().__init__(delay = delay, current_index = current_index)

        self.num_sexes = 2
        self.num_genotypes = 6
        self.num_genotypes2 = 6 * 6

        self.num_species = num_species
        self.num_species2 = num_species * num_species
        self.num_populations = self.num_sexes * self.num_genotypes * self.num_species * (self.delay + 1)

        self.m_w = m_w
        self.m_c = m_c
        
        # Probabilities used in the inheritance cube
        self.w_prob = 0.5 * (1 - self.m_w)
        self.c_prob = 0.5 * (1 - self.m_c)
        self.r_prob = 0.5 * (self.m_w + self.m_c)
        # create the inheritance cube
        self.setInheritance()

        self.setFecundityP(sex_ratio, female_bias)

        self.num_diffusing = self.num_sexes * self.num_genotypes * self.num_species
        self.num_advecting = self.num_sexes * self.num_genotypes * self.num_species
        self.diffusing_indices = array.clone(array.array('I', []), self.num_diffusing, zero = False)
        self.advecting_indices = array.clone(array.array('I', []), self.num_advecting, zero = False)
        self.advection_class = array.array('I', [0] * self.num_advecting)
        cdef unsigned ind = 0
        for sex in range(self.num_sexes):
            for genotype in range(self.num_genotypes):
                for species in range(self.num_species):
                    offset = species + genotype * self.num_species + sex * self.num_species * self.num_genotypes + self.current_index * self.num_species * self.num_genotypes * self.num_sexes
                    self.diffusing_indices.data.as_uints[ind] = offset
                    self.advecting_indices.data.as_uints[ind] = offset
                    ind = ind + 1

        self.setDeathRate(death_rate)
        self.setCompetition(competition)
        self.setEmergenceRate(emergence_rate)
        self.setActivity(activity)
        self.setReduction(reduction)
        self.setHybridisation(hybridisation)
        self.setOffspringModifier(offspring_modifier)

        self.have_precalculated = 0

    cdef void incrementCurrentIndex(self):
        # don't know why i can't super().incrementCurrentIndex()
        self.current_index = (self.current_index + 1) % (self.delay + 1)
        cdef unsigned sex, genotype, species, offset
        cdef unsigned ind = 0
        for sex in range(self.num_sexes):
            for genotype in range(self.num_genotypes):
                for species in range(self.num_species):
                    offset = species + genotype * self.num_species + sex * self.num_species * self.num_genotypes + self.current_index * self.num_species * self.num_genotypes * self.num_sexes
                    self.diffusing_indices.data.as_uints[ind] = offset
                    self.advecting_indices.data.as_uints[ind] = offset
                    ind = ind + 1

    cdef setDeathRate(self, list death_rate):
        self.death_rate = array.clone(array.array('f', []), self.num_sexes * self.num_genotypes * self.num_species, zero = False)
        if len(death_rate) != self.num_sexes:
            raise ValueError("size of death_rate, " + str(len(death_rate)) + ", must be equal to " + str(self.num_sexes))
        for s in range(self.num_sexes):
            if len(death_rate[s]) != self.num_genotypes:
                raise ValueError("size of death_rate[" + str(s) + "], " + str(len(death_rate[s])) + ", must be equal to " + str(self.num_genotypes))
            for g in range(self.num_genotypes):
                if len(death_rate[s][g]) != self.num_species:
                    raise ValueError("size of death_rate[" + str(s) + "][" + str(g) + "], " + str(len(death_rate[s][g])) + ", must be equal to " + str(self.num_species))
                for m in range(self.num_species):
                    if death_rate[s][g][m] <= 0.0:
                        raise ValueError("all death rates must be positive")
                    self.death_rate[m + g * self.num_species + s * self.num_species * self.num_genotypes] = death_rate[s][g][m]

    cdef list getDeathRate(self):
        return [[[self.death_rate[m + g * self.num_species + s * self.num_species * self.num_genotypes] for m in range(self.num_species)] for g in range(self.num_genotypes)] for s in range(self.num_sexes)]
                     
    cdef setCompetition(self, list competition):
        self.competition = array.clone(array.array('f', []), self.num_species * self.num_species, zero = False)
        if len(competition) != self.num_species:
            raise ValueError("size of competition, " + str(len(competition)) + ", must be equal to " + str(self.num_species))
        for m in range(self.num_species):
            if len(competition[m]) != self.num_species:
                raise ValueError("size of competition[" + str(m) + "], " + str(len(competition[m])) + ", must be equal to " + str(self.num_species))
            for mprime in range(self.num_species):
                self.competition[mprime + m * self.num_species] = competition[m][mprime]
            
    cdef list getCompetition(self):
        return [[self.competition[mprime + m * self.num_species] for mprime in range(self.num_species)] for m in range(self.num_species)]

    cdef setEmergenceRate(self, list emergence_rate):
        if len(emergence_rate) != self.num_species:
            raise ValueError("size of emergence_rate, " + str(len(emergence_rate)) + ", must be equal to " + str(self.num_species))
        for er in emergence_rate:
            if er < 0.0:
                raise ValueError("all emergence rates must be non-negative")
        self.emergence_rate = array.array('f', emergence_rate)
        self.have_precalculated = 0

    cdef list getEmergenceRate(self):
        return [self.emergence_rate[i] for i in range(self.num_species)]
                     
    cdef setActivity(self, list activity):
        self.activity = array.clone(array.array('f', []), self.num_species * self.num_species, zero = False)
        if len(activity) != self.num_species:
            raise ValueError("size of activity, " + str(len(activity)) + ", must be equal to " + str(self.num_species))
        for mF in range(self.num_species):
            if len(activity[mF]) != self.num_species:
                raise ValueError("size of activity[" + str(mF) + "], " + str(len(activity[mF])) + ", must be equal to " + str(self.num_species))
            for mM in range(self.num_species):
                if activity[mF][mM] < 0.0:
                    raise ValueError("all activity values must be non-negative")
                self.activity[mM + mF * self.num_species] = activity[mF][mM]

    cdef list getActivity(self):
        return [[self.activity[mM + mF * self.num_species] for mM in range(self.num_species)] for mF in range(self.num_species)]
                     
    cdef void evolveTrial(self, float timestep, float[:] pops_and_params):
        """This just implements dx/dt = -death_rate * x + lambdah * x[t - delay * dt], which
        discretises to x[t + dt] = x[t - delay * dt] / death_rate + (x[t] - x[t - delay * dt] / death_rate) * exp(-death_rate * dt)
        This function is not optimised!"""
        
        cdef float lambdah = 1.1
        cdef array.array new_pop = array.clone(array.array('f', []), self.num_sexes * self.num_genotypes * self.num_species, zero = False)

        cdef unsigned adult_base = self.current_index * self.num_species * self.num_genotypes * self.num_sexes
        cdef unsigned delayed_base = (self.current_index + 1) % (self.delay + 1) * self.num_species * self.num_genotypes * self.num_sexes

        cdef unsigned current_index = 0
        cdef unsigned delayed_index = 0
        cdef unsigned ind = 0
        cdef float dr = 0.0
        for sex in range(self.num_sexes):
            for genotype in range(self.num_genotypes):
                for species in range(self.num_species):
                    current_index = adult_base + ind
                    delayed_index = delayed_base + ind
                    dr = self.death_rate[ind]
                    new_pop[ind] = lambdah * pops_and_params[delayed_index] / dr + (pops_and_params[current_index] - lambdah * pops_and_params[delayed_index] / dr) * exp(- dr * timestep)
                    ind += 1

        ind = 0
        for sex in range(self.num_sexes):
            for genotype in range(self.num_genotypes):
                for species in range(self.num_species):
                    delayed_index = delayed_base + ind
                    pops_and_params[delayed_index] = new_pop[ind]
                    ind += 1

    cdef void setInheritance(self):
        self.inheritance_cube = array.clone(array.array('f', []), self.num_genotypes2 * self.num_genotypes, zero = False)

        # inherited allele frequencies based on parental genotype
        allele_list = [[2.*self.w_prob, 0., (1 - 2.*self.w_prob)], # allele freqs from ww parent
        [self.w_prob, self.c_prob, self.r_prob], # wc
        [self.w_prob, 0., (1. - self.w_prob)], # wr
        [0., 2. * self.c_prob, (1. - 2. * self.c_prob)], # cc
        [0., self.c_prob, (1 - self.c_prob)], # cr
        [0., 0., 1.]] # rr
        
        # Mendelian inheritance of parental alleles
        cdef unsigned gt_father, gt_mother, gt_offspring
        for gt_father in range(self.num_genotypes):
            for gt_mother in range(self.num_genotypes):
                self.inheritance_cube.data.as_floats[gt_father + gt_mother * self.num_genotypes + 0 * self.num_genotypes2] = allele_list[gt_father][0] * allele_list[gt_mother][0] # ww offspring
                self.inheritance_cube.data.as_floats[gt_father + gt_mother * self.num_genotypes + 1 * self.num_genotypes2] = allele_list[gt_father][0] * allele_list[gt_mother][1] +  allele_list[gt_father][1] * allele_list[gt_mother][0] # wc offspring (w from mother, c from father or vice versa)               
                self.inheritance_cube.data.as_floats[gt_father + gt_mother * self.num_genotypes + 2 * self.num_genotypes2] = allele_list[gt_father][0] * allele_list[gt_mother][2] + allele_list[gt_father][2] * allele_list[gt_mother][0] # wr offspring                
                self.inheritance_cube.data.as_floats[gt_father + gt_mother * self.num_genotypes + 3 * self.num_genotypes2] = allele_list[gt_father][1] * allele_list[gt_mother][1] # cc offspring
                self.inheritance_cube.data.as_floats[gt_father + gt_mother * self.num_genotypes + 4 * self.num_genotypes2] = allele_list[gt_father][1] * allele_list[gt_mother][2] + allele_list[gt_father][2] * allele_list[gt_mother][1] # cr offspring
                self.inheritance_cube.data.as_floats[gt_father + gt_mother * self.num_genotypes + 5 * self.num_genotypes2] = allele_list[gt_father][2] * allele_list[gt_mother][2] # rr offspring
        self.have_precalculated = 0

    cdef array.array getInheritance(self):
        return self.inheritance_cube

    cdef void setFecundityP(self, float sex_ratio, float female_bias):
        self.sex_ratio = sex_ratio
        self.female_bias = female_bias
        self.fecundity_p = array.array('f', [0.5] * self.num_sexes * self.num_genotypes2)
        cdef unsigned gM, gF, s
        cdef unsigned ind

        for gM in [1, 3, 4]: # wc, cc, cr
            for gF in range(self.num_genotypes):
                s = 0 # Male
                self.fecundity_p[gM + gF * self.num_genotypes + s * self.num_genotypes2] = self.sex_ratio
                s = 1 # Female
                self.fecundity_p[gM + gF * self.num_genotypes + s * self.num_genotypes2] = 1.0 - self.sex_ratio

        for gM in [0]: # ww
            for gF in [1, 3, 4]: # wc, cc, cr
                s = 0 # Male
                self.fecundity_p[gM + gF * self.num_genotypes + s * self.num_genotypes2] = 1.0 - self.female_bias
                s = 1 # Female
                self.fecundity_p[gM + gF * self.num_genotypes + s * self.num_genotypes2] = self.female_bias

        self.have_precalculated = 0

    cdef float getSexRatio(self):
        return self.sex_ratio

    cdef float getFemaleBias(self):
        return self.female_bias

    cdef list getFecundityP(self):
        return [[[self.fecundity_p[gM + gF * self.num_genotypes + s * self.num_genotypes * self.num_genotypes] for s in range(self.num_sexes)] for gF in range(self.num_genotypes)] for gM in range(self.num_genotypes)]

    cdef setReduction(self, list reduction):
        self.reduction = array.clone(array.array('f', []), self.num_genotypes2, zero = False)
        if len(reduction) != self.num_genotypes:
            raise ValueError("size of reduction, " + str(len(reduction)) + ", must be equal to " + str(self.num_genotypes))
        for gM in range(self.num_genotypes):
            if len(reduction[gM]) != self.num_genotypes:
                raise ValueError("size of reduction[" + str(gM) + "], " + str(len(reduction[gM])) + ", must be equal to " + str(self.num_genotypes))
            for gF in range(self.num_genotypes):
                self.reduction[gF + gM * self.num_genotypes] = reduction[gM][gF]
        self.have_precalculated = 0
        

    cdef list getReduction(self):
        return [[self.reduction[gF + gM * self.num_genotypes] for gF in range(self.num_genotypes)] for gM in range(self.num_genotypes)]

    cdef setHybridisation(self, list hybridisation):
        self.hybridisation = array.clone(array.array('f', []), self.num_species2 * self.num_species, zero = False)
        if len(hybridisation) != self.num_species:
            raise ValueError("size of hybridisation, " + str(len(hybridisation)) + ", must be equal to " + str(self.num_species))
        for mM in range(self.num_species):
            if len(hybridisation[mM]) != self.num_species:
                raise ValueError("size of hybridisation[" + str(mM) + "], " + str(len(hybridisation[mM])) + ", must be equal to " + str(self.num_species))
            for mF in range(self.num_species):
                if len(hybridisation[mM][mF]) != self.num_species:
                    raise ValueError("size of hybridisation[" + str(mM) + "][" + str(mF) + "], " + str(len(hybridisation[mM][mF])) + ", must be equal to " + str(self.num_species))
                for m in range(self.num_species):
                    self.hybridisation[m + mF * self.num_species + mM * self.num_species2] = hybridisation[mM][mF][m]
        self.have_precalculated = 0

    cdef list getHybridisation(self):
        return [[[self.hybridisation[m + mF * self.num_species + mM * self.num_species * self.num_species] for m in range(self.num_species)] for mF in range(self.num_species)] for mM in range(self.num_species)]

    cdef setOffspringModifier(self, list offspring_modifier):
        self.offspring_modifier = array.clone(array.array('f', []), self.num_species2 * self.num_sexes, zero = False)
        if len(offspring_modifier) != self.num_sexes:
            raise ValueError("size of offspring_modifier, " + str(len(offspring_modifier)) + ", must be equal to " + str(self.num_sexes))
        for s in range(self.num_sexes):
            if len(offspring_modifier[s]) != self.num_species:
                raise ValueError("size of offspring_modifier[" + str(s) + "], " + str(len(offspring_modifier[s])) + ", must be equal to " + str(self.num_species))
            for mM in range(self.num_species):
                if len(offspring_modifier[s][mM]) != self.num_species:
                    raise ValueError("size of offspring_modifier[" + str(s) + "][" + str(mM) + "], " + str(len(offspring_modifier[s][mM])) + ", must be equal to " + str(self.num_species))
                for mF in range(self.num_species):
                    self.offspring_modifier[mF + mM * self.num_species + s * self.num_species2] = offspring_modifier[s][mM][mF]
        self.have_precalculated = 0

    cdef list getOffspringModifier(self):
        return [[[self.offspring_modifier[mF + mM * self.num_species + s * self.num_species2] for mF in range(self.num_species)] for mM in range(self.num_species)] for s in range(self.num_sexes)]

    cdef void setMinCarryingCapacity(self, float value):
        self.min_cc = value

    cdef float getMinCarryingCapacity(self):
        return self.min_cc

    cdef precalculate(self):
        self.have_precalculated = 1

    cdef unsigned getNumSpecies(self):
        return self.num_species
    
cdef class CellDynamicsMosquitoLogistic26Delay(CellDynamics26DelayBase):
    """Mosquito lifecycle dynamics as governed by a delay differential equation using logistic growth"""
    
    def __init__(self, num_species = 3, delay = 1, current_index = 0, death_rate = [[[1.0] * 3] * 6] * 2, competition = [[1, 0, 0], [0, 1, 0], [0, 0, 1]], emergence_rate = [1.0] * 3, activity = [[1, 0, 0], [0, 1, 0], [0, 0, 1]], reduction = [[1.0] * 6] * 6, hybridisation = [[[1, 0, 0], [0, 1, 0], [0, 0, 1]]] * 3, offspring_modifier = [[[1.0] * 3] * 3] * 2, sex_ratio = 0.5, female_bias = 0.5, m_w = 1E-6, m_c = 1E-6, min_cc = 1E-6):
        """Constructor
        Note that num_sexes = 2 and num_genotypes = 6.  These two parameters could be arguments in the constructor, since all methods use self.num_sexes and self.num_genotypes (ie, no methods hardcode 2 and 6) but no tests exist for different num_sexes and num_genotypes.

        Parameters
        ----------
        num_species : unsigned
            number of mosquito subspecies (default = 3)
        delay : unsigned
            number of timesteps involved in the delay, so the total lag = delay * dt (default = 1)
        current_index: unsigned
            defines the generation that have most recently emerged as adults.  0 <= current_index <= delay.  (default = 0)
        death_rate: list
            death_rate[sex][genotype][mosquito_species].  All elements must be positive (default = 1.0)
        competition: list
            competition[species1][species2].  This is called alpha in the documentation (default = identity)
        emergence_rate: list
            emergence_rate[species].  Larvae per wildtype female at the end of the larval duration period, per day, per adult female, in the absence of density dependence (default = 0.0)
        activity: list
            activity[female_of_species1][male_of_species2] is activity level in the proportionate mixing (default = identity)
        reduction: list
            reduction[gM][gF] is the reduction in adults due to the construct (default = 1.0)
        hybridisation: list
            hybridisation[mM][mF][m] = probability that offspring of species m results from male of species mM and female of species mF (default = delta_{m, mF})
        offspring_modifier: list
            offspring_modifier[s][mM][mF] = suppression (if <1, or increased vigour if >1) of offspring of sex s that arises from male of species mM and female of species mF.
        sex_ratio : float
            probability that offspring of wc or cc fathers are male (paternal male bias has sex_ratio > 0.5) (default = 0.5)
        female_bias : float
            probability that offspring of (mother wc or cc + father ww) is female (usually female_bias > 0.5) (default = 0.5)
        m_w : float
            description.  Default value in report based on Beighton assuming spontaneous resistance (default = 1E-6)
        m_c : float
            description.  Default value in report based on Beighton assuming spontaneous resistance (default = 1E-6)
        min_cc : float
            minimum carrying capacity.  If carrying capacity in pops_and_params is less than this quantity, no births will occur
        """
        
        super().__init__(num_species = num_species, delay = delay, current_index = current_index, death_rate = death_rate, competition = competition, emergence_rate = emergence_rate, activity = activity, reduction = reduction, hybridisation = hybridisation, offspring_modifier = offspring_modifier, sex_ratio = sex_ratio, female_bias = female_bias, m_w = m_w, m_c = m_c)

        self.num_parameters = self.num_species # carrying capacities
        self.new_pop = array.clone(array.array('f', []), self.num_sexes * self.num_genotypes * self.num_species, zero = False)
        self.xprimeM = array.clone(array.array('f', []), self.num_species * self.num_genotypes * self.num_species, zero = False)
        self.yy = array.clone(array.array('f', []), self.num_sexes * self.num_genotypes * self.num_species, zero = False)
        self.yyp = array.clone(array.array('f', []), self.num_sexes * self.num_genotypes * self.num_species, zero = False)
        self.precalc = array.clone(array.array('f', []), self.num_sexes * self.num_genotypes2 * self.num_genotypes * self.num_species2 * self.num_species, zero = False)
        self.precalcp = array.clone(array.array('f', []), self.num_sexes * self.num_genotypes2 * self.num_genotypes * self.num_species2 * self.num_species, zero = False)
        self.comp = array.clone(array.array('f', []), self.num_species, zero = False)
        self.min_cc = min_cc
        self.precalculate()

    cdef void evolve(self, float timestep, float[:] pops_and_params):
        """This function is not optimised"""

        if self.have_precalculated == 0:
            self.precalculate()
        
        cdef unsigned mF, mM, gM, gF, mprime, gprime, sex, ind, s, g, m
        cdef unsigned current_index, delayed_index, yy_ind, f_ind, precalc_ind
        cdef float denom, xF, bb

        cdef unsigned adult_base = self.current_index * self.num_species * self.num_genotypes * self.num_sexes
        cdef unsigned delayed_base = (self.current_index + 1) % (self.delay + 1) * self.num_species * self.num_genotypes * self.num_sexes

        for mF in range(self.num_species):
            denom = 0.0
            sex = 0 # male
            for gprime in range(self.num_genotypes):
                for mprime in range(self.num_species):
                    delayed_index = mprime + gprime * self.num_species + delayed_base # + sex * num_species * num_genotypes
                    denom += self.activity[mprime + mF * self.num_species] * pops_and_params[delayed_index]
            for gM in range(self.num_genotypes):
                for mM in range(self.num_species):
                    delayed_index = mM + gM * self.num_species + delayed_base # + sex * num_species * num_genotypes
                    ind = mF + mM * self.num_species + gM * self.num_species2
                    if denom <= 0.0:
                        self.xprimeM[ind] = 0.0
                    else:
                        self.xprimeM[ind] = self.activity[mM + mF * self.num_species] * pops_and_params[delayed_index] / denom

        array.zero(self.yy)
        array.zero(self.yyp)

        for s in range(self.num_sexes):
            for g in range(self.num_genotypes):
                for m in range(self.num_species):
                    yy_ind = m + g * self.num_species + s * self.num_species * self.num_genotypes
                    for mF in range(self.num_species):
                        for gF in range(self.num_genotypes):
                            f_ind = mF + gF * self.num_species + 1 * self.num_species * self.num_genotypes + delayed_base
                            xF = pops_and_params[f_ind]
                            for mM in range(self.num_species):
                                for gM in range(self.num_genotypes):
                                    xprime_ind = mF + mM * self.num_species + gM * self.num_species2
                                    precalc_ind = s + self.num_sexes * (g + self.num_genotypes * (gF + self.num_genotypes * (gM + self.num_genotypes * (m + self.num_species * (mF + self.num_species * mM)))))
                                    self.yy[yy_ind] += self.precalc[precalc_ind] * xF * self.xprimeM[xprime_ind]
                                    self.yyp[yy_ind] += self.precalcp[precalc_ind] * xF * self.xprimeM[xprime_ind]

        array.zero(self.comp)
        for m in range(self.num_species):
            for mprime in range(self.num_species):
                alpha_ind = mprime + m * self.num_species
                alpha = self.competition[alpha_ind]
                for s in range(self.num_sexes):
                    for g in range(self.num_genotypes):
                        yy_ind = m + g * self.num_species + s * self.num_species * self.num_genotypes
                        self.comp[m] += alpha * self.yy[yy_ind]
        # calculate max(0, 1 - comp/carrying_capacity) and shove into self.comp[m]
        for m in range(self.num_species):
            if pops_and_params[self.num_populations + m] < self.min_cc:
                self.comp[m] = 0.0
            else:
                self.comp[m] = max(0.0, 1.0 - self.comp[m] / pops_and_params[self.num_populations + m])

        for g in range(self.num_genotypes):
            for m in range(self.num_species):
                for s in range(self.num_sexes):
                    ind = m + g * self.num_species + s * self.num_species * self.num_genotypes
                    dr = self.death_rate[ind]
                    bb = self.comp[m] * self.yyp[ind]
                    current_index = adult_base + ind
                    self.new_pop[ind] = bb / dr + (pops_and_params[current_index] - bb / dr) * exp(- dr * timestep)

        # put the new_pop in the correct slots in pops_and_params in readyness for incrementCurrentIndex
        for g in range(self.num_genotypes):
            for m in range(self.num_species):
                for s in range(self.num_sexes):
                    ind = m + g * self.num_species + s * self.num_species * self.num_genotypes
                    delayed_ind = delayed_base + ind
                    pops_and_params[delayed_ind] = self.new_pop[ind]

    cdef array.array calcQm(self, float[:] eqm_pops_and_params):
        """Not defined for this class"""
        return

    cdef void setMinCarryingCapacity(self, float value):
        self.min_cc = value

    cdef float getMinCarryingCapacity(self):
        return self.min_cc

    cdef precalculate(self):
        cdef unsigned mM, mF, m, gM, gF, g, s, ind
        array.zero(self.precalc)
        array.zero(self.precalcp)
        for mM in range(self.num_species):
            for mF in range(self.num_species):
                for m in range(self.num_species):
                    for gM in range(self.num_genotypes):
                        for gF in range(self.num_genotypes):
                            for g in range(self.num_genotypes):
                                for s in range(self.num_sexes):
                                    ind = s + self.num_sexes * (g + self.num_genotypes * (gF + self.num_genotypes * (gM + self.num_genotypes * (m + self.num_species * (mF + self.num_species * mM)))))
                                    self.precalc[ind] += self.hybridisation[m + mF * self.num_species + mM * self.num_species2] * self.emergence_rate[mF] * self.inheritance_cube[gM + gF * self.num_genotypes + g * self.num_genotypes2] * self.fecundity_p[gM + gF * self.num_genotypes + s * self.num_genotypes2] * self.reduction[gF + gM * self.num_genotypes]
                                    self.precalcp[ind] += self.offspring_modifier[mF + mM * self.num_species + s * self.num_species2] * self.hybridisation[m + mF * self.num_species + mM * self.num_species2] * self.emergence_rate[mF] * self.inheritance_cube[gM + gF * self.num_genotypes + g * self.num_genotypes2] * self.fecundity_p[gM + gF * self.num_genotypes + s * self.num_genotypes2] * self.reduction[gF + gM * self.num_genotypes]
        self.have_precalculated = 1
        
cdef class CellDynamicsMosquitoBH26Delay(CellDynamics26DelayBase):
    """Mosquito lifecycle dynamics as governed by a delay differential equation using Beverton-Holt growth"""
    
    def __init__(self, num_species = 3, delay = 1, current_index = 0, death_rate = [[[1.0] * 3] * 6] * 2, competition = [[1, 0, 0], [0, 1, 0], [0, 0, 1]], emergence_rate = [1.0] * 3, activity = [[1, 0, 0], [0, 1, 0], [0, 0, 1]], reduction = [[1.0] * 6] * 6, hybridisation = [[[1, 0, 0], [0, 1, 0], [0, 0, 1]]] * 3, offspring_modifier = [[[1.0] * 3] * 3] * 2, sex_ratio = 0.5, female_bias = 0.5, m_w = 1E-6, m_c = 1E-6):
        """Constructor
        Note that num_sexes = 2 and num_genotypes = 6.  These two parameters could be arguments in the constructor, since all methods use self.num_sexes and self.num_genotypes (ie, no methods hardcode 2 and 6) but no tests exist for different num_sexes and num_genotypes.

        Parameters
        ----------
        num_species : unsigned
            number of mosquito subspecies (default = 3)
        delay : unsigned
            number of timesteps involved in the delay, so the total lag = delay * dt (default = 1)
        current_index: unsigned
            defines the generation that have most recently emerged as adults.  0 <= current_index <= delay.  (default = 0)
        death_rate: list
            death_rate[sex][genotype][mosquito_species].  All elements must be positive (default = 1.0)
        competition: list
            competition[species1][species2].  This is called alpha in the documentation (default = identity)
        emergence_rate: list
            emergence_rate[species].  Larvae per wildtype female at the end of the larval duration period, per day, per adult female, in the absence of density dependence (default = 0.0)
        activity: list
            activity[female_of_species1][male_of_species2] is activity level in the proportionate mixing (default = identity)
        reduction: list
            reduction[gM][gF] is the reduction in adults due to the construct (default = 1.0)
        hybridisation: list
            hybridisation[mM][mF][m] = probability that offspring of species m results from male of species mM and female of species mF (default = delta_{m, mF})
        offspring_modifier: list
            offspring_modifier[s][mM][mF] = suppression (if <1, or increased vigour if >1) of offspring of sex s that arises from male of species mM and female of species mF.
        sex_ratio : float
            probability that offspring of wc or cc fathers are male (paternal male bias has sex_ratio > 0.5) (default = 0.5)
        female_bias : float
            probability that offspring of (mother wc or cc + father ww) is female (usually female_bias > 0.5) (default = 0.5)
        m_w : float
            description.  Default value in report based on Beighton assuming spontaneous resistance (default = 1E-6)
        m_c : float
            description.  Default value in report based on Beighton assuming spontaneous resistance (default = 1E-6)
        """
        
        super().__init__(num_species = num_species, delay = delay, current_index = current_index, death_rate = death_rate, competition = competition, emergence_rate = emergence_rate, activity = activity, reduction = reduction, hybridisation = hybridisation, offspring_modifier = offspring_modifier, sex_ratio = sex_ratio, female_bias = female_bias, m_w = m_w, m_c = m_c)

        self.num_parameters = self.num_species # "q" parameters
        self.new_pop = array.clone(array.array('f', []), self.num_sexes * self.num_genotypes * self.num_species, zero = False)
        self.xprimeM = array.clone(array.array('f', []), self.num_species * self.num_genotypes * self.num_species, zero = False)
        self.yy = array.clone(array.array('f', []), self.num_sexes * self.num_genotypes * self.num_species, zero = False)
        self.yyp = array.clone(array.array('f', []), self.num_sexes * self.num_genotypes * self.num_species, zero = False)
        self.precalc = array.clone(array.array('f', []), self.num_sexes * self.num_genotypes2 * self.num_genotypes * self.num_species2 * self.num_species, zero = False)
        self.precalcp = array.clone(array.array('f', []), self.num_sexes * self.num_genotypes2 * self.num_genotypes * self.num_species2 * self.num_species, zero = False)
        self.eqm_pops_and_params = array.clone(array.array('f', []), self.num_populations + self.num_parameters, zero = False)
        self.comp = array.clone(array.array('f', []), self.num_species, zero = False)
        self.death_terms = array.clone(array.array('f', []), self.num_species, zero = False)
        self.qm_vals = array.clone(array.array('f', []), self.num_species, zero = False)
        self.num_genotypes_to_calc = self.num_genotypes
        self.num_sexes_to_calc = self.num_sexes
        self.use_qm = 1
        self.species_present = array.clone(array.array('B', []), self.num_species, zero = False) 
        self.genotype_present = array.clone(array.array('B', []), self.num_genotypes, zero = False)
        self.precalculate()

    cdef setNumGenotypesToCalc(self, unsigned num_genotypes_to_calc):
        if num_genotypes_to_calc > self.num_genotypes:
            raise ValueError("setNumGenotypesToCalc: num_genotypes_to_calc must be <= " + str(self.num_genotypes))
        self.num_genotypes_to_calc = num_genotypes_to_calc

    cdef unsigned getNumGenotypesToCalc(self):
        return self.num_genotypes_to_calc

    cdef setNumSexesToCalc(self, unsigned num_sexes_to_calc):
        if num_sexes_to_calc > self.num_sexes:
            raise ValueError("setNumSexesToCalc: num_sexes_to_calc must be <= " + str(self.num_sexes))
        self.num_sexes_to_calc = num_sexes_to_calc

    cdef unsigned getNumSexesToCalc(self):
        return self.num_sexes_to_calc

    cdef setUseQm(self, unsigned use_qm):
        if not (use_qm == 0 or use_qm == 1):
            raise ValueError("setUseQm: use_qm must be 0 or 1")
        self.use_qm = use_qm

    cdef unsigned getUseQm(self):
        return self.use_qm

    cdef void evolve(self, float timestep, float[:] pops_and_params):

        if self.have_precalculated == 0:
            self.precalculate()
        
        cdef unsigned ind, s, g, m, current_index, delayed_index, some_males
        cdef float bb, one_over_dr, expdrdt, qm

        cdef unsigned adult_base = self.current_index * self.num_species * self.num_genotypes * self.num_sexes
        cdef unsigned delayed_base = (self.current_index + 1) % (self.delay + 1) * self.num_species * self.num_genotypes * self.num_sexes
        

        if self.use_qm == 0:
            # populate eqm_pops_and_params with carrying capacity in relevant places (others don't matter as are not used)
            g = 0 # wild-type only
            for m in range(self.num_species):
                for s in range(self.num_sexes):
                    ind = m + g * self.num_species + s * self.num_species * self.num_genotypes
                    current_index = adult_base + ind
                    self.eqm_pops_and_params[current_index] = 0.5 * pops_and_params[self.num_populations + m] # 0.5 comes from half males and half females

            self.qm_vals = self.calcQm(self.eqm_pops_and_params)

        # Calculate indicators whether a given species or genotype is present
        # (based on delayed data, as this is where births come from)
        # (The above calcQm() function call might have set species_present = 1 (based on eqm_pops_and_param) and genotype_present[0] = 1.  That is probably irrelevant in most calculations, but worth knowing.)
        array.zero(self.species_present)
        for m in range(self.num_species):
            for g in range(self.num_genotypes):
                for s in range(self.num_sexes):
                    ind = m + g * self.num_species + s * self.num_species * self.num_genotypes
                    delayed_ind = delayed_base + ind
                    if pops_and_params[delayed_ind] > 0:
                        self.species_present.data.as_uchars[m] = 1
                        break
                if pops_and_params[delayed_ind] > 0:
                    break

        array.zero(self.genotype_present)
        for g in range(self.num_genotypes):
            for m in range(self.num_species):
                for s in range(self.num_sexes):
                    ind = m + g * self.num_species + s * self.num_species * self.num_genotypes
                    delayed_ind = delayed_base + ind
                    if pops_and_params[delayed_ind] > 0:
                        self.genotype_present.data.as_uchars[g] = 1
                        break
                if pops_and_params[delayed_ind] > 0:
                    break

        # calculate xprimeM
        some_males = self.calcXprimeM(delayed_base, pops_and_params)  # zero if there are no delayed males whatsoever

        # calculate Y
        self.calcYYprime(some_males, delayed_base, pops_and_params)

        # calculate C
        self.calcCompetition(some_males)

        # calculate birth rate and new populations
        for g in range(self.num_genotypes):
            for m in range(self.num_species):
                if self.use_qm == 1:
                    ind = self.num_populations + m
                    qm = pops_and_params[ind]
                else:
                    qm = self.qm_vals.data.as_floats[m]
                if qm <= self.small_value:
                    qm = 0.0 # this accounts for problems in double-to-float conversion, eg, if qm = 1E-40
                for s in range(self.num_sexes):
                    ind = m + g * self.num_species + s * self.num_species * self.num_genotypes
                    # note that death_rate > 0 by the constructor and setter
                    # unoptimised:
                    #one_over_dr = self.death_rate[ind]
                    #expdrdt = exp(- self.death_rate[ind] * timestep)
                    # better:
                    one_over_dr = 1.0 / self.death_rate.data.as_floats[ind]
                    expdrdt = exp(- self.death_rate.data.as_floats[ind] * timestep)
                    # unoptimised:
                    #if self.yyp[ind] <= 0.0 or qm <= self.small.value:
                    #    bb = 0
                    #else:
                    #    bb = qm * self.yyp[ind] / (qm + self.comp[m])
                    #current_index = adult_base + ind
                    #self.new_pop[ind] = bb * one_over_dr + (pops_and_params[current_index] - bb * one_over_dr) * expdrdt
                    # better:
                    if self.yyp.data.as_floats[ind] <= 0.0 or qm <= self.small_value:
                        bb = 0
                    else:
                        bb = qm * self.yyp.data.as_floats[ind] / (qm + self.comp.data.as_floats[m])
                    current_index = adult_base + ind
                    self.new_pop.data.as_floats[ind] = bb * one_over_dr + (pops_and_params[current_index] - bb * one_over_dr) * expdrdt
                    if self.new_pop.data.as_floats[ind] <= self.small_value:
                        self.new_pop.data.as_floats[ind] = 0.0

        # put the new_pop in the correct slots in pops_and_params in readiness for incrementCurrentIndex
        for s in range(self.num_sexes):
            for g in range(self.num_genotypes):
                for m in range(self.num_species):
                    ind = m + g * self.num_species + s * self.num_species * self.num_genotypes
                    delayed_ind = delayed_base + ind
                    # unoptimised:
                    #pops_and_params[delayed_ind] = self.new_pop[ind]
                    # better:
                    pops_and_params[delayed_ind] = self.new_pop.data.as_floats[ind]

    cdef precalculate(self):
        cdef unsigned mM, mF, m, gM, gF, g, s, ind
        array.zero(self.precalc)
        array.zero(self.precalcp)
        for mM in range(self.num_species):
            for mF in range(self.num_species):
                for m in range(self.num_species):
                    for gM in range(self.num_genotypes):
                        for gF in range(self.num_genotypes):
                            for g in range(self.num_genotypes):
                                for s in range(self.num_sexes):
                                    ind = s + self.num_sexes * (g + self.num_genotypes * (gF + self.num_genotypes * (gM + self.num_genotypes * (m + self.num_species * (mF + self.num_species * mM)))))
                                    self.precalc[ind] += self.hybridisation[m + mF * self.num_species + mM * self.num_species2] * self.emergence_rate[mF] * self.inheritance_cube[gM + gF * self.num_genotypes + g * self.num_genotypes2] * self.fecundity_p[gM + gF * self.num_genotypes + s * self.num_genotypes2] * self.reduction[gF + gM * self.num_genotypes]
                                    self.precalcp[ind] += self.offspring_modifier[mF + mM * self.num_species + s * self.num_species2] * self.hybridisation[m + mF * self.num_species + mM * self.num_species2] * self.emergence_rate[mF] * self.inheritance_cube[gM + gF * self.num_genotypes + g * self.num_genotypes2] * self.fecundity_p[gM + gF * self.num_genotypes + s * self.num_genotypes2] * self.reduction[gF + gM * self.num_genotypes]
        self.have_precalculated = 1
        
    cdef array.array calcQm(self, float[:] eqm_pops_and_params):
        if self.have_precalculated == 0:
            self.precalculate()

        cdef unsigned m, g, s, some_males
        cdef unsigned adult_base = self.current_index * self.num_species * self.num_genotypes * self.num_sexes

        self.setNumGenotypesToCalc(1)
        self.setNumSexesToCalc(1)

        # NOTE: the following lines potentially alter genotype_present and species_present
        self.genotype_present.data.as_uchars[0] = 1 # this function assumes there are nonzero genotype=0 populations
        for m in range(self.num_species):
            if eqm_pops_and_params[m + adult_base] > 0 or eqm_pops_and_params[m + 1 * self.num_species * self.num_genotypes + adult_base] > 0:
                self.species_present.data.as_uchars[m] = 1

        # calculate xprimeM
        some_males = self.calcXprimeM(adult_base, eqm_pops_and_params)  # zero if there are no adult males whatsoever

        # calculate Y
        self.calcYYprime(some_males, adult_base, eqm_pops_and_params)

        # calculate C
        self.calcCompetition(some_males)

        self.setNumGenotypesToCalc(self.num_genotypes)
        self.setNumSexesToCalc(self.num_sexes)
        # calculate sum over genotypes and sexes of the death terms and the Y' terms
        array.zero(self.death_terms)
        for m in range(self.num_species):
            self.death_terms.data.as_floats[m] = self.death_rate.data.as_floats[m] * eqm_pops_and_params[m + adult_base]

        for m in range(self.num_species):
            if self.death_terms.data.as_floats[m] <= 0.0:
                raise ValueError("Sum_{g, s}d_{g, s, m}X_{g, s, m} for m = " + str(m) + " is zero")
        # calculate the resulting qm values
        for m in range(self.num_species):
            self.qm_vals.data.as_floats[m] = self.comp.data.as_floats[m] / (self.yyp.data.as_floats[m] / self.death_terms.data.as_floats[m] - 1.0)

        return self.qm_vals

    cdef unsigned calcXprimeM(self, unsigned delayed_base, float[:] pops_and_params):
        cdef unsigned some_males = 0 # zero if there are no delayed males whatsoever
        cdef unsigned mF, sex, gprime, mprime, delayed_index, ind, gM, mM, inda
        cdef float denom, one_over_denom
        for mF in range(self.num_species):
            denom = 0.0
            sex = 0 # male
            for gprime in range(self.num_genotypes_to_calc):
                if self.genotype_present.data.as_uchars[gprime] == 1: # otherwise, pops_and_params will be zero
                    for mprime in range(self.num_species):
                        if self.species_present.data.as_uchars[mprime] == 1: # otherwise, pops_and_params will be zero
                            delayed_index = mprime + gprime * self.num_species + delayed_base # + sex * num_species * num_genotypes
                            # unoptimised:
                            # denom += self.activity[mprime + mF * self.num_species] * pops_and_params[delayed_index]
                            # better:
                            ind = mprime + mF * self.num_species
                            denom = denom + self.activity.data.as_floats[ind] * pops_and_params[delayed_index]
            if denom > self.small_value:
                one_over_denom = 1.0 / denom
                some_males = 1
            else:
                one_over_denom = 0.0 # this ensures self.xprimeM = 0 below
            for gM in range(self.num_genotypes_to_calc):
                if self.genotype_present.data.as_uchars[gM] == 1: # if not present then x'M[ind] = 0
                    for mM in range(self.num_species):
                        if self.species_present.data.as_uchars[mM] == 1: # otherwise, pops_and_params will be zero
                            delayed_index = mM + gM * self.num_species + delayed_base # + sex * num_species * num_genotypes
                            ind = mF + mM * self.num_species + gM * self.num_species2
                            # unoptimised:
                            #self.xprimeM[ind] = self.activity[mM + mF * self.num_species] * pops_and_params[delayed_index] / denom
                            # better:
                            inda = mM + mF * self.num_species
                            self.xprimeM.data.as_floats[ind] = self.activity.data.as_floats[inda] * pops_and_params[delayed_index] * one_over_denom
        return some_males

    cdef void calcYYprime(self, unsigned some_males, unsigned delayed_base, float[:] pops_and_params):
        cdef unsigned s, g, m, yy_ind, mF, gF, f_ind, mM, gM, xprime_ind, precalc_ind
        cdef float xF
        array.zero(self.yy)
        array.zero(self.yyp)
        if some_males > 0:
            for s in range(self.num_sexes_to_calc):
                for g in range(self.num_genotypes_to_calc):
                    for m in range(self.num_species):
                    #if self.species_present.data.as_uchars[m] == 1: # can never generate a species if it does not exist
                        yy_ind = m + g * self.num_species + s * self.num_species * self.num_genotypes
                        for mF in range(self.num_species):
                            if self.species_present.data.as_uchars[mF] == 1: # can never generate a species if it does not exist
                                for gF in range(self.num_genotypes_to_calc):
                                    if self.genotype_present.data.as_uchars[gF] == 1: # if not present then x[ind] = 0
                                        f_ind = mF + gF * self.num_species + 1 * self.num_species * self.num_genotypes + delayed_base
                                        xF = pops_and_params[f_ind]
                                        if xF > self.small_value:
                                            for mM in range(self.num_species):
                                                if self.species_present.data.as_uchars[mM] == 1: # can never generate a species if it does not exist
                                                    for gM in range(self.num_genotypes_to_calc):
                                                        if self.genotype_present.data.as_uchars[gM] == 1: # if not present then x[ind] = 0
                                                            xprime_ind = mF + mM * self.num_species + gM * self.num_species2
                                                            precalc_ind = s + self.num_sexes * (g + self.num_genotypes * (gF + self.num_genotypes * (gM + self.num_genotypes * (m + self.num_species * (mF + self.num_species * mM)))))
                                                            # unoptimised:
                                                            #self.yy[yy_ind] += self.precalc[precalc_ind] * xF * self.xprimeM[xprime_ind]
                                                            #self.yyp[yy_ind] += self.precalcp[precalc_ind] * xF * self.xprimeM[xprime_ind]
                                                            # better:
                                                            self.yy.data.as_floats[yy_ind] = self.yy.data.as_floats[yy_ind] + self.precalc.data.as_floats[precalc_ind] * xF * self.xprimeM.data.as_floats[xprime_ind]
                                                            self.yyp.data.as_floats[yy_ind] = self.yyp.data.as_floats[yy_ind] + self.precalcp.data.as_floats[precalc_ind] * xF * self.xprimeM.data.as_floats[xprime_ind]

    cdef void calcCompetition(self, unsigned some_males):
        cdef unsigned m, mprime, alpha_ind, s, g, yy_ind
        cdef float alpha
        array.zero(self.comp)
        if some_males > 0:
            for m in range(self.num_species):
                for mprime in range(self.num_species):
                    alpha_ind = mprime + m * self.num_species
                    # unoptimised:
                    #alpha = self.competition[alpha_ind]
                    # better:
                    alpha = self.competition.data.as_floats[alpha_ind]
                    for s in range(self.num_sexes_to_calc):
                        for g in range(self.num_genotypes_to_calc):
                            if self.genotype_present.data.as_uchars[g] == 1:
                                yy_ind = m + g * self.num_species + s * self.num_species * self.num_genotypes
                                # unoptimised:
                                #self.comp[m] += alpha * self.yy[yy_ind]
                                # better:
                                self.comp.data.as_floats[m] = self.comp.data.as_floats[m] + alpha * self.yy.data.as_floats[yy_ind]
                if self.num_sexes_to_calc == 1: # assume male=female populations
                    self.comp.data.as_floats[m] = 2.0 * self.comp.data.as_floats[m]
