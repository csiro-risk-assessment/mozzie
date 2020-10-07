import sys
import array
cimport cpython.array as array
import numpy as np
from scipy.integrate import solve_ivp
from math import exp, ceil, log, cos, sqrt
from libc.stdlib cimport rand, RAND_MAX

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
        self.small = 0.1

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
        cpdef float f = - self.mux + (1.0 - (x + self.axy * y) / (kx + self.small)) * (x / (x + self.w * y + self.small)) * self.gax
        cpdef float g = - self.muy + (1.0 - (self.ayx * x + y) / (ky + self.small)) * (self.gay + self.w * x / (x + self.w * y + self.small) * self.gax)
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
        self.num_ages = 2
        self.num_species = 1
        self.accuracy = 0.95

        self.num_species2 = 1

        # Set default carrying capacity
        self.one_over_kk = 1.0

        # size inheritance correctly
        self.inheritance_cube = array.clone(array.array('f', []), self.num_genotypes * self.num_genotypes * self.num_genotypes, zero = False)
        self.setInheritance()

        # allocate alpha array correctly, and set to the identity
        self.alpha = array.clone(array.array('f', []), 1, zero = True)
        self.setAlphaComponent(0, 0, 1.0)

        # allocate hyb array correctly, and set it to the identity
        self.hyb = array.clone(array.array('f', []), 1, zero = True)
        self.setHybridisationRate(0, 0, 0, 1.0)

        # allocate the mating array correctly, and set it to the identiy
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
        self.accuracy = accuracy
        
        self.num_populations = self.num_ages * self.num_sexes * self.num_genotypes * self.num_species

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

    cpdef setMinimumDt(self, float value):
        self.min_dt = value

    cpdef float getMinimumDt(self):
        return self.min_dt

    cpdef setAdaptive(self, unsigned value):
        self.adaptive = value

    cpdef unsigned getAdaptive(self):
        return self.adaptive

    cpdef setZeroCutoff(self, float value):
        self.zero_cutoff = value

    cpdef float getZeroCutoff(self):
        return self.zero_cutoff

    cpdef setMinCarryingCapacity(self, float value):
        self.min_cc = value
        self.one_over_min_cc = 1.0 / self.min_cc

    cpdef float getMinCarryingCapacity(self):
        return self.min_cc

    cpdef setAlphaComponent(self, unsigned sp0, unsigned sp1, float value):
        if sp0 >= self.num_species or sp1 >= self.num_species:
            raise ValueError("sp0 " + str(sp0) + " and sp1 " + str(sp1) + " must be less than the number of species, " + str(self.num_species))
        self.alpha.data.as_floats[sp0 + sp1 * self.num_species] = value

    cpdef setHybridisationRate(self, unsigned species_father, unsigned species_mother, unsigned species_offspring, float value):
        if species_father >= self.num_species or species_mother >= self.num_species or species_offspring >= self.num_species:
            raise ValueError("All species numbers, " + str(species_father) + ", " + str(species_mother) + ", " + str(species_offspring) + " must be less than the number of species, " + str(self.num_species))
        self.hyb.data.as_floats[species_father + species_mother * self.num_species + species_offspring * self.num_species2] = value

    def getHybridisationRateFromPython(self, unsigned species_father, unsigned species_mother, unsigned species_offspring):
        """Returns the hybridisation rate for given father, mother and offspring.  This is a slow interface: use getAlphaComponent from all cython code"""
        if species_father >= self.num_species or species_mother >= self.num_species or species_offspring >= self.num_species:
            raise ValueError("All species numbers, " + str(species_father) + ", " + str(species_mother) + ", " + str(species_offspring) + " must be less than the number of species, " + str(self.num_species))
        return self.getHybridisationRate(species_father, species_mother, species_offspring)

    cpdef setMatingComponent(self, unsigned species_father, unsigned species_mother, float value):
        if species_father >= self.num_species or species_mother >= self.num_species:
            raise ValueError("species_father " + str(species_father) + " and species_mother " + str(species_mother) + " must be less than the number of species, " + str(self.num_species))
        self.mating.data.as_floats[species_father + species_mother * self.num_species] = value


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

    cpdef unsigned getNumSpecies(self):
        return self.num_species
    
    cpdef void setAccuracy(self, float accuracy):
        self.setInternalParameters(self.num_ages, self.num_species, accuracy)

    cpdef float getAccuracy(self):
        return self.accuracy

    cdef void computeRHS(self, float[:] x):
        """Evaluates d(populations)/dt"""

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

        # newborn larvae
        if self.one_over_kk < self.one_over_min_cc:

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
                                        self.mat.data.as_floats[ind_mat] = self.mat.data.as_floats[ind_mat] + self.getHybridisationRate(spm, spf, sp) * self.getInheritance(gtm, gtf, gt) * self.getMatingComponent(spm, spf) * x[ind] * self.fecundity_proportion(sex, gtf, gtm)
                                # multiply mat by things that don't depend on gtm or spm
                                self.mat.data.as_floats[ind_mat] = self.mat.data.as_floats[ind_mat] * self.comp.data.as_floats[sp] * self.fecundity * self.denom.data.as_floats[spf]
            

        # mortality, and aging into/from neighbouring age brackets
        for sex in range(self.num_sexes):
            for gt in range(self.num_genotypes):
                for sp in range(self.num_species):
                    age = 0
                    row = self.getIndex(sp, gt, sex, age)
                    ind = row + self.num_populations * row # diagonal entry
                    if self.num_ages > 1:
                        self.mat.data.as_floats[ind] = self.mat.data.as_floats[ind] - self.mu_larvae # mortality
                        self.mat.data.as_floats[ind] = self.mat.data.as_floats[ind] - self.aging_rate # aging to older bracket
                    for age in range(1, self.num_ages - 1):
                        row = self.getIndex(sp, gt, sex, age)
                        col = self.getIndex(sp, gt, sex, age)
                        ind = col + self.num_populations * row # diagonal entry
                        self.mat.data.as_floats[ind] = self.mat.data.as_floats[ind] - self.mu_larvae - self.aging_rate # mortality and aging to older bracket
                        col = self.getIndex(sp, gt, sex, age - 1)
                        ind = col + self.num_populations * row # below diagonal
                        self.mat.data.as_floats[ind] = self.mat.data.as_floats[ind] + self.aging_rate # contribution from younger age bracket
                    age = self.num_ages - 1
                    row = self.getIndex(sp, gt, sex, age)
                    ind = row + self.num_populations * row # diagonal entry
                    self.mat.data.as_floats[ind] = self.mat.data.as_floats[ind] - self.mu_adult # mortality
                    if self.num_ages > 1:
                        col = self.getIndex(sp, gt, sex, age - 1)
                        ind = col + self.num_populations * row # below diagonal
                        self.mat.data.as_floats[ind] = self.mat.data.as_floats[ind] + self.aging_rate # contribution from younger age bracket


        #for row in range(self.num_populations):
        #    for col in range(self.num_populations):
        #        sys.stdout.write(str(self.mat[col + self.num_populations * row]) + " ")
        #    sys.stdout.write("\n")

        array.zero(self.rhs)
        for row in range(self.num_populations):
            for col in range(self.num_populations):
                self.rhs.data.as_floats[row] += self.mat.data.as_floats[col + self.num_populations * row] * x[col]

    cpdef void evolve(self, float timestep, float[:] pops_and_params):
        cdef unsigned ind

        if pops_and_params[self.num_populations] <= self.min_cc:
            self.one_over_kk = 2 * self.one_over_min_cc # signal to computeRHS that the carrying capacity has fallen below min_cc, so no newborns will be produced
        else:
            self.one_over_kk = 1.0 / pops_and_params[self.num_populations]        

        cdef float time_done = 0.0
        cdef float dt = timestep
        cdef float new_dt = timestep
        while time_done < timestep:
            # try timestep dt
            self.popChange(dt, pops_and_params, self.cchange)
            new_dt = timestep # biggest value that ever needs to be used
            for ind in range(self.num_populations):
                if pops_and_params[ind] + self.cchange[ind] < 0.0:
                    # the change will send the population negative, which is deemed to be erroneous
                    # and due to dt being too big.  Using explicit-Euler timestepping, if we choose
                    # new_dt = dt * pops_and_params[ind] / (-cchange[ind]), then running popChange again
                    # will give exactly cchange[ind] = -pops_and_params[ind].  So multiply this
                    # new_dt by 0.9 to give a factor of safety
                    new_dt = min(new_dt, - 0.9 * dt * pops_and_params[ind] / self.cchange[ind])
            if self.adaptive == 0 or new_dt == timestep:
                # not doing adaptivity, or none of the populations went negative, so success:
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
            self.crky[self.num_populations] = current_pops_and_params[self.num_populations] # the carrying capacity
            self.computeRHS(self.crky)
            for ind in range(self.num_populations):
                self.rk2.data.as_floats[ind] = timestep * self.rhs.data.as_floats[ind]
            # step 3
            for ind in range(self.num_populations):
                self.crky[ind] = current_pops_and_params[ind] + 0.5 * self.rk2.data.as_floats[ind]
            self.crky[self.num_populations] = current_pops_and_params[self.num_populations] # the carrying capacity
            self.computeRHS(self.crky)
            for ind in range(self.num_populations):
                self.rk3.data.as_floats[ind] = timestep * self.rhs.data.as_floats[ind]
            # step 4
            for ind in range(self.num_populations):
                self.crky[ind] = current_pops_and_params[ind] + self.rk3.data.as_floats[ind]
            self.crky[self.num_populations] = current_pops_and_params[self.num_populations] # the carrying capacity
            self.computeRHS(self.crky)
            for ind in range(self.num_populations):
                self.rk4.data.as_floats[ind] = timestep * self.rhs.data.as_floats[ind]
            # put it all together
            for ind in range(self.num_populations):
                cchange[ind] = (1.0 / 6.0) * (self.rk1.data.as_floats[ind] + 2 * self.rk2.data.as_floats[ind] + 2 * self.rk3.data.as_floats[ind] + self.rk4.data.as_floats[ind])


    cpdef setTimeIntegrationMethod(self, str method):
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
        
        # newborn larvae
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
        self.num_genotypes = 6 # ww, wc, wr, cc, cr, rr always
        
        self.k_c = 0.995
        self.k_j = 0.02
        self.k_ne = 0.0001
        self.w_prob = 0.5*(1. - self.k_c)
        self.c_prob = 0.5*(1 + self.k_c*(1. - self.k_j)*(1. - self.k_ne))
        self.r_prob = 1. - self.w_prob - self.c_prob
        
        # size inheritance correctly
        self.inheritance_cube = array.clone(array.array('f', []), self.num_genotypes * self.num_genotypes * self.num_genotypes, zero = False)
        self.setInheritance()
        
    cdef void setInheritance(self):
        inheritance_list = [[[1., 0., 0., 0., 0., 0.], # ww x ww
                             [self.w_prob, self.c_prob, self.r_prob, 0., 0., 0.], # ww x wc
                             [1., 0., 1., 0., 0., 0.], # ww x wr
                             [0., 1., 0., 0., 0., 0.], # ww x cc
                             [0., 0.5, 0.5, 0., 0., 0.], # ww x cr
                             [0., 0., 1., 0., 0., 0.]],# ww x rr
                            [[self.w_prob, self.c_prob, self.r_prob, 0., 0., 0.], # wc x ww
                             [self.w_prob*self.w_prob, 2.*self.w_prob*self.c_prob, 2.*self.w_prob*self.r_prob, self.c_prob*self.c_prob, 2.*self.c_prob*self.r_prob., self.r_prob*self.r_prob], # wc x wc
                             [0.5*self.w_prob, 0.5*self.c_prob, 0.5*(self.w_prob + self.r_prob), 0., 0.5*self.c_prob, 0.5*self.r_prob], # wc x wr
                             [0., self.w_prob, 0., self.c_prob, self.r_prob, 0.], # wc x cc
                             [0., 0.5*self.w_prob, 0.5*self.w_prob, 0.5*self.c_prob, 0.5*(self.c_prob + self.r_prob), 0.5*self.r_prob], # wc x cr
                             [0., 0., self.w_prob, 0., self.c_prob, self.r_prob]],# wc x rr
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
                             [0., 0., 0., 0., 0., 1.]]]# rr x rr
        cdef unsigned gt_father, gt_mother, gt_offspring
        for gt_father in range(self.num_genotypes):
            for gt_mother in range(self.num_genotypes):
                for gt_offspring in range(self.num_genotypes):
                    self.inheritance_cube.data.as_floats[gt_father + gt_mother * self.num_genotypes + gt_offspring * self.num_genotypes2] = inheritance_list[gt_father][gt_mother][gt_offspring]
