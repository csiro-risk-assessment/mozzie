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
    """Solves Mosquito ODE with 2 sexes and 3 genotypes.
    The number of populations is
    num_ages * num_species * num_genotypes * num_sexes = num_ages * num_species * 6.
    For species M, genotype G, sex S and age A, the index into pops_and_params is
    index = M + G * num_species + S * num_species * num_genotypes + A * num_species * num_genotypes * num_sexes
    Ages 0, 1, ..., num_ages - 2 are all larval stages.  Age = num_ages - 1 is the adult stage.
    Sex 0 is male.  Sex 1 is female.
    Genotype 0 is ww, Genotype 1 is Gw, Genotype 2 is GG.
    There is one spatially-varying parameter, that is the carrying capacity
    """

    # male, female, always.  So num_sexes = 2
    cdef unsigned num_sexes

    # ww, Gw, GG, always.  So num_genotypes = 3
    cdef unsigned num_genotypes

    # square of num_genotypes, viz 9
    cdef unsigned num_genotypes2

    # death rate of larvae
    cdef float mu_larvae

    # death rate of adults
    cdef float mu_adult

    # doco
    cdef float fecundity

    # rate of transferral from one age bracket to the next-eldest age bracket
    cdef float aging_rate

    # 1/(carrying capacity).  This is spatially-varying and is set upon entry to evolve()
    cdef float one_over_kk

    # inheritance_cube[i, j, k] = probability of mother genotype i, father genotype j producing offspring genotype k
    # where index 0, 1, 2 = ww, Gw, GG respectively
    cdef array.array inheritance_cube

    # inter-specific competition "matrix"
    cdef array.array alpha

    # hybridisation rate "matrix".  Default value is 1 if species_father==species_mother==species_offspring, and 0 otherwise"""
    cdef array.array hyb

    # relative-mating "matrix".  Default value is 1 if species_father==species_mother, and 0 otherwise
    cdef array.array mating

    # age categories are: larvae0, larvae1, larvae2, ..., larvaeN, adults
    cdef unsigned num_ages

    # number of species
    cdef unsigned num_species

    # square of num_species
    cdef unsigned num_species2

    # accuracy of PMB
    cdef float accuracy

    # competition vector, defined here for efficiency (so don't have to keep allocating memory)
    cdef array.array comp

    # denominator vector, defined here for efficiency (so don't have to keep allocating memory)
    cdef array.array denom

    # "matrix" array, defined here for efficiency (so don't have to keep allocating memory)
    cdef array.array mat

    # X array, defined here for efficiency (so don't have to keep allocating memory).  This holds the value of the populations in the "fun_for_scipy" method.  This is only used if solve_ivp is the time_integration_method.  If we remove the latter then we can remove these array.
    cdef array.array Xarray
    cdef float[:] cXarray

    # Runge-Kutta4 arrays that hold the values of populations (not parameters).  Sized to num_populations.
    cdef array.array rk1, rk2, rk3, rk4

    # Runge-Kutta "y" array that holds the value of populations and the carrying capacity: the C version is passed to computeRHS.  Sized to num_populations + num_parameters
    cdef array.array rky
    cdef float[:] crky

    # the RHS arrays in dX/dt = rhs.  Sized to num_populations.
    cdef array.array rhs

    # the change arrays.  For dX/dt = rhs using explicit-Euler timestepping, this is just dt * rhs.  Sized to num_populations.
    cdef array.array change
    cdef float[:] cchange

    # time-integration method (explicit_euler=0, solve_ivp=1, runge_kutta4=2,...)
    cdef unsigned time_integration_method

    # the minimum timestep allowed before the code crashes
    cdef float min_dt

    # if adaptive == 0 then do not do adaptive timestepping
    cdef unsigned adaptive

    # if a population < zero_cutoff at the end of a time-step then it is set to zero
    cdef float zero_cutoff

    # if the carrying capacity < min_cc, then no new larvae are produced
    cdef float min_cc
    # reciprocal of min_cc, for efficiency
    cdef float one_over_min_cc

    cdef inline unsigned getIndex(self, unsigned species, unsigned genotype, unsigned sex, unsigned age):
        """gets the index in pops_and_params corresponding to the given age, sex, genotype and species"""
        return species + self.num_species * (genotype + self.num_genotypes * (sex + self.num_sexes * age))

    cdef void setInheritance(self)
    """Sets the inheritance cube"""

    cdef inline float getInheritance(self, unsigned gt_father, unsigned gt_mother, unsigned gt_offspring):
        """Gets inhericance_cube[gt_father][gt_mother][gt_offpring]"""
        return self.inheritance_cube.data.as_floats[gt_father + gt_mother * self.num_genotypes + gt_offspring * self.num_genotypes2]

    cdef void setInternalParameters(self, unsigned num_ages, unsigned num_species, float accuracy)
    """Given num_ages, num_species and accuracy, set all internal parameters (num_populations, diffusing_indices, fecundity_proportion, etc) that depend on these"""

    cpdef setMinimumDt(self, float value)
    """Sets the minimum allowed timestep to value"""

    cpdef float getMinimumDt(self)
    """Gets the minimum allowed timestep"""

    cpdef setAdaptive(self, unsigned value)
    """Sets adaptive to value.  If value == 0 then do no adaptive timestepping"""

    cpdef unsigned getAdaptive(self)
    """Gets the minimum allowed timestep"""

    cpdef setZeroCutoff(self, float value)
    """Sets zero_cutoff to value.  If any population is less than this value at the end of the time-step, it is set to zero"""

    cpdef float getZeroCutoff(self)
    """Gets the value of zero_cutoff"""

    cpdef setMinCarryingCapacity(self, float value)
    """Sets min_cc to value.  If the carrying capacity is less than this value, no new larvae are produced"""

    cpdef float getMinCarryingCapacity(self)
    """Gets the value of min_cc"""

    cpdef setAlphaComponent(self, unsigned sp0, unsigned sp1, float value)
    """Sets alpha[sp0][sp1] = value (for inter-specific competition).  Note, if you setNumSpecies, alpha will be re-initialised to the identity"""

    cdef inline float getAlphaComponent(self, unsigned sp0, unsigned sp1):
        """Gets alpha[sp0][sp1]: the inter-specific competition"""
        return self.alpha.data.as_floats[sp0 + sp1 * self.num_species]

    cpdef setMatingComponent(self, unsigned species_father, unsigned species_mother, float value)
    """Sets mating[species_father][species_father] = value (for relative probability of inter-specific mating).  Note, if you setNumSpecies, mating will be re-initialised to the identity"""

    cdef inline float getMatingComponent(self, unsigned species_father, unsigned species_mother):
        """Gets mating[species_father][species_mother]: the inter-specific mating relative probability"""
        return self.mating.data.as_floats[species_father + species_mother * self.num_species]

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
    """Set number of species.  This calls setInternalParameters to appropriately set other parameters, given num_species.  It also reinitialises the alpha, hybridisation and mating matrices"""

    cpdef unsigned getNumSpecies(self)
    """Get number of species"""
    
    cpdef void setAccuracy(self, float accuracy)
    """Set accuracy of PMB.  This calls setInternalParameters to appropriately set other parameters, given accuracy"""

    cpdef float getAccuracy(self)
    """Get accuracy"""

    cdef float fecundity_proportion(self, unsigned offspring_sex, unsigned mother_gt, unsigned father_gt)
    """Returns the proportion of total fecundity for: mother genotype and father genotype to produce offspring sex"""

    cdef inline float getHybridisationRate(self, unsigned species_father, unsigned species_mother, unsigned species_offspring):
        """Returns the hybridisation rate for given father, mother and offspring"""
        return self.hyb.data.as_floats[species_father + species_mother * self.num_species + species_offspring * self.num_species2]

    cpdef setHybridisationRate(self, unsigned species_father, unsigned species_mother, unsigned species_offspring, float value)
    """Sets the hybridisation rate for the given father, mother and offspring.  Note, if you setNumSpecies, this will be reinitialised to its default value of 1 if species_father=species_mother=species_offspring, and 0 otherwise"""

    cdef void computeRHS(self, float[:] x)
    """Compute the rhs in dX/dt = rhs.  Here rhs is a funtion of x.  The result is put into self.rhs"""

    cdef void popChange(self, float timestep, float[:] current_pops_and_params, float[:] cchange)
    """Given the timestep and current_pops_and_params, compute the change in populations according to the ODE.
    For explicit-Euler timestepping, this is just timestep * rhs.
    For Runge-Kutta4 timestepping, this is given by the usual formula (1/6)(k_1 + 2k_2 + 2k_3 + k_4)"""

    cpdef setTimeIntegrationMethod(self, str method)
    """Sets the time integration method to be explicit_euler, solve_ivp, runge_kutta4"""


cdef class CellDynamicsMosquito23F(CellDynamicsMosquito23):
    """CellDynamicsMosquito23 but with modified fecundity"""

    cdef float fecundity_proportion(self, unsigned offspring_sex, unsigned mother_gt, unsigned father_gt)
    """Returns the proportion of total fecundity for: mother genotype and father genotype to produce offspring sex"""

cdef class CellDynamicsMosquito23G(CellDynamicsMosquito23F):
    """CellDynamicsMosquito23F but stochastic"""

    cdef void computeRHS_stoc(self, float[:] x)
    """Evaluates change in populations for a timestep (assuming dt = 1 day)"""
    
    cdef void popChange(self, float timestep, float[:] current_pops_and_params, float[:] cchange)
    """Given the timestep and current_pops_and_params, compute the change in populations according to the ODE.
    For this class we add stochastic timestepping, this is just rhs.

