# Copyright (c) 2024 Commonwealth Scientific and Industrial Research Organisation (CSIRO) ABN 41 687 119 230.
import array
cimport cpython.array as array

cdef class CellDynamicsBase:
    """Manipulates information at a single cell, in particular this class solves lifecycle ODEs"""

    # number of populations (maleGw, femaleGG, larvaeMaleGw, whatever) in the Cell dynamics
    cdef unsigned num_populations

    # number of populations that diffuse
    cdef unsigned num_diffusing

    # diffusing_indices[i] is the i^th population that is diffusing
    # 0 <= diffusing_indices[i] < num_populations.  0 <= i < num_diffusing
    cdef array.array diffusing_indices

    # number of populations that advect
    cdef unsigned num_advecting

    # advecting_indices[i] is the i^th population that is advecting
    # 0 <= advecting_indices[i] < num_populations.  0 <= i < num_advecting
    cdef array.array advecting_indices

    # advection_class[i] is the class of the i^th population that is advecting
    # 0 <= advection_class[i], with 0 <= i < num_advecting
    # This allows for different types of advection, for instance, males (class = 0) advecting with a different probability than females (class = 1)
    cdef array.array advection_class

    # number of parameters (carrying capacity, mortality rate, etc) in the Cell dynamics
    cdef unsigned num_parameters

    # A value that can be used as a lower bound on population numbers, carrying capacities, or anything else, to eliminate underflows or to zero populations when they reach a threshold
    cdef float small_value

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

    cpdef array.array getAdvectionClass(self)
    """Returns: advection_class[i], which is the class of the i^th population that is advecting
    0 <= advection_class[i], with 0 <= i < num_advecting.  This allows for different types of
    advection, for instance, males (class = 0) advecting with a different probability than
    females (class = 1)"""

    cpdef setAdvectionClass(self, unsigned population_index, unsigned new_class)
    """Sets the advection class of the given population_index to new_class.  0 <= population_index < num_populations.
    This will raise a ValueError if the population_index has not been defined to be advecting"""



    cpdef unsigned getNumberOfParameters(self)
    """Returns the number of parameters (carrying capacity, mortality rate, etc) in the cell dynamics"""

    cpdef void setSmallValue(self, float small_value)
    """Sets the small_value"""

    cpdef float getSmallValue(self)
    """Gets self.small_value"""

    cpdef evolve(self, float timestep, float[:] pops_and_params)
    """Performs one timestep of evolution
    Note: pops_and_params is a pointer to an array of floats, so any changes made to its values will be evident to the calling function
    Note: pops_and_params will be of size num_populations + num_parameters, and should not be resized
    Note: the first num_populations entries of pops_and_params will be the population values, while the remainder are the parameter values"""

    cpdef array.array calcQm(self, float[:] eqm_pops_and_params)
    """This currently only works for CellDynamicsMosquitoBH26Delay.  Given the eqm_pops_and_params, which is an array containing the populations at equilibrium, return the qm values.  Only the 'current_index' and wild-type entries of eqm_pops_and_params are used in this calculation.  This function assumes the equilibrium populations for male equal those for females.  You must make sure the entries corresponding to adults of species M, genotype=0, sex S are correct (they have index = M + 0 * num_species + S * num_species * num_genotypes + current_index * num_species * num_genotypes * num_sexes).  This function does not work if the equilibrium populations contain genotype!=0 mosquitoes, and only works if male=female=carrying_capacity/2.  This function does not check that eqm_pops_and_params is actually an equilibrium: if you feed it garbage, it will produce garbage!"""

    cpdef unsigned getNumSpecies(self)
    """Get number of species"""

cdef class CellDynamicsStatic15_9_3_2(CellDynamicsBase):
    """No dynamics within the cell (all populations are static as far as the cell is concerned)
    15 populations, 9 of them are diffusing and 3 advecting (all with advection_class=0), with 2 parameters"""
    

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
    cdef float mux # mu_x
    cdef float muy # mu_y
    cdef float gax # gamma_x
    cdef float gay # gamma_y
    cdef float axy # alpha_{xy}
    cdef float ayx # alpha_{yx}
    cdef float w   # w
    cdef float small

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

    cpdef void setSmall(self, float small)
    """Sets small"""

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

    cpdef float getSmall(self)
    """Gets small"""

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

    # 1/(carrying capacity).  This is of size num_parameters.  This is spatially-varying and is set upon entry to evolve()
    cdef array.array one_over_kk

    # inheritance_cube[i, j, k] = probability of mother genotype i, father genotype j producing offspring genotype k
    # where index 0, 1, 2 = ww, Gw, GG respectively
    cdef array.array inheritance_cube

    # inter-specific competition "matrix"
    cdef array.array alpha

    # hybridisation rate "matrix".  Default value is 1 if species_father==species_mother==species_offspring, and 0 otherwise"""
    cdef array.array hyb

    # relative-mating "matrix".  Default value is 1 if species_father==species_mother, and 0 otherwise
    cdef array.array mating

    # relative fitness "vector".  Default value is 1 for all genotypes
    cdef array.array fitness

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
    
    # Dummy floats using evolve
    cdef float species_stuff
    cdef float genotype_stuff
    cdef float tmp_float
    cdef float xcol

    # Array to hold intermediate values of sums, sized to num_sexes * num_genotypes**2
    ###################################################################################
    #
    # NOTE: this array makes the code much more fragile to reorderings of function calls
    #
    # genotypeRapidAccess is sized whenever the following functions are called:
    # setNumGenotypes
    # during construction (via setNumGenotypes)
    # Hence, self.num_sexes must be set appropriately before setNumGenotypes in the constructor
    #
    # genotypeRapidAccess is computed during construction and whenever the following functions are called:
    # setInternalParameters
    # setFitnessComponent
    # setNumGenotypes (via setFitnessComponent)
    # during construction (via setNumGenotypes and setInternalParameters)
    # setNumAges (which calls setInternalParameters)
    # setNumSpecies (which calls setInternalParameters)
    # setAccuracy (which calls setInternalParameters)
    # Hence, before any of these things, self.num_sexes must be set appropriately, and self.accuracy mustn't be zero (setAccuracy and setInternalParameters set self.accuracy before calling setGenotypeRapidAccess
    #
    ###################################################################################
    cdef array.array genotypeRapidAccess

    # Array to indicate presence or absence of species, sized to num_species
    cdef array.array species_present

    # Array to indicate presence or absence of genotype, sized to num_genotypes
    cdef array.array genotype_present

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

    cpdef setFitnessComponent(self, unsigned genotype, float value)
    """Sets fitness[genotype] = value (for relative reproductive fitness of males of each genotype).  Note, if you setNumGenotypes, fitness will be re-initialised to the identity"""

    cdef inline float getMatingComponent(self, unsigned species_father, unsigned species_mother):
        """Gets mating[species_father][species_mother]: the inter-specific mating relative probability"""
        return self.mating.data.as_floats[species_father + species_mother * self.num_species]

    cdef inline float getFitnessComponent(self, unsigned genotype):
        """Gets fitness[genotype]: the relative reproductive fitness of males of each genotype"""
        return self.fitness.data.as_floats[genotype]

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

    cpdef setNumGenotypes(self, unsigned num_sexes, unsigned num_genotypes)
    """Set number of genotypes. This calls setInheritance and setFitnessComponent to appropriately initialise the inheritance tensor and fitness vector.  This depends on num_sexes because it also sizes and calculates genotypeRapidAccess, and you should ensure that self.num_sexes is set appropriately before calling setNumGenotypes"""

    cpdef void setGenotypeRapidAccess(self)
    """Sets the values in self.genotypeRapidAccess.  Before calling this, self.num_sexes should have been set, since self.genotypeRapidAccess depends on that parameter.  Also, self.num_genotypes should have been set using setNumGenotypes(...).  Also, self.accuracy should not be zero.  Please see extended comments associated with the genotypeRapidAccess array"""

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
    """Compute the rhs in dX/dt = rhs.  Here rhs is a function of x.  The result is put into self.rhs"""

    cdef void popChange(self, float timestep, float[:] current_pops_and_params, float[:] cchange)
    """Given the timestep and current_pops_and_params, compute the change in populations according to the ODE.
    For explicit-Euler timestepping, this is just timestep * rhs.
    For Runge-Kutta4 timestepping, this is given by the usual formula (1/6)(k_1 + 2k_2 + 2k_3 + k_4)"""

    cpdef setTimeIntegrationMethod(self, str method)
    """Sets the time integration method to be explicit_euler, solve_ivp, runge_kutta4"""


cdef class CellDynamicsMosquito23F(CellDynamicsMosquito23):
    """CellDynamicsMosquito23 but with modified fecundity"""

    cdef float fecundity_proportion(self, unsigned offspring_sex, unsigned mother_gt, unsigned father_gt)
    """Returns the proportion of total fecundity for: mother genotype and father genotype to produce offspring sex.  This is overwritten by a special from to the 23F class"""

    cpdef float getFecundityProportion(self, unsigned offspring_sex, unsigned mother_gt, unsigned father_gt)
    """Returns the fecundity proportion for this 23F class"""

cdef class CellDynamicsMosquito26(CellDynamicsMosquito23):
    """Solves Mosquito ODE with
      num_sexes = 2 (cannot be changed)
      num_genotypes = 6 (should not be changed)
      num_species = 2 (by default.  Should not be changed unless you reinterpret the carrying-capacities)
    The number of populations is
    num_ages * num_species * num_genotypes * num_sexes = num_ages * 24
    For species M, genotype G, sex S and age A, the index into pops_and_params is
    index = M + G * num_species + S * num_species * num_genotypes + A * num_species * num_genotypes * num_sexes
    Ages 0, 1, ..., num_ages - 2 are all larval stages.  Age = num_ages - 1 is the adult stage.
    Sex 0 is male.  Sex 1 is female.
    Genotypes:
      0 is ww
      1 is wc
      2 is wr
      3 is cc
      4 is cr
      5 is rr
    There are two spatially-varying parameters, which are the carrying capacities for each species
    See the constructor for more details such as the default values for:
      inheritance
      alpha
      accuracy
      fitness
      hybridisation
      mating components
    """

    cdef float w_prob
    cdef float c_prob
    cdef float r_prob

    cdef void setInheritance(self)
    """Version of setInheritance for 6 genotypes"""

    cpdef setFitnessComponents26(self, float h_e, float h_n, float s_e, float s_n)
    """sets FitnessParams assuming there are 6 genotypes: fitness(0) = fitness(2) = fitness(5) = 1; fitness(1) = fitness(3) = (1 - h_e * s_e) * (1 - h_n * s_n); fitness(4) = (1 - s_e) * (1 - s_n)"""

    cpdef setInheritance26(self, float k_c, float k_j, float k_ne)
    """sets the inheritance cube based on k_c, k_j and k_ne, assuming there are 6 genotypes"""

cdef class CellDynamicsDelayBase(CellDynamicsBase):
    """Base class for lifecycle dynamics as governed by a delay differential equation

    An example of a delay equation is
    dx/dt = x(t - delay * dt)
    Here "delay" is the number of dt involved in the delay equation
    So, if delay = 0 then the equations just revert to ODEs.

    An important variable is current_index, which is an integer between 0 and delay.
    This defines where the current populations are located in pops_and_params.
    It is initialised to 0 in the constructor.
    Examples that will appear in derived classes:
      - Adults of species M, genotype G, sex S have index = M + G * num_species + S * num_species * num_genotypes + current_index * num_species * num_genotypes * num_sexes
      - The adult population at the previous dt has index = M + G * num_species + S * num_species * num_genotypes + (current_index - 1)%(delay + 1) * num_species * num_genotypes * num_sexes
      - The adult population at N dt-steps ago has index = M + G * num_species + S * num_species * num_genotypes + (current_index - N)%(delay + 1) * num_species * num_genotypes * num_sexes
      - The population at delay dt-steps ago has index = M + G * num_species + S * num_species * num_genotypes + (current_index - delay)%(delay + 1) * num_species * num_genotypes * num_sexes = = M + G * num_species + S * num_species * num_genotypes + (current_index + 1)%(delay + 1) * num_species * num_genotypes * num_sexes

    Hence, for the simple process equation of the form dx/dt = x(t - delay * dt), evolve could simply set:
    delayed_index = (current_index + 1)%(delay + 1)
    new_pop = pops_and_params[current_index] + dt * pops_and_params[delayed_index]
    pops_and_params[delayed_index] = new_pop

    !!! IMPORTANT NOTE !!!
    After evolve is called for all grid cells, current_index must be incremented!!  This is the responsibility of the calling code, for the cellDynamics class has no way of knowing evolve has been called for all grid cells.  (An alternative would be to make current_index one of the spatially-varying parameters, and evolve could update it at every call of evolve, but this would consume memory.)
    """
    # number of dt in the delay equation.  (Set to 1 in constructor)
    cdef unsigned delay

    # index of the current adult population.  (Set to zero in constructor)
    cdef unsigned current_index

    cpdef unsigned getDelay(self)
    """Returns delay"""

    cpdef void incrementCurrentIndex(self)
    """Sets current_index = (current_index + 1) % (delay + 1).  In derived classes, this will also update advecting_indices and diffusing_indices"""

    cpdef unsigned getCurrentIndex(self)
    """returns current_index"""

cdef class CellDynamics26DelayBase(CellDynamicsDelayBase):
    """Base class for any ODE with 2 sexes, 6 genotypes, using a delay equation

    An example of a delay equation is
    dx/dt = x(t - delay * dt)
    Here "delay" is the number of dt involved in the delay equation
    So, if delay = 0 then the equations just revert to ODEs.

    The number of populations is
    (delay + 1) * num_species * num_geotypes * num_sexes = (delay + 1) * num_species * 12.

    An important variable is current_index, which is an integer between 0 and delay.
    This defines where the current populations are located in pops_and_params.
    It is initialised to 0 in the constructor.
    Examples:
      - Adults of species M, genotype G, sex S have index = M + G * num_species + S * num_species * num_genotypes + current_index * num_species * num_genotypes * num_sexes
      - The adult population at the previous dt has index = M + G * num_species + S * num_species * num_genotypes + (current_index - 1)%(delay + 1) * num_species * num_genotypes * num_sexes
      - The adult population at N dt-steps ago has index = M + G * num_species + S * num_species * num_genotypes + (current_index - N)%(delay + 1) * num_species * num_genotypes * num_sexes
      - The population at delay dt-steps ago has index = M + G * num_species + S * num_species * num_genotypes + (current_index - delay)%(delay + 1) * num_species * num_genotypes * num_sexes = = M + G * num_species + S * num_species * num_genotypes + (current_index + 1)%(delay + 1) * num_species * num_genotypes * num_sexes

    Hence, for the simple process equation of the form dx/dt = x(t - delay * dt), evolve could simply set:
    delayed_index = (current_index + 1)%(delay + 1)
    new_pop = pops_and_params[current_index] + dt * pops_and_params[delayed_index]
    pops_and_params[delayed_index] = new_pop

    !!! IMPORTANT NOTE !!!
    After evolve is called for all grid cells, current_index must be incremented!!  This is the responsibility of the calling code, for the cellDynamics class has no way of knowing evolve has been called for all grid cells.  (An alternative would be to make current_index one of the spatially-varying parameters, and evolve could update it at every call of evolve, but this would consume memory.)
    
    Sex 0 is male.  Sex 1 is female.
    Genotypes:
      0 is ww
      1 is wc
      2 is wr
      3 is cc
      4 is cr
      5 is rr
    """
    # male, female, always.  So num_sexes = 2
    cdef unsigned num_sexes

    # 0 is ww, 1 is wc, 2 is wr, 3 is cc, 4 is cr, 5 is rr.  So num_genotypes = 3
    cdef unsigned num_genotypes

    # square of num_genotypes, viz 9
    cdef unsigned num_genotypes2

    # inheritance_cube[gF, gM, g] = probability of mother genotype gF, father genotype gM producing offspring genotype g.  This is arranged as a vector with index = gM + gF * num_genotypes + g * num_genotypes * num_genotypes
    cdef array.array inheritance_cube

    # number of species (set to 3 in constructor)
    cdef unsigned num_species

    # square of num_species
    cdef unsigned num_species2

    # death rate of mosquito type M and genotype G and sex S  has index M + G * num_species + S * num_species * num_genotypes (all set to 1.0 in constructor)
    cdef array.array death_rate

    # competition (alpha) between mosquito type M and type M' has index M' + M * num_species
    cdef array.array competition

    # emergence rate[mF] (lambda) of type mF.  Larvae per wildtype female at the end of the larval duration period, per day, per adult female, in the absence of density dependence
    cdef array.array emergence_rate

    # activity level (a) of female type mF and male type mM has index mM + mF * num_species
    cdef array.array activity

    # parameters involved in the inheritance cube (set in the constructor)
    cdef float m_w
    cdef float m_c
    cdef float w_prob
    cdef float c_prob
    cdef float r_prob

    # probability that offspring of wc or cc fathers are male (paternal male bias has sex_ratio > 0.5)
    cdef float sex_ratio

    # probability that offspring of (mother wc or cc + father ww) is female (usually female_bias > 0.5)
    cdef float female_bias

    # fecundity[gM][gF][s] = proportion (male gM + female gF) producing offspring of sex s.  This is a vector with index = gM + gF * num_genotypes + s * num_genotypes * num_genotypes
    cdef array.array fecundity_p

    # reduction[gM][gF] = reduced number of adults because of construct.  This is a vector with index = gF + gM * num_genotypes
    cdef array.array reduction

    # hybridisation[mM][mF][m] = probability that offspring of species m results from male of species mM and female of species mF.  This is a vector with index m + mF * num_species + mM * num_species * num_species
    cdef array.array hybridisation

    # offspring_modifier[s][mM][mF] = suppression (if <1, or increased vigour if >1) of offspring of sex s that arises from male of species mM and female of species mF.  This is a vector with index mF + mM * num_species + s * num_species * num_species
    cdef array.array offspring_modifier

    # whether precalculate() has been called with the most up-to-date information concerning number of species, emergence rates, inheritance, fecundity, reduction, hybridisation, offspring_modifier (have_precalculated = 0 means precalculate() has not been called with the most up-to-date information)
    cdef int have_precalculated

    cpdef setDeathRate(self, list death_rate)
    """sets self.death_rate to death_rate.
    The death_rate list must be of the form death_rate[sex][genotype][mosquito_species]
    All elements must be positive.
    Note that evolve uses exp(-death_rate * dt).  You must ensure that death_rate and dt are set to this does not overflow."""

    cpdef list getDeathRate(self)
    """Returns death_rate[sex][genotype][mosquito_species]"""

    cpdef setCompetition(self, list competition)
    """sets self.competition to competition[species][species_prime]
    This is called "alpha" is the documentation"""

    cpdef list getCompetition(self)
    """Returns competition[species][species_prime].  This is called "alpha" is the documentation"""

    cpdef setEmergenceRate(self, list emergence_rate)
    """sets self.emergence_rate to emergence_rate.
    The emergence_rate list must be num_species in length, and must be a list of non-negative floats."""

    cpdef list getEmergenceRate(self)
    """Returns emergence_rate[species]"""

    cpdef setActivity(self, list activity)
    """sets self.activity to activity[speciesFemale][speciesMale].
    Each element must be a list of non-negative floats."""

    cpdef list getActivity(self)
    """Returns activity[speciesFemale][speciesMale]"""

    cpdef setReduction(self, list reduction)
    """sets self.reduction to reduction[gM][gF]"""

    cpdef list getReduction(self)
    """returns reduction[gM][gF]"""

    cpdef setHybridisation(self, list hybridisation)
    """sets self.hybridisation to hybridisation[mM][mF][m]"""

    cpdef list getHybridisation(self)
    """returns hybridisation[mM][mF][m]"""

    cpdef setOffspringModifier(self, list offspring_modifier)
    """sets self.offspring_modifier to offspring_modifier[s][mM][mF]"""

    cpdef list getOffspringModifier(self)
    """returns offspring_modifier[s][mM][mF]"""

    cdef void setInheritance(self)
    """Version of setInheritance for 6 genotypes"""

    cpdef array.array getInheritance(self)
    """returns inheritance_cube"""

    cpdef evolveTrial(self, float timestep, float[:] pops_and_params)
    """This just implements dx/dt = -death_rate * x + lambdah * x[t - delay * dt]: used for testing
    If desired, it can be removed in production code"""

    cpdef void setFecundityP(self, float sex_ratio, float female_bias)
    """Sets sex_ratio, female_bias and fecundity_p"""

    cpdef float getSexRatio(self)
    """Returns sex_ratio"""

    cpdef float getFemaleBias(self)
    """Returns female_bias"""        

    cpdef list getFecundityP(self)
    """Returns fecundity[gM][gF][s] = proportion (male gM + female gF) producing offspring of sex s.  This is a vector with index = gM + gF * num_genotypes + s * num_genotypes * num_genotypes
"""

    cpdef precalculate(self)
    """Sets have_precalculated = 1.  This method may be over-ridden by derived classes in order to calculate spatially-independent things after number of species, emergence rates, inheritance, fecundity, reduction, hybridsation or offspring_modifier have changed"""

cdef class CellDynamicsMosquitoLogistic26Delay(CellDynamics26DelayBase):
    """Solves Mosquito ODE with 2 sexes, 6 genotypes, using a logistic delay equation

    An example of a delay equation is
    dx/dt = x(t - delay * dt)
    Here "delay" is the number of dt involved in the delay equation
    So, if delay = 0 then the equations just revert to ODEs.

    The number of populations is
    (delay + 1) * num_species * num_geotypes * num_sexes = (delay + 1) * num_speces * 12.

    An important variable is current_index, which is an integer between 0 and delay.
    This defines where the current populations are located in pops_and_params.
    It is initialised to 0 in the constructor.
    Examples:
      - Adults of species M, genotype G, sex S have index = M + G * num_species + S * num_species * num_genotypes + current_index * num_species * num_genotypes * num_sexes
      - The adult population at the previous dt has index = M + G * num_species + S * num_species * num_genotypes + (current_index - 1)%(delay + 1) * num_species * num_genotypes * num_sexes
      - The adult population at N dt-steps ago has index = M + G * num_species + S * num_species * num_genotypes + (current_index - N)%(delay + 1) * num_species * num_genotypes * num_sexes
      - The population at delay dt-steps ago has index = M + G * num_species + S * num_species * num_genotypes + (current_index - delay)%(delay + 1) * num_species * num_genotypes * num_sexes = = M + G * num_species + S * num_species * num_genotypes + (current_index + 1)%(delay + 1) * num_species * num_genotypes * num_sexes

    Hence, for the simple process equation of the form dx/dt = x(t - delay * dt), evolve could simply set:
    delayed_index = (current_index + 1)%(delay + 1)
    new_pop = pops_and_params[current_index] + dt * pops_and_params[delayed_index]
    pops_and_params[delayed_index] = new_pop

    !!! IMPORTANT NOTE !!!
    After evolve is called for all grid cells, current_index must be incremented!!  This is the responsibility of the calling code, for the cellDynamics class has no way of knowing evolve has been called for all grid cells.  (An alternative would be to make current_index one of the spatially-varying parameters, and evolve could update it at every call of evolve, but this would consume memory.)
    
    Sex 0 is male.  Sex 1 is female.
    Genotypes:
      0 is ww
      1 is wc
      2 is wr
      3 is cc
      4 is cr
      5 is rr
    There are num_species spatially-varying parameters, which are the carrying-capacities of each species
    """

    # this is used in evolve to hold the new population, new_pop[s, g, m] with index m + g * self.num_species + s * self.num_species * self.num_genotypes
    cdef array.array new_pop
    
    # this is used in evolve to hold xprimeM (proportionate mixing quantity).  xprimeM[gM][mM][mF] has index mF + mM * num_species + gM * num_species * num_species
    cdef array.array xprimeM

    # Y used in evolve.  yy[sex][genotype][m] has index m + genotype * num_species + sex * num_species * num_genotypes
    cdef array.array yy

    # Y' used in evolve.  yyp[sex][genotype][m] has index m + genotype * num_species + sex * num_species * num_genotypes
    cdef array.array yyp

    # comp used in evolve.  comp[m]
    cdef array.array comp

    # precalc is used in evolve.  precalc[mM, mF, m, gM, gF, g, s] = hybridisation[mM, mF, m] * emergence_rate[mF] * inheritance_cube[gM, gF, g] * fecundity[gM, gF, s] * reduction[gM, gF].  It has index s + num_sexes * (g + num_genotypes * (gF + num_genotypes * (gM + num_genotypes * (m + num_species * (mF + num_species * mM))))), so is of size num_sexes * num_genotypes^3 * num_species^3 = 11664, which is probably tiny compared to other things.  Since this is constant for all grid cells at all times, it should be precalculated, using precalculate() prior to evolve
    cdef array.array precalc

    # precalcp is used in evolve.  precalcp[mM, mF, m, gM, gF, g, s] = offspring_modifier[s, mM, mF] * hybridisation[mM, mF, m] * emergence_rate[mF] * inheritance_cube[gM, gF, g] * fecundity[gM, gF, s] * reduction[gM, gF].  It has index s + num_sexes * (g + num_genotypes * (gF + num_genotypes * (gM + num_genotypes * (m + num_species * (mF + num_species * mM))))), so is of size num_sexes * num_genotypes^3 * num_species^3 = 11664, which is probably tiny compared to other things.  Since this is constant for all grid cells at all times, it should be precalculated, using precalculate() prior to evolve
    cdef array.array precalcp

    cdef float min_cc
    
    cpdef void setMinCarryingCapacity(self, float value)
    """Sets min_cc to value.  If the carrying capacity is less than this value, no new larvae are produced"""

    cpdef float getMinCarryingCapacity(self)
    """Gets the value of min_cc"""

cdef class CellDynamicsMosquitoBH26Delay(CellDynamics26DelayBase):
    """Solves Mosquito ODE with 2 sexes, 6 genotypes, using a Beverton-Holt delay equation.

    Carrying capacity is wrapped up in the Beverton-Holt "qm" parameter, as documented in TODO (Hence, min_cc is not used).

    An example of a delay equation is
    dx/dt = x(t - delay * dt)
    Here "delay" is the number of dt involved in the delay equation
    So, if delay = 0 then the equations just revert to ODEs.

    The number of populations is
    (delay + 1) * num_species * num_geotypes * num_sexes = (delay + 1) * num_speces * 12.

    An important variable is current_index, which is an integer between 0 and delay.
    This defines where the current populations are located in pops_and_params.
    It is initialised to 0 in the constructor.
    Examples:
      - Adults of species M, genotype G, sex S have index = M + G * num_species + S * num_species * num_genotypes + current_index * num_species * num_genotypes * num_sexes
      - The adult population at the previous dt has index = M + G * num_species + S * num_species * num_genotypes + (current_index - 1)%(delay + 1) * num_species * num_genotypes * num_sexes
      - The adult population at N dt-steps ago has index = M + G * num_species + S * num_species * num_genotypes + (current_index - N)%(delay + 1) * num_species * num_genotypes * num_sexes
      - The population at delay dt-steps ago has index = M + G * num_species + S * num_species * num_genotypes + (current_index - delay)%(delay + 1) * num_species * num_genotypes * num_sexes = = M + G * num_species + S * num_species * num_genotypes + (current_index + 1)%(delay + 1) * num_species * num_genotypes * num_sexes

    Hence, for the simple process equation of the form dx/dt = x(t - delay * dt), evolve could simply set:
    delayed_index = (current_index + 1)%(delay + 1)
    new_pop = pops_and_params[current_index] + dt * pops_and_params[delayed_index]
    pops_and_params[delayed_index] = new_pop

    !!! IMPORTANT NOTE !!!
    After evolve is called for all grid cells, current_index must be incremented!!  This is the responsibility of the calling code, for the cellDynamics class has no way of knowing evolve has been called for all grid cells.  (An alternative would be to make current_index one of the spatially-varying parameters, and evolve could update it at every call of evolve, but this would consume memory.)
    
    Sex 0 is male.  Sex 1 is female.
    Genotypes:
      0 is ww
      1 is wc
      2 is wr
      3 is cc
      4 is cr
      5 is rr
    There are num_species spatially-varying parameters, which are the carrying-capacities of each species
    """
    # number of genotypes to be calculated in calcXprimeM, calcYYprime and calcCompetition.  If there are only wild-type (that is, genotype=0) mosquitoes present, then by setting num_genotypes_to_calc=1, eliminates loops over the genotype (only genotype=0 is used).  Otherwise, num_genotypes_to_calc should be set to num_genotypes.  num_genotypes_to_calc is currently only used when computing qm from carrying capacities."""
    cdef unsigned num_genotypes_to_calc

    # number of sexes to be calculated in calcYYprime and calcCompetition.  If both sexes have identical properties, including identical populations (this is probably only relevant for purely wild-type populations that are at equilibrium at their carrying capacity) then computational time may be saved by setting num_sexes_to_calc=1, which eliminates loops over the sexes (only sex=0 is used) in the aforementioned functions.  Otherwise, num_sexes_to_calc should be set to 2.  num_sexes_to_calc is currently only used when computing qm from carrying capacities."""
    cdef unsigned num_sexes_to_calc

    # determines whether the parameters in pops_and_params represent the carrying capacity or q_m."""
    cdef unsigned use_qm

    # defines whether a discrete approach is used (Geoff method) or a continuous approximation (Andy method)
    cdef unsigned Geoff_method

    # Array to indicate presence or absence of species, sized to num_species
    cdef array.array species_present

    # Array to indicate presence or absence of genotype, sized to num_genotypes
    cdef array.array genotype_present

    # this is used in evolve to hold the new population, new_pop[s, g, m] with index m + g * self.num_species + s * self.num_species * self.num_genotypes
    cdef array.array new_pop

    # this is used in evolve to hold xprimeM (proportionate mixing quantity).  xprimeM[gM][mM][mF] has index mF + mM * num_species + gM * num_species * num_species
    cdef array.array xprimeM

    # Y used in evolve.  yy[sex][genotype][m] has index m + genotype * num_species + sex * num_species * num_genotypes
    cdef array.array yy

    # Y' used in evolve.  yyp[sex][genotype][m] has index m + genotype * num_species + sex * num_species * num_genotypes
    cdef array.array yyp

    # comp used in evolve.  comp[m]
    cdef array.array comp

    # birth_terms passed to output via evolveCells.  birth_terms[sex][genotype][m] has index m + genotype * num_species + sex * num_species * num_genotypes
    cdef array.array birth_terms

    # death_terms used in calcQm.  death_terms[m]
    cdef array.array death_terms

    # qm_vals used in calcQm.  qm_vals[m]
    cdef array.array qm_vals

    # precalc is used in evolve.  precalc[mM, mF, m, gM, gF, g, s] = hybridisation[mM, mF, m] * emergence_rate[mF] * inheritance_cube[gM, gF, g] * fecundity[gM, gF, s] * reduction[gM, gF].  It has index s + num_sexes * (g + num_genotypes * (gF + num_genotypes * (gM + num_genotypes * (m + num_species * (mF + num_species * mM))))), so is of size num_sexes * num_genotypes^3 * num_species^3 = 11664, which is probably tiny compared to other things.  Since this is constant for all grid cells at all times, it should be precalculated, using precalculate() prior to evolve
    cdef array.array precalc

    # precalcp is used in evolve.  precalcp[mM, mF, m, gM, gF, g, s] = offspring_modifier[s, mM, mF] * hybridisation[mM, mF, m] * emergence_rate[mF] * inheritance_cube[gM, gF, g] * fecundity[gM, gF, s] * reduction[gM, gF].  It has index s + num_sexes * (g + num_genotypes * (gF + num_genotypes * (gM + num_genotypes * (m + num_species * (mF + num_species * mM))))), so is of size num_sexes * num_genotypes^3 * num_species^3 = 11664, which is probably tiny compared to other things.  Since this is constant for all grid cells at all times, it should be precalculated, using precalculate() prior to evolve
    cdef array.array precalcp

    # Used in evolve if use_qm = 0.  This array stores the equilibrium populations and parameters
    cdef array.array eqm_pops_and_params

    cpdef unsigned calcXprimeM(self, unsigned delayed_base, float[:] pops_and_params)
    """Calculates X'_M (= self.xprimeM), given the pops_and_params, and the delayed_base, which usually = (self.current_index + 1) % (self.delay + 1) * self.num_species * self.num_genotypes * self.num_sexes.  Returns 0 if there are no male mosquitoes whatsoever in the delayed pops_and_params slots, and returns 1 otherwise.  This is used in evolve() and probably isn't much use elsewhere"""

    cpdef void calcYYprime(self, unsigned some_males, unsigned delayed_base, float[:] pops_and_params)
    """Calculates Y (= self.yy) and Y' (= self.yyp), given some_males (whether there are any males whatsoever in the delayed pops_and_params), and pops_and_params, and the delayed_base, which usually = (self.current_index + 1) % (self.delay + 1) * self.num_species * self.num_genotypes * self.num_sexes.  This is used in evolve() and calcQm() and probably isn't much use elsewhere.  Note that if num_sexes_to_calc!=2 and num_genotypes_to_calc!=num_genotypes, then Y and Y' will not be fully populated (this occurs during calcQm, for instance)"""

    cpdef void calcCompetition(self, unsigned some_males)
    """Calculates C = competition (= self.comp), given some_males (whether there are any males whatsoever in the delayed pops_and_params).  This uses Y internally, so calcYYprime() should usually precede the call to this function.  This is used in evolve() and probably isn't much use elsewhere"""

    cpdef setNumGenotypesToCalc(self, unsigned num_genotypes_to_calc)
    """Sets number of genotypes to be calculated in calcXprimeM, calcYYprime and calcCompetition. Saves computation when calculating these with no genotypes and equal sexes as in carrying capacity."""

    cpdef unsigned getNumGenotypesToCalc(self)
    """Gets the value of num_genotypes_to_calc"""

    cpdef setNumSexesToCalc(self, unsigned num_sexes_to_calc)
    """Sets number of sexes to be calculated in calcXprimeM, calcYYprime and calcCompetition. Saves computation when calculating these with no genotypes and equal sexes as in carrying capacity."""
    
    cpdef unsigned getNumSexesToCalc(self)
    """Gets the value of num_sexes_to_calc"""

    cpdef setUseQm(self, unsigned use_qm)
    """Sets the use_qm constant to determine whether the parameters in pops_and_params represent the carrying capacity or q_m."""
    
    cpdef unsigned getUseQm(self)
    """Gets the value of use_qm"""

    cpdef setGeoffMethod(self, unsigned Geoff_method)
    """Sets the Geoff_method constant to determine whether Geoff's method is used """
    
    cpdef unsigned getGeoffMethod(self)
    """Gets the value of Geoff_method"""

    cpdef array.array getYYprime(self)
    """Gets the value of yyp"""
    
    cpdef array.array getB(self)
    """Gets the value of B"""    
    