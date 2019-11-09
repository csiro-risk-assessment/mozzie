# mozzie
Mosquito lifecycle, diffusion and advection

## Use

### Building the code

The core code is written in `cython`, which is a mix of python (ease of development) and C (performance).  To our knowledge, all python distributions come 
with `cythonize` which converts the cython code to C code, which may then be compiled and run.  The actual process of doing this is different on different computers: in the `code` directory 
we provide a few different build scripts (`build_easy.sh`, `build_pearcey.sh`, etc).

### Simulating

The core code consists of python objects that you must instantiate in a "runner" python script that defines your mathematical model.  An example is `example1/runner.py`.  Generally, each "runner" python script contains the following.

- An `import` block in which you import all the core libraries (described below) and any other python libraries you need
- A block in which you set up the grid and active cells, defining the spatial extents and discretisation of the model.  You will use `Grid` and `Grid.setActiveAndInactive`.
- A block that defines `Wind` (multiple input files corresponding to different times).  This is different from the spatially-varying parameters (next item) because it usually requires some sort of pre-processing, or is read from specially-preprocessed files.
- A block in which you read files corresponding to spatially-varying parameters in your model, such as carrying capacities.  You'll typically have to read multiple of these per parameter, as the parameters will also vary with time.  You will use `SpatialDependence` and `SpatialDependence.restrictToActive`
- A definition of your cell dynamics (the ODE model), such as `CellDynamicsLogistic1_1`
- A block that defines the populations and parameters, and sets initial conditions.  You'll use `PopulationsAndParameters`.
- A block that gathers the grid structure, the populations and the cell dynamics together into a `SpatialDynamics` object.
- A block that controls the time-stepping of the simulation, involving `diffuse`, `evolveCells`, and `advect`.  There can be multiple of these, depending on the complexity of your model, as different carrying capacities, wind vectors, etc, are used at different times.

Because of this "block" structure, you can run partial simulations, for instance, just processing wind files, or loading files, processing in some way and outputting to produce figures.

### Directory layout

- `tests` directory contains tests (.py files) and associated files (all other files).  Run the tests by using, for example, `python TestGrid.py`.  Much of the code is tested by multiple tests.
- `code` directory contains the core code (described below) that numerically simulates mosquito population dynamics.
- `code/auxillary` directory contains python scripts that perform auxillary functions, such as plotting results.  These scripts are useful but are not necessary for the numerical simulation of mosquitoes.


## Units

The units used should be consistent throughout.  For instance, if the spatial grid is defined using km, and the time-step is measured in days, then the diffusivity should have units km*km/day, and the wind velocity should be in km/day.

## Spatial structure and spatially-varying quantities

Mosquitoes are assumed to advect and diffuse over a grid of square cells, defined by:
- (xmin, ymin): the lower left-hand corner
- the cell side-length
- (nx, ny): the number of cells in the x and y directions

The cells can be "active" or "inactive".  This is specified through a CSV file.  Indeed, all parameters that vary spatially, such as wind vectors, or carrying capacities, or outputted populations numbers, are held in similarly-formatted CSV files.  All such CSV files must contain a header that begins with `#xmin=...` that specifies the quantities mentioned 
above.  This facilitates error checking.  Some files must contain other headers (defined below) and the code will complain if the headers are incorrect.  Data is arranged in rows and must be defined over the entire grid, not just the active cells (data on the inactive cells is never used, however, except for wind data).

Here is an example "active/inactive" file:
```
#xmin=1.0,ymin=2.0,cell_size=3.0,nx=4,ny=3
1,0,1,1
1,0,1,0
0,0,1,1
```
The lines following the header lines (preceeded by "#") always correspond to rows of cells.  The rows always appear in *upside-down order*, viz, in the above example:
- the first line, `1,0,1,1`, corresponds to the cells at `y=ymin=2.0`.  The cell at `xmin=1.0` is active; the cell at `xmin=4.0` is inactive; the cell at `xmin=7.0` is active; the cell at `xmin=10.0` is active.
- the second line, `1,0,1,0`, corresponds to the cells at `y=ymin+cell_size=5.0`.  The cell at `xmin=1.0` is active; the cell at `xmin=4.0` is inactive; the cell at `xmin=7.0` is active; the cell at `xmin=10.0` is inactive.
- the third line, `0,0,1,1`, corresponds to the cells at `y=ymin+2*cell_size=8.0`.  The cell at `xmin=1.0` is inactive; the cell at `xmin=4.0` is inactive; the cell at `xmin=7.0` is active; the cell at `xmin=10.0` is active.

You may include extra headers, preceded by "#", and they will be ignored.


## Core code descriptions

### `Grid`

This is a cython class that may be imported or cimported into other classes.  It contains geometric information (`xmin`, `ymin`, `cell_size`, `nx` and `ny`), but its main purpose is to handle the following.

- Active and inactive cells.  All other classes need only work with active cells.  These methods are useful:
  - `setActiveAndInactive(filename)` sets active and inactive cells with file format specified in the "Spatial structure" section, above.
  - `getNumActiveCells()` returns the number of active cells
  - `getGlobalIndex()` returns `a`, where `a[i]` is the global index of the active cell `i`
  - `getActiveIndex()` returns `a`, where `a[i]` is the active cell index of the global cell `i`
- Adjacency information.  This may be retrieved from two methods: `getConnetionsFrom()` and `getConnectionsTo()`, which both return arrays.  Call the arrays `f` and `t`.  Then `f[i]` is the active cell index, and it is connected to the active cell with index `t[i]`.  So this is a sparse representation of a symmetric adjacency matrix.  Only active cells are considered, so other classes need not know about inactive cells.

`Grid` also contains various other utility methods, and an `outputActiveCSV(filename)` method.

### `SpatialDependence`

This is a cython class that may be imported or cimported into other classes.  It's purpose is to load and hold spatially-dependent parameters, such as carrying capacity.  Most simulations would involve many of these objects.  It is instantiated with `xmin`, etc, information, which facilitates error checking in its other methods.  The methods are:

- `parse`, which parses a CSV file, which typically has format described in the "Spatial structure" section above.

- `restrictToActive` which restricts the data read by `parse` to active cells only.  After this, the data may be used

- `getData0` provides a pointer to the data (in certain cases `getData1` and `getData2` provide pointers to additional data).

- `outputCSV` outputs the data to a CSV file


### `Wind`

This is a cython class that may be imported or cimported into other classes.  Its purpose is to define wind advection, for a single time-step over the entire spatial grid defined by a `Grid` object.  If wind advection does not change in time, just one `Wind` object is needed, while if wind advection changes daily or seasonally, multiple `Wind` objects will be needed.

A `Wind` object is instantiated with:
- `raw_wind_fn`: a filename describing the raw wind data.  This must be defined over the entire grid, not just at the active cells.  These assumptions are made:
  - the file contains a header of the form `#xmin=...`, exactly the same as defined in the "Spatial structure" section, above.  This facilitates error checking.
  - the file contains `ny` lines of data.  Each line contains CSV entries defining the wind velocities in the x and y direction.  The line has the format `vx0,vy0,vx1,vy1,vx2,vy2,...`, where `vx0` is the x-component of velocity at the first cell, etc.  So, each line contains `2*nx` pieces of data.  The first line corresponds to the velocities at `ymin`, the second to velocities at `ymin + cell_size`, etc.
  - the units for velocity must be consistent with the remainder of the simulation (for instance km/day)
- `processed_wind_fn`: a filename describing the processed wind data.  Each data line contains three numbers: `from_index,to_index,probability`.  The method `outputProcessedCSV()` writes this file for later (or immediate) use.  Examples of files are in the test suite, and some more details are mentioned below (see `getAdvectionFrom()`, etc).
- `pdf`: a probability distribution of times an advecting mosquito will stay in the air and be advected by the wind.  This is specified as a list of the form `[[time0, prob0], [time1, prob1], [time2, prob2], ...]`.  The times must be in the same time unit as the wind data (for instance, days).  The probabilities should sum to 1 (prob0 + prob1 + ... = 1).  Some examples are:
  - If all mosquitoes that are being advected by the wind stay in the wind for 0.5 days, then pdf = [[0.5, 1.0]]
  - If 30% of mosquitoes stay in the wind for 0.1 days, and 70% stay in the wind for 0.4 days, then pdf = [[0.1, 0.3], [0.4, 0.7]]
  - If all mosquitoes stay in the wind for 0.5 days, BUT an internal timestep of 0.1 days is required to accurately track advective movements, then pdf = [[0.1, 0], [0.2, 0], [0.3, 0], [0.4, 0], [0.5, 1.0]].  See the discussion in `parseRawFile()` to understand this better.
It is important to make the `pdf` commensurate with the time-step size used in `SpatialDynamics.diffuse` and `SpatialDynamics.evolveCells` described below.
- `grid`: a `Grid` instance.  This is to facilitate error checking (eg, that the wind is defined on the same grid) and to retrieve active cells

A `Wind` object is only of use if `getProcessedDataComputed()==1`.  Upon construction, this is not true.  To make this true, one (or both) of the following methods must be called.
- `parseRawFile()`.  This parses the raw file defined in the constructor, to extract the velocity at each cell.  This data is then processed in the following way.
  - For each active cell, a mosquito is advected by that cell's velocity for `time0` (specified in the `pdf`).  Whichever cell it reaches (can be the same cell), there is a probability `prob0` of it exiting the airstream to that cell.  If that cell is active, that cell's ID and the probability is recorded.
  - Then it is advected by that new cell's velocity for time `time1 - time0`.  It may end up in another new cell (or the same cell), and if that cell is active, information is recorded again.
  - This process is repeated until the `pdf` is exhausted.
  - This process is repeated for all active cells.
  - This has the consequence that if a mosquito ends its advection in a inactive cell it is assumed to instantly die, but mosquitoes may in principal pass over inactive cells, using the wind in those cells to reach active cells.  If a mosquito ever advects outside the grid's boundary they are assumed to instantly die.
  - The result may be written to a file using `outputProcessedCSV()`.
- `parseProcessedFile()`.  This parses the processed file defined in the constructor, to build the data mentioned in the previous bullet-point.

After this, the `Wind` object is fully operational, and the useful methods are:
- `getAdvectionFrom()`, returns the array `f`.
- `getAdvectionTo()`, returns the array `t`.
- `getAdvectionP()`, returns the array `p`.
- `getNumAdvection()`, returns the length of `f` (which equals the length of `t`, which also equals the length of `p`).

For any `i`, `f[i]` is the active cell index from which a mosquito is advecting, `t[i]` is the active cell index to which it is advecting, and `p[i]` is the probability of this occuring.  This is stored in the processed data file and may be outputted using `outputProcessedCSV()`.

Note that `f` may not contain all active cell indices.  For instance, for a cell on the grid boundary, wind may instantly advect all mosquitoes out of the domain.  Hence, `p` is the probability of moving from `f` to `t`, given that the mosquito is indeed advecting.  It is not simply the probability of advecting away from `f`.  This latter probability is specified in `spatialDynamics.advect`.

### `CellDynamicsX`

This is a cython class that may be imported or cimported into other classes.  Its purpose is to define the mosquito lifecycle dynamics at the grid-cell level.  Different dynamics may be easily defined by inheriting from the base class (the `X` is replaced by something like `Logistic`).  The following methods are important

- `getNumberOfPopulations` and `getNumberOfParameters` define the number of populations (maleGG, femaleGW, etc) and parameters (carrying capacity, mortality rate, etc) in the ODEs describing the lifecycle dynamics

- `getNumberOfDiffusingPopulations` and `getDiffusingIndices` define which of the populations diffuse over the spatial grid (eg, mosquito eggs do not diffuse, but adults do).

- `getNumberOfAdvectingPopulations` and `getAdvectingIndices` define which of the populations advect over the spatial grid (eg, mosquito eggs do not advect, but adults do).

- `evolve(dt, pp)` is the most important method, for it solves the ODEs.  Here `dt` is the time-step size, and `pp` is an array holding the populations and parameters.  When another class (eg the "runner" of the simulation) calls this method, it will have placed the current populations and parameters into `pp`, and will expect `evolve` to modify the populations so they're correct for the next time-step.

### `PopulationsAndParameters`

This is a cython class that may be imported or cimported into other classes.  Its purpose is to hold the populations and parameters for the entire active grid.  It is constructed using the `Grid` and `CellDynamicsX` so that it knows the number of cells and populations+parameters per cell.  It has some useful `set` methods for setting the populations and parameters, and a `get` method to retrieve these.

### `SpatialDynamics`

This is a cython class that may be imported into other classes.  It provides methods to `diffuse`, `advect` and `evolveCells` the cell populations.  It is constructed using the `Grid` and `PopulationsAndParameters` objects.  Its most important methods are the following

- `diffuse(dt, diffusion_coeff)`.  This performs one time-step of diffusion, updating the grid-cell populations contained in `PopulationsAndParameters`.  Only the populations labelled as "diffusing" by the `CellDynamicsX` object are updated.  `diffuse` uses the nearest-neighbour finite-difference approximation to the Laplacian.  Hence, with time-step size `dt` the fraction of mosquitoes removed from one cell is `diffusion_coeff * 4 * dt / (dx)^2`, where `dx` is the cell-size.  One quarter of this amount is added to each of the 4 neighbours.  If the neighbours happen to be inactive, the mosquitoes are assumed to die instantly.

- `advect(frac, wind)`.  This performs one time-step of advection using the given `wind` object.  The grid-cell populations contained in `PopulationsAndParameters` are updated.  Only the populations labelled as "advecting" by the `CellDynamicsX` object are updated.  The `wind` object will have processed the raw-wind data using its `pdf`, and mosquitoes are advected as described in the `Wind` section above.  For instance, if the `pdf` contains times up to `timeN`, then `advect` will advect mosquitoes up to that time, which is completely independent of the time-step size in `diffuse` or `evolveCells`.

- `evolveCells(dt)`.  This calls `CellDynamics.evolve(dt)` for all cells in the grid.  Thus it performs one lifecycle time-step over the entire grid

- `outputCSV(filename, pop, inactive_val, header)` outputs information for population number `pop` in `PopulationsAndParameters` for the entire grid to the `filename`.
