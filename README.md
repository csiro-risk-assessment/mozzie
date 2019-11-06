# mozzie
Mosquito lifecycle, diffusion and advection

## Use

The core code is written in `cython`, which is a mix of python (ease of development) and C (performance).  To our knowledge, all python distributions come 
with `cythonize` which converts the cython code to C code, which may then be compiled and run.  The actual process of doing this is different on different computers: in the `code` directory 
we provide a few different build scripts (`build_easy.sh`, `build_pearcey.sh`, etc).

## Units

The units used should be consistent throughout.  For instance, if the spatial grid is defined using km, and the time-step is measured in days, then the diffusivity should have units km*km/day, and the wind velocity should be in km/day.

## Spatial structure

Mosquitoes are assumed to advect and diffuse over a grid of square cells, defined by:
- (xmin, ymin): the lower left-hand corner
- the cell side-length
- (nx, ny): the number of cells in the x and y directions

The cells can be "active" or "inactive".  This is specified through a CSV file.  The CSV file must contain a header (that begins with `#xmin`, etc) that specifies the quantities mentioned 
above (this facilitates error-checking) and data arranged in rows.  Here is an example file:
```
#xmin=1.0,ymin=2.0,cell_size=3.0,nx=4,ny=3
1,0,1,1
1,0,1,0
0,0,1,1
```
The lines following the header correspond to rows of cells.  The rows appear in *upside-down order*, viz, in the above example:
- the first line, `1,0,1,1`, corresponds to the cells at `y=ymin=2.0`.  The cell at `xmin=1.0` is active; the cell at `xmin=4.0` is inactive; the cell at `xmin=7.0` is active; the cell at `xmin=10.0` is active.
- the second line, `1,0,1,0`, corresponds to the cells at `y=ymin+cell_size=5.0`.  The cell at `xmin=1.0` is active; the cell at `xmin=4.0` is inactive; the cell at `xmin=7.0` is active; the cell at `xmin=10.0` is inactive.
- the third line, `0,0,1,1`, corresponds to the cells at `y=ymin+2*cell_size=8.0`.  The cell at `xmin=1.0` is inactive; the cell at `xmin=4.0` is inactive; the cell at `xmin=7.0` is active; the cell at `xmin=10.0` is active.


## Directory layout

- `tests` directory contains tests (.py files) and associated files (all other files).  Run the tests by using, for example, `python TestGrid.py`
- `code` directory contains the core code that numerically simulates mosquito population dynamics
- `code/auxillary` directory contains python scripts that perform auxillary functions, such as plotting results.  These scripts depend on lots of python libraries, and, while useful, are not necessary for the numerical simulation of mosquitoes.


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
  - If all mosquitoes stay in the wind for 0.5 days, BUT an internal timestep of 0.1 days is required to accurately track advective movements, then pdf = [[0.1, 0], [0.2, 0], [0.3, 0], [0.4, 0], [0.5, 1.0]]
- `grid`: a `Grid` instance.  This is to facilitate error checking (eg, that the wind is defined on the same grid) and to retrieve active cells

A `Wind` object is only of use if `getProcessedDataComputed()==1`.  Upon construction, this is not true.  To make this true, one (or both) of the following methods must be called.
- `parseRawFile()`.  This parses the raw file defined in the constructor, to extract the velocity at each cell.  This data is then processed in the following way.
  - For each active cell, a mosquito is advected by that cell's velocity for `time0` (specified in the `pdf`).  Whichever cell it reaches (can be the same cell), there is a probability `prob0` of it exiting the airstream to that cell.  If that cell is active, that cell's ID and the probability is recorded.
  - Then it is advected by that new cell's velocity for time `time1 - time0`.  It may end up in another new cell (or the same cell), and if that cell is active, information is recorded again.
  - This process is repeated until the `pdf` is exhausted
  - This process is repeated for all active cells.
  - This has the consequence that if a mosquito ends its advection in a inactive cell it is assumed to instantly die, but mosquitoes may in principal pass over inactive cells, using the wind in those cells to reach active cells.  If a mosquito ever advects outside the grid's boundary they are assumed to instantly die.
  - The result may be written to a file using `outputProcessedCSV()`.
- `parseProcessedFile()`.  This parses the processed file defined in the constructor, to build the data mentioned in the previous bullet-point.

After this, the `Wind` object is fully operational, and the useful methods are:
- `getAdvectionFrom()`, returns the array `f`.
- `getAdvectionTo()`, returns the array `t`.
- `getAdvectionP()`, returns the array `p`.
- `getNumAdvection()`, returns the length of `f` (which equals the length of `t`, which also equals the length of `p`).

For any `i`, `f[i]` is the active cell index from which a mosquito is advecting, `t[i]` is the active cell index to which it is advecting, and `p[i]` is the probability of this occuring.  This is stored in the processed data file and may be outputted using outputProcessedCSV().

Note that `f` may not contain all active cell indices.  For instance, for a cell on the grid boundary, wind may instantly advect all mosquitoes out of the domain.  Hence, `p` is the probability of moving from `f` to `t`, given that the mosquito is indeed advecting, and not simply the probability of advecting away from `f`.  This latter probability is specified elsewhere in the code.

