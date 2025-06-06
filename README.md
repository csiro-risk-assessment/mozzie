# Mosquito lifecycle, diffusion and advection

## Description

`Mozzie` enables simulation of the lifecycle and spatial spread of mosquitoes.  `Mozzie` can be used to assess risks associated with disease-control strategies at local, regional or continental scales.  Strategies involving genetic alterations of mosquitoes to eliminate malaria, are of prime interest.

More technically, `Mozzie` simulates a population-dynamics model that uses differential equations or delay differential equations to describe the spread and persistence of mosquitoes that may be genetic altered.  Genetic alterations are flexibly modelled: these can involve any number of alleles; Mendelian or non-Mendelian inheritance, including gene drives; they can be self-limiting or self-sustaining; and can include the emergence of resistant allelles.  The model allows simulation of any number of mosquito species.   It incorporates mate-choice, hybridisation and intra-specific competition that occur within complexes of mosquito species.  This fills a gap that currently exists among similar models, allowing researchers to assess potential transfer of the genetic alterations between (sub-)species.

`Mozzie` supports spatial and temporal variations in lifecycle parameters, and local diffusion and wind-assisted, long range, advection.  For example, wind patterns and the capacity of the landscape to support mosquitoes can vary spatially and temporally, reflecting daily variations, seasonality, and local conditions.

Conversely, `Mozzie` does not contain human agents, nor does it consider the effect of genetic control strategies on the prevalence of pathogens such as the malaria parasite, among human or animal populations.


## Contributing

Australia's Commonwealth Scientific and Industrial Research Organisation (CSIRO) welcomes contributions to this code.  The recommended method is to:
 
1. Submit an [Issue](https://github.com/csiro-risk-assessment/mozzie/issues) on GitHub, describing the reasons for your potential contribution;
2. Engage in discussions with CSIRO staff on GitHub regarding the best architecture, etc, for your contribution;
3. When your contribution is ready, issue a Pull Request on GitHub.
 
CSIRO does not guarantee to inspect such contributions, nor accept them into the code.


## Licensing and disclaimers

This code is released under a GPLv3 license. To discuss other licenses, please contact CSIRO via GitHub.
 
Copyright (c) 2024 Commonwealth Scientific and Industrial Research Organisation (CSIRO) ABN 41 687 119 230.
 
CSIRO does not make any guarantees of accuracy or relevance of the code.  CSIRO does not guarantee the code has no security vulnerabilities.  CSIRO does not guarantee to bug-fix or security-fix the code.  CSIRO reserves the right to "retire" the code and make no further changes.

## Reporting problems and seeking support

Please report problems using a new GitHub [issue](https://github.com/csiro-risk-assessment/mozzie/issues).  Please tag someone (eg @WilkAndy) in your issue.  CSIRO will endeavour to fix any problems, but does not guarantee to do so.

To seek help, please also use a new GitHub [issue](https://github.com/csiro-risk-assessment/mozzie/issues), and please tag someone (eg @WilkAndy).  Using a GitHub issue means that future users may benefit from the solution to your problem.  CSIRO will endeavour to help, but does not guarantee to do so.


## How to build and test the code

The core code is written in `cython`, which is a mix of python (ease of development) and C (performance).  Your computer system possibly has all the necessary features already installed, but a vanilla system will need various items.

You must be familiar with the terminal (linux/mac) or command prompt (windows) to use the mozzie code.  In the remainder of this document, the python interpreter will be denoted by `python3`, even though the python interpreter may be invoked by `python` or `py`, etc, on your system.

### Step 0

Obtain the `mozzie` code.  This will probably be via a `git clone` command, and probably you have already done this step.

### Step 1a

Ensure your system has python3 with pip and venv, cmake and a C compiler.  You may check by entering the following commands at the terminal:

- `python3 --version`: the version number should be 3.10 or greater
- `pip --version` shouldn't return an error
- `python3 -m venv -h` shouldn't return an error
- `cmake --version` shouldn't return an error.

On Windows computers, replace `python3` with `python` or `py` here and elsewhere.

If one or more of these aren't available, one of the following could help.

- On ubuntu and Mac: `sudo apt install python3-dev python3-pip python3-venv cmake`
- On redhat: TODO
- On an HPC supercomputer: `module load XXXX` (ask your systems administrators for what XXXX should be)
- On Windows:
  - Install python by downloading from the internet (eg [python.org](https://python.org)).  Ensure you add python to your PATH, and install the pip optional package.
  - Install cmake by downloading from cmake.org.  Ensure you add cmake to your PATH.
  - Install the C compiler (eg Visual Studio Build Tools 2022) from [Microsoft's website](https://visualstudio.microsoft.com/visual-cpp-build-tools/), ticking the "Desktop development with C++".

### Step 1b

Create a virtual environment and activate it.  This is so you can `pip install` the necessary packages without conflicting with other things on your system.

On non-windows computers:

```
python3 -m venv ~/mozzie_venv
. ~/mozzie_venv/bin/activate
```

On windows computers using the command prompt:

```
python -m venv mozzie_venv
mozzie_venv\Scripts\activate.bat
```

In the above commands, `mozzie_venv` can be any path you desire.  You should remember it for later use of the `mozzie` software.  **Whenever** you want to work with the `mozzie` software, you should first `. ~/mozzie_venv/bin/activate` (on non-windows computers) or `mozzie_venv\Scripts\activate.bat` (on windows computers using the command prompt).


### Step 1c

Install all the required python packages.  Navigate to the mozzie repository, and then:

```
pip install -r requirements.txt
```

Because you have created and activated your `mozzie_venv`, this will install the required python packages into that virtual environment.  Check all is well by typing `coverage` into your terminal or command prompt.  If an error is returned, you will have to modify your PATH variable.


### Step 2

Compile the code.  Navigate to the mozzie repository, and then

```
cd code
python3 setup.py build_ext --inplace
```

There may be some warnings from your Cython installation.  Don't worry for now, but maybe later you can contribute changes to the code to avoid these warnings!

If you are going to use the code in earnest, you will very likely want to manipulate plaintext and binary input files, for which you need the `ab_convert` program (see below for documentation details).  To create this program, navigate to the mozzie repository, and then

```
cd code/auxillary
cmake .
cmake --build .
```

The code is built in the `code` directory, so when you want to `import` the modules you have to ensure your `PATH` is set correctly.  See the primers and examples below for how we do this.

### Step 3

Test the code.  Navigate to the mozzie repository and then

```
cd tests
coverage run -m unittest -v
coverage report
```

You may get some errors associated with MAX_FILE_LENGTH being too large if you are running on a computer with not enough memory.

You are now ready to start using or developing the code!

## Primers

Before attempting a full-scale simulation with complicated lifecycle dynamics, you might want to inspect the samples found in the `primers` directory.  In order of increasing complexity, these are:

- `logistic.py` which describes how to simulate lifecycle dynamics at a single location with no spatial spread.  In this case, the logistic equation is used.  The output is plotted in `logistic.pdf`.
- `diffusion.py` which describes the diffusion of a single species with no lifecycle dynamics (the species does not breed, evolve or die).  The output is plotted in `diffusion.pdf` and `diffusion_line.pdf`.
- `advection.py` which describes the advection of a single species with no lifecycle dynamics (the species does not breed, evolve or die).  The output is plotted in `advection.pdf`.
- `island.py` which describes lifecycle dynamics, diffusion and advection of a single species.  The lifecycle is governed by the logistic equation.  An island is separated from the mainland.  The species initially exists only on the mainland, but advection carries it to the island.  The output is plotted in `island.pdf`
- `spatial.py` which is the same as `island.py`, but the wind velocity and carrying capacity are spatio-temporally dependent.  The result is plotted in `spatial.pdf`.

Each of these python files contains extensive in-code documentation.  Run each of them using, for instance,

```
cd primers
python3 logistic.py
```

The `cd primers` is necessary because of the way we manipulate `PATH` to pick up the modules (see the start of each of the aforementioned files).


## How to set up a simulation

The core code consists of python objects that you must instantiate in a "runner" python script that defines your mathematical model.  The aforementioned primers contain simple examples, and now we describe `example1/runner.py` that contains all the components of a full-scale simulation.  There are other more sophisticated models in the other `example*` directories, and in the demo directory.  Run the simulations using, for example,

```
cd example1
python3 runner.py
```

### runner.py: the `import` block

The `import` block is used to import all the core libraries (described below) and any other python libraries you need.  `example1/runner.py` looks like:

```
import os
import sys
import array
import timeit

# so we can find our ../code no matter how we are called
findbin = os.path.dirname(os.path.realpath(sys.argv[0]))
sys.path.append(findbin + "/../code")

from wind import Wind
from grid import Grid
from cellDynamics import CellDynamicsLogistic1_1
from spatialDynamics import SpatialDynamics
from spatialDependence import SpatialDependence
from populationsAndParameters import PopulationsAndParameters
```

The first imports are standard python libraries.  In order of import, the others will enable: a definition of wind direction; a definition of the spatial grid; lifecycle dynamics governed by the logistic equation; spatial dynamics such as diffusion and wind advection; spatially-dependent lifecycle parameters such as the carrying capacity; populations and parameter quantities.  These are discussed in more detail below.

### runner.py: the grid block

This is a block of code in which you set up the grid and active cells, defining the spatial extents and discretisation of the model.  You will use `Grid` and `Grid.setActiveAndInactive`.  In `example1/runner.py` this looks like:

```
g1 = Grid(-4614.0, -3967.0, 5.0, 1517, 1667, False)
g1.setActiveAndInactive("active.csv")
```

These lines set up a spatial grid of 5km x 5km, with 1517 x 1667 grid cells in the horizontal x vertical direction, with corner coordinate (-4614, -3967).  Then, some of these cells are made inactive (no mosquitoes will live there: in this case because it an ocean area) using the `setActiveAndInactive` method.

### runner.py: the spatially-dependent lifecycle parameters

In `example1/runner.py`, the carrying-capacity, `cc`, is spatially-dependent:

```
cc_parser = SpatialDependence(-4614.0, -3967.0, 5.0, 1517, 1667)
cc_parser.parse("carrying.csv", "generic_float", [])
cc_parser.restrictToActive(g1.getGlobalIndex())
cc = cc_parser.getData0()
```

Generally, you'll typically have to read multiple of these per lifecycle parameter, as the parameters will also vary with time.  Each time, you will use `SpatialDependence`, `SpatialDependence.restrictToActive`, and `getData0`.

### runner.py: lifecycle definition

This is typically one line.  In `example1/runner.py`, we want to use a logistic-growth model, so it is:

```
cell = CellDynamicsLogistic1_1()
```

### runner.py: setting populations, lifecycle parameters and initial conditions

First, define an object that contains all the mosquito populations over the grid, and all the lifecycle parameters (such as carrying capacity) associated with those populations:

```
all_pops = PopulationsAndParameters(g1, cell)
```

In this case, when the lifecycle parameter (carrying capacity) is time-independent, it is convenient to set it now.  In `CellDynamicsLogistic1_1`, the second slot in the populations and parameters array contains the carrying capacity at the cell, so:

```
pop_and_param_array = all_pops.getQuantities()
for i in range(g1.getNumActiveCells()):
   pop_and_param_array[2 * i + 1] = max(1.0, cc[i])
```

(The `max` ensures there is a nonzero carrying capacity everywhere).  If carrying capacity were time-dependent then `pop_and_param_array` would be set during time-stepping, not during this initialisation phase (the primers contain examples of this).

There are various ways of defining initial conditions, and the `example*/runner*` provide examples. In this case, 10000 mosquitoes are introduced at (x, y) = (-2000, 700), and the carrying capacity is set to 150000 at that cell:

```
all_pops.setPopulationAndParametersFromXY(-2000, 700, [10000.0, 150000.0])
```

### runner.py: readying for simulation

The populations and the cell dynamics together into a `SpatialDynamics` object, to be ready for lifecycle dynamics, diffusion and advection:

```
spatial = SpatialDynamics(g1, all_pops)
```

### runner.py: defining wind

Wind advection is special in `mozzie`.  It is always spatially dependent, and usually temporally dependent.  It is not deemed a "parameter" such as carrying capacity because it is independent of the lifecycle dynamics.  The `mozzie` Wind class is described in detail below: in `example1/runner.py`, the definition is:

```
import math
# define the probability function:
xa = [0.02, 0.04, 0.06, 0.08, 0.1, 0.12, 0.14, 0.17, 0.2, 0.25, 0.3, 0.35, 0.4, 0.45, 0.5]
ya = [math.exp(-x * 6) for x in xa]
su = sum(ya)
ya = [y / sum(ya) for y in ya]
pdf = [[xa[i], ya[i]] for i in range(len(xa))]
# read the wind velocity vectors, and specify the desired proability function to use
wind = Wind(os.path.join(findbin, "raw_wind.csv"), os.path.join(findbin, "wind_8_subdivisions.csv"), pdf, g1)
# Parse the wind velocity vectors, and do particle tracking based on the probability function
wind.parseRawFile()
```

Typically, you'll read multiple input files corresponding to different times.

### runner.py: time-stepping the simulation

After all this set-up, finally we come to the simulation!  Often, there are multiple blocks similar to the following, as different carrying capacities, wind vectors, etc are used at different times.  But they will all use `diffuse`, `evolveCells`, and `advect`.  The primers and `example*/runner*.py` contain examples.  For `example1/runner.py` the block is simple:

```
timestep_size = 0.5
for i in range(501):
   sys.stdout.write("Time step " + str(i) + "\n")
   # lifecycle (logistically grow) using evolveCells for 0.5 days
   spatial.evolveCells(0.5 * 10) # the multiplication here is just because Logistic1_1 has a hard-coded 0.01 growth rate which is too tiny to notice anything interesting
   # now diffuse with diffusion coefficient 0.1 km^2/day
   spatial.diffuse(timestep_size, 0.1)
   # anow dvect 1% of the population at each cell using the defined wind
   spatial.advect(array.array('f', [0.01]), wind)
   if (i%100 == 0):
      spatial.outputCSV("runner_" + str(i) + "_days.csv", 0, "0", "")
```

Note that you will have to **think carefully about the biological reality of your time-stepping**.  In the case above, logistic growth occurs only for half a day (during the daylight), and then diffusion and wind advection occur (eg, during the night).

### runner.py: plotting results and general comments

Notice the `spatial.outputCSV` line above.  This outputs populations to a CSV file (the arguments depend on the lifecycle chosen, such as logistic growth, or any of the other CellDynamics classes).  Matplotlib can be used to plot the results.  The primers mentioned above contain examples of this.

Because of the "block" structure of each `runner.py` script, you can run partial simulations, for instance, just processing wind files, or loading files, processing in some way and outputting to produce figures.  You can also chain together multiple simulations that rely on just one initial block of file reading, to avoid reading data files for each and every simulation.

## Directory layout

- `tests` directory contains tests (.py files) and associated files (all other files).  Run the tests by using, for example, `python TestGrid.py`.  Much of the code is tested by multiple tests.
- `code` directory contains the core code (described below) that numerically simulates mosquito population dynamics.
- `code/auxillary` directory contains python scripts that perform auxiliary functions, such as plotting results and the important `ab_convert` program that converts between plaintext and binary data files.  These scripts are useful but are not necessary for the numerical simulation of mosquitoes.

## File I/O and memory requirements

Full-sized simulations require substantial file I/O and memory.

For this reason, you will typically want to read *binary* versions of the input files (that define wind, carrying-capacity, etc) rather than plaintext files.  The only exception to this is the file that sets inactive/active cells: that must be plaintext.  Reading a binary file using the optimised C parser (this is called `csvparser` in the code) is approximately *30 times* faster than reading a plaintext CSV file using python (`SpatialDependence.parse` vs `SpatialDependence.parseWithPython`).

To instruct the code to read binary, rather than plaintext files, you need to ensure the `SpatialDependence` file-type is `generic_float_binary` (rather than just `generic_float`).  Similarly, to instruct the code that your wind files ("raw" or "processed") are binary, use `Wind.setBinaryFileFormat(1)`.  The `tests` directory has examples of this.

To create a binary version of a plaintext file use the `ab_convert` program contained in the `auxillary` directory.  `ab_convert` stands for "ascii-binary converter".  Eg
```
./code/auxillary/ab_convert ascii2binary 1517 1667 generic_float plaintext.csv binary_version.bin
```
The `ab_convert` program has in-built documentation: just use `./code/auxillary/ab_convert` without arguments to retrieve it.  `ab_convert` can also convert binary files to plain-text.  For instance:
```
./code/auxillary/ab_convert binary2ascii 1517 1667 wind_raw wind_file.bin wind_file.csv
```

The other choice that must be made regarding file I/O is whether to read "raw" wind files (describing the velocity at each grid cell) or "processed" wind files (describing the probability of a mosquito advecting between cells).  Optimally, these are binary.

- The "raw" files are typically much smaller, so the file I/O is much quicker, but they require in-code processing to make them usable.  That means that reading+processing binary files is only about 2.7-times faster than plaintext versions.
- The "processed" files are generated in a pre-processing script by parsing "raw" files (`Wind.parseRawFile`) and then outputting the processed result (`Wind.outputProcessedCSV`).  The `tests` have examples of this.  Reading binary "processed" files is much faster than plaintext versions, but it may be better to read+process "raw" files instead because the "processed" versions are typically substantially bigger (if a detailed `pdf` is used).

It is simply a matter of experimentation to determine whether reading "raw" or "processed" files is faster.

Remember that there are hardware limitations when considering file I/O.  For example, in our experiments on our HPC system, the SSDs have a read-speed of about 1GB/s (this is impacted by caching and other users).  Processed, binary wind data can be read, checked and organised into data structures at a speed of about 0.7GB/s, indicating the code is close to optimal.  `Generic_float_binary` data (describing carrying-capacity, etc) can be read and organised into data structures at greater than 0.8GB/s.

Lots of memory is needed when reading lots of data describing wind, etc.  A lower bound on the program's memory requirements is the size of the binary files (sum of the "processed" wind and the `generic_float_binary` files).  For instance, in a recent simulation, 1 year of processed wind files used 84GB on disk, so the simulation is going to use at least this amount of memory, assuming all the wind files (either "raw" or "processed") are read and stored.  In reality, in this case this simulation uses almost exactly 84GB.


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
The lines following the header lines (preceded by "#") always correspond to rows of cells.  The rows always appear in *upside-down order*, viz, in the above example:
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

- `parse`, which parses a CSV file (or a binary version of it) which typically has format described in the "Spatial structure" section above.  One argument to `parse` is the file type, which is most frequently `generic_float` or `generic_float_binary`.  (In addition, `SpatialDependence` is also used by `Grid` to define inactive/active cells, and `Wind` to define wind vectors, but these are not of the `generic_float` type.)  A `binary` file will have been created by the program `ab_convert`.  (See above for a discussion of binary files: they can be read about 30-times faster than plaintext files.)

- `restrictToActive` which restricts the data read by `parse` to active cells only.  After this, the data may be used.

- `getData0` provides a pointer to the data (in certain cases `getData1` and `getData2` provide pointers to additional data).

- `outputCSV` outputs the data to a CSV file


### `Wind`

This is a cython class that may be imported or cimported into other classes.  Its purpose is to define wind advection, for a single time-step over the entire spatial grid defined by a `Grid` object.  If wind advection does not change in time, just one `Wind` object is needed, while if wind advection changes daily or seasonally, multiple `Wind` objects will be needed.

A `Wind` object is instantiated with:
- `raw_wind_fn`: a filename describing the raw wind data.  This may be in plaintext (default) or binary form (use `setBinaryFileFormat(1)`).  This must be defined over the entire grid, not just at the active cells.  These assumptions are made:
  - the file contains a header of the form `#xmin=...`, exactly the same as defined in the "Spatial structure" section, above.  This facilitates error checking.
  - the file contains `ny` lines of data.  Each line contains CSV entries defining the wind velocities in the x and y direction.  The line has the format `vx0,vy0,vx1,vy1,vx2,vy2,...`, where `vx0` is the x-component of velocity at the first cell, etc.  So, each line contains `2*nx` pieces of data.  The first line corresponds to the velocities at `ymin`, the second to velocities at `ymin + cell_size`, etc.
  - the units for velocity must be consistent with the remainder of the simulation (for instance km/day)
- `processed_wind_fn`: a filename describing the processed wind data.  This may be in plaintext (default) or binary form (use `setBinaryFileFormat(1)`).  Each data line contains three numbers: `from_index,to_index,probability`.  The method `outputProcessedCSV()` writes this file for later (or immediate) use.  Examples of files are in the test suite, and some more details are mentioned below (see `getAdvectionFrom()`, etc).
- `pdf`: a probability distribution of times an advecting mosquito will stay in the air and be advected by the wind.  This is specified as a list of the form `[[time0, prob0], [time1, prob1], [time2, prob2], ...]`.  The times must be in the same time unit as the wind data (for instance, days).  It is important to make the times commensurate with the time-step size used in `SpatialDynamics.diffuse` and `SpatialDynamics.evolveCells` described below.  The probabilities should sum to 1 (prob0 + prob1 + ... = 1).  Some examples are:
  - If all mosquitoes that are being advected by the wind stay in the wind for 0.5 days, then pdf = [[0.5, 1.0]]
  - If 30% of mosquitoes stay in the wind for 0.1 days, and 70% stay in the wind for 0.4 days, then pdf = [[0.1, 0.3], [0.4, 0.7]]
  - If all mosquitoes stay in the wind for 0.5 days, BUT an internal timestep of 0.1 days is required to accurately track advective movements, then pdf = [[0.1, 0], [0.2, 0], [0.3, 0], [0.4, 0], [0.5, 1.0]].  See the discussion in `parseRawFile()` to understand this better.
- `grid`: a `Grid` instance.  This is to facilitate error checking (eg, that the wind is defined on the same grid) and to retrieve active cells

A `Wind` object is only of use if `getProcessedDataComputed()==1`.  Upon construction, this is not true.  To make this true, one (or both) of the following methods must be called.
- `parseRawFile()`.  This parses the raw file defined in the constructor, to extract the velocity at each cell.  This data is then processed in the following way.
  - For each active cell, a mosquito is advected by that cell's velocity for `time0` (specified in the `pdf`).  Whichever cell it reaches (can be the same cell), there is a probability `prob0` of it exiting the airstream to that cell.  If that cell is active, that cell's ID and the probability is recorded.
  - Then it is advected by that new cell's velocity for time `time1 - time0`.  It may end up in another new cell (or the same cell), and if that cell is active, information is recorded again.
  - This process is repeated until the `pdf` is exhausted.
  - This process is repeated for all active cells.
  - This has the consequence that if a mosquito ends its advection in an inactive cell it is assumed to instantly die, but mosquitoes may in principal pass over inactive cells, using the wind in those cells to reach active cells.  If a mosquito ever advects outside the grid's boundary they are assumed to instantly die: there is an "instant-kill zone" outside the grid's boundary
  - The result may be written to a file using `outputProcessedCSV()`.  Usually you will want to `setBinaryFileFormat(1)` before the writing, so a binary file is written.
- `parseProcessedFile()`.  This parses the processed file defined in the constructor, to build the data mentioned in the previous bullet-point.

After this, the `Wind` object is fully operational, and the useful methods are:
- `getAdvectionFrom()`, returns the array `f`.
- `getAdvectionTo()`, returns the array `t`.
- `getAdvectionP()`, returns the array `p`.
- `getNumAdvection()`, returns the length of `f` (which equals the length of `t`, which also equals the length of `p`).

For any `i`, `f[i]` is the active cell index from which a mosquito is advecting, `t[i]` is the active cell index to which it is advecting, and `p[i]` is the probability of this occurring.  This is stored in the processed data file and may be outputted using `outputProcessedCSV()` (as plaintext, or, if `setBinaryFileFormat(1), binary format).

Note that `f` may not contain all active cell indices.  For instance, for a cell on the grid boundary, wind may instantly advect all mosquitoes out of the domain.  Hence, `p` is the probability of moving from `f` to `t`, given that the mosquito is indeed advecting.  It is not simply the probability of advecting away from `f`.  This latter probability is specified in `spatialDynamics.advect`.

### `CellDynamicsX`

This is a cython class that may be imported or cimported into other classes.  Its purpose is to define the mosquito lifecycle dynamics (ODEs) at the grid-cell level.  Different dynamics may be easily defined by inheriting from the base class (the `X` is replaced by something like `Logistic`).  Look at `doc/mosquito23.pdf` for a description of the `CellDynamicsMosquito23` class.

The following methods are important

- `getNumberOfPopulations` and `getNumberOfParameters` define the number of populations (maleGG, femaleGW, etc) and parameters (carrying capacity, mortality rate, etc) in the ODEs describing the lifecycle dynamics.

- `getNumberOfDiffusingPopulations` and `getDiffusingIndices` define which of the populations diffuse over the spatial grid (eg, mosquito eggs do not diffuse, but adults do).

- `getNumberOfAdvectingPopulations` and `getAdvectingIndices` define which of the populations advect over the spatial grid (eg, mosquito eggs do not advect, but adults do).  `setAdvectionClass` and `getAdvectionClass` controls the advection class of the advecting populations: different classes can advect with different probability.

- `evolve(dt, pp)` is the most important method, for it solves the ODEs.  Here `dt` is the time-step size, and `pp` is an array holding the populations and parameters.  When another class (eg `SpatialDynamics.evolveCells`, below) calls this method, it will have placed the current populations and parameters into `pp`, and will expect `evolve` to modify the populations so they're correct for the next time-step.

### `PopulationsAndParameters`

This is a cython class that may be imported or cimported into other classes.  Its purpose is to hold the populations and parameters for the entire active grid.  It is constructed using the `Grid` and `CellDynamicsX` so that it knows the number of cells and populations+parameters per cell.  It has some useful `set` methods for setting the populations and parameters, and a `get` method to retrieve these.

### `SpatialDynamics`

This is a cython class that may be imported into other classes.  It provides methods to `diffuse`, `advect` and `evolveCells` the cell populations.  It is constructed using the `Grid` and `PopulationsAndParameters` objects.  Its most important methods are the following

- `diffuse(dt, diffusion_coeff)`.  This performs one time-step of diffusion, updating the grid-cell populations contained in `PopulationsAndParameters`.  Only the populations labelled as "diffusing" by the `CellDynamicsX` object are updated.  `diffuse` uses the nearest-neighbour finite-difference approximation to the Laplacian.  Hence, with time-step size `dt` the fraction of mosquitoes removed from one cell is `diffusion_coeff * 4 * dt / (dx)^2`, where `dx` is the cell-size.  One quarter of this amount is added to each of the 4 neighbours.  If the neighbours happen to be inactive, the mosquitoes are assumed to die instantly.  If the mosquitoes diffuse outside the grid, they are instantly killed.

- `advect(frac, wind)`.  This performs one time-step of advection using the given `wind` object.  The grid-cell populations contained in `PopulationsAndParameters` are updated.  Only the populations labelled as "advecting" by the `CellDynamicsX` object are updated.  The `Wind` object will have processed the raw-wind data using its `pdf`, and mosquitoes are advected as described in the `Wind` section above.  For instance, if the `pdf` contains times up to `timeN`, then `advect` will advect mosquitoes up to that time, which is completely independent of the time-step size in `diffuse` or `evolveCells`.  This is why the `time` quantities in the `pdf` must be commensurate with `dt`, as mentioned above in the section on the `Wind` class.  The `frac` input to `advect` is an array that defines the fraction of advecting mosquitoes of each advection class that is advected.  For instance, if there are two advection classes (males and females, for instance) then `frac = (0.1, 0.2)` means 10% of the first class advects, while 20% of the second class advects.

- `evolveCells(dt)`.  This calls `CellDynamics.evolve(dt)` for all cells in the grid.  Thus it performs one lifecycle time-step over the entire grid

- `outputCSV(filename, pop, inactive_val, header)` outputs information for population number `pop` in `PopulationsAndParameters` for the entire grid to the `filename`.
