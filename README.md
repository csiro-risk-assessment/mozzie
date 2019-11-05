# mozzie
Mosquito lifecycle, diffusion and advection

## Use

The core code is written in `cython`, which is a mix of python (ease of development) and C (performance).  To our knowledge, all python distributions come with `cythonize` which converts the cython code to C code, which may then be compiled and run.  The actual process of doing this is different on different computers: in the `code` directory we provide a few different build scripts (`build_easy.sh`, `build_pearcey.sh`, etc).

## Spatial structure

Mosquitoes are assumed to advect and diffuse over a grid of square cells, defined by:
- (xmin, ymin): the lower left-hand corner
- the cell side-length
- (nx, ny): the number of cells in the x and y directions

The cells can be "active" or "inactive".  This is specified through a CSV file.  The CSV file must contain a header (that begins with `#xmin`, etc) that specifies the quantities mentioned above (this facilitates error-checking) and data arranged in rows.  Here is an example file:
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
- Adjacency information.  This is contained in two methods: `getConnetionsFrom()` and `getConnectionsTo()`, which both return arrays.  Call the arrays `f` and `t`.  Then `f[i]` is the active cell index, and it is connected to the active cell with index `t[i]`.  So this is a sparse representation of a symmetric adjacency matrix.  Only active cells are considered, so other classes need not know about inactive cells.

`Grid` also contains various other utility methods, and an `outputActiveCSV(filename)` method.


