# mozzie
Mosquito lifecycle, diffusion and advection

## Spatial structure

Mosquitoes are assumed to advect and diffuse over a grid of square cells, defined by:
- (xmin, ymin): the lower left-hand corner
- the cell side-length
- (nx, ny): the number of cells in the x and y directions
The cells can be "active" or "inactive".  This is specified through a CSV file.

## Directory layout

- tests directory contains tests (.py files) and associated files (all other files).  Run the tests by using, for example, `python TestGrid.py`
- code directory contains the core code that numerically simulates mosquito population dynamics
- code/auxillary directory contains python scripts that perform auxillary functions, such as plotting results.  These scripts depend on lots of python libraries, and, while useful, are not necessary for the numerical simulation of mosquitoes.

