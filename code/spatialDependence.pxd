import time
import array
cimport cpython.array as array

from grid cimport Grid

cdef class SpatialDependence:
    """Parses files that define spatial dependence, building appropriate array or arrays.
    If required, the arrays can be restricted to the active cells only"""

    # xmin
    cdef float xmin

    # ymin
    cdef float ymin

    # cell side length
    cdef float cell_size

    # number of columns in the grid
    cdef unsigned nx

    # number of rows in the grid
    cdef unsigned ny

    # total number of cells
    cdef unsigned num_cells

    # required header in all files: #xmin=...
    cdef str required_header

    # filetype (active_inactive, wind_raw, wind_processed, etc)
    cdef str filetype

    # the imported data
    cdef array.array data0

    # more imported data
    cdef array.array data1

    # event more imported data
    cdef array.array data2

    cpdef parse(self, str filename, str filetype, list required_additional_headers)
    """Parses filename and places data into self.data internal data structures
    Filename must have the following format:
    - any empty lines are ignored
    - Data must be preceded by a line of the form '#xmin=a,ymin=b,cell_size=c,nx=d,ny=e' where a,b,c,d,e must match that defined by the Grid.  This allows for error checking
    - The file may also have required_additional_headers (a list of strings, which must all begin with '#')
    - Data is CSV form.
    The first row corresponds to the ymin set of cells
    The second row corresponds to the ymin+cell_size set of cells
    The third row corresponds to the ymin+2*cell_size set of cells, etc
    - Each row must typically have nx comma-separated entries (depending on filetype)
    - The entries must sometimes obey bounds (depending on filetype)"""

    cdef checkHeader(self, list data, str filename, list required_headers)
    """Checks the header.  If any errors, then a ValueError exception is raised.  Otherwise the header is fine and data reading may proceed"""

    cpdef array.array getData0(self)

    cpdef array.array getData1(self)

    cpdef array.array getData2(self)

    cpdef void restrictToActive(self, array.array global_index)
    """Restrict the data so it only occurs on active cells.
    Usually you should set global_index = grid.getGlobalIndex()
    Here global_index[i] = the global cell index of active cell i
    
    NOTE: YOU WILL HAVE TO getData AFTER CALLING THIS METHOD
    """

    cpdef outputCSV(self, str filename, array.array active, str inactive_value, str additional_header_lines)
    """Outputs data0 to a CSV file.
    Usually you should set active = grid.getActiveIndex()
    Usually the data will have been restricted to the active cells, so
    the value inactive_value (as a string) is used in the CSV file for inactive cells.
    A header line containing the current time will be written to the file
    A header line of the form #xmin... will be written to the file
    additional_header_lines will be written verbatim (including any \n, etc) into the file"""
