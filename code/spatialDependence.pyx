import time
import array
cimport cpython.array as array
cimport csvparser
from libc.stdlib cimport free

cdef class SpatialDependence:

    def __init__(self, float xmin, float ymin, float cell_size, unsigned nx, unsigned ny):

        self.xmin = xmin
        self.ymin = ymin
        self.cell_size = cell_size
        self.nx = nx
        self.ny = ny
        self.num_cells = nx * ny
        self.required_header = "#xmin=" + str(xmin) + ",ymin=" + str(ymin) + ",cell_size=" + str(cell_size) + ",nx=" + str(nx) + ",ny=" + str(ny)
        self.filetype = "NO_FILETYPE"
        self.uint_template = array.array('I', [])
        self.float_template = array.array('f', [])


    cpdef parse(self, str filename, str filetype, list required_additional_headers):
        if not (filetype == "active_inactive" or filetype == "wind_raw" or filetype == "wind_processed" or filetype == "generic_float" or filetype == "active_inactive_binary" or filetype == "wind_raw_binary" or filetype == "wind_processed_binary" or filetype == "generic_float_binary"):
            raise ValueError("filetype not recognized")
        self.filetype = filetype

        cdef int binary = (1 if filetype.endswith("binary") else 0)
        cdef char* header = NULL
        cdef size_t header_length = 0
        cdef size_t wind_length = 0
        cdef unsigned* udata = NULL
        cdef float* fdata = NULL
        cdef int error_code = 0
        cdef size_t expected_num_in_row = self.nx
        cdef size_t expected_num_data = self.num_cells
        if filetype.startswith("active_inactive"):
            error_code = csvparser.parseBool(filename.encode(), &header, &header_length, &udata, expected_num_in_row, expected_num_data, binary)
        elif filetype.startswith("wind_raw"):
            expected_num_in_row = 2 * self.nx
            expected_num_data = 2 * self.num_cells
            error_code = csvparser.parseFloat(filename.encode(), &header, &header_length, &fdata, expected_num_in_row, expected_num_data, binary)
        elif filetype.startswith("wind_processed"):
            expected_num_in_row = 3
            expected_num_data = 0
            error_code = csvparser.parseProcessedWind(filename.encode(), &header, &header_length, &udata, &fdata, &wind_length, expected_num_in_row, binary)
        elif filetype.startswith("generic_float"):
            error_code = csvparser.parseFloat(filename.encode(), &header, &header_length, &fdata, expected_num_in_row, expected_num_data, binary)

        if error_code == 1:
            raise IOError('Cannot open or read ' + filename)
        elif error_code == 2:
            raise MemoryError('C buffer not big enough to read ' + filename + '  You will have to increase MAX_FILE_LENGTH')
        elif error_code == 3:
            raise MemoryError('C memory allocation problem when reading ' + filename)

        # file was opened and header successfully read, so check header
        try:
            py_bytes_header = header[:header_length]
        finally:
            free(header)
        try:
            self.checkHeader(py_bytes_header.decode(), filename, required_additional_headers + [self.required_header])
        except:
            raise

        # header is OK, so check for poor data
        if error_code == 4:
            raise ValueError("There must be " + str(expected_num_data // expected_num_in_row) + " data lines in " + filename)
        elif error_code == 5:
            raise ValueError("There must be " + str(expected_num_in_row) + " entries per line in " + filename)
        elif error_code == 6:
            raise ValueError("The data entries in " + filename + " must be either 0 or 1")

        # shove the data into the data0, etc, arrays
        # i can't figure out how to just do self.data0 = uint, so have to copy
        cdef unsigned ind = 0
        if filetype.startswith("active_inactive"):
            self.data0 = array.clone(self.uint_template, self.num_cells, zero=False)
            for ind in range(self.num_cells):
                self.data0.data.as_uints[ind] = udata[ind]
            free(udata)
            return
        elif filetype.startswith("wind_raw"):
            self.data0 = array.clone(self.float_template, self.num_cells, zero=False)
            self.data1 = array.clone(self.float_template, self.num_cells, zero=False)
            for ind in range(self.num_cells):
                self.data0.data.as_floats[ind] = fdata[2 * ind]
                self.data1.data.as_floats[ind] = fdata[2 * ind + 1]
            free(fdata)
            return
        elif filetype.startswith("wind_processed"):
            self.data0 = array.clone(self.uint_template, wind_length, zero=False)
            self.data1 = array.clone(self.uint_template, wind_length, zero=False)
            self.data2 = array.clone(self.float_template, wind_length, zero=False)
            for ind in range(wind_length):
                self.data0.data.as_uints[ind] = udata[2 * ind]
                self.data1.data.as_uints[ind] = udata[2 * ind + 1]
                self.data2.data.as_floats[ind] = fdata[ind]
            free(udata)
            free(fdata)
            return
        elif filetype.startswith("generic_float"):
            self.data0 = array.clone(self.float_template, self.num_cells, zero=False)
            for ind in range(self.num_cells):
                self.data0.data.as_floats[ind] = fdata[ind]
            free(fdata)
            return

    cpdef parseWithPython(self, str filename, str filetype, list required_additional_headers):
        if not (filetype == "active_inactive" or filetype == "wind_raw" or filetype == "wind_processed" or filetype == "generic_float"):
            raise ValueError("filetype not recognized")
        self.filetype = filetype

        # open and read file
        try:
            with open(filename, 'r') as f:
                data = f.readlines()
        except:
            raise IOError('Cannot open or read ' + filename)

        # check the header
        try:
            self.checkHeader(data, filename, required_additional_headers + [self.required_header])
        except:
            raise

        cdef unsigned lendata
        # Size data with full-sized arrays
        if filetype == "active_inactive":
            self.data0 = array.clone(self.uint_template, self.num_cells, zero=False)
        elif filetype == "wind_raw":
            self.data0 = array.clone(self.float_template, self.num_cells, zero=False)
            self.data1 = array.clone(self.float_template, self.num_cells, zero=False)
        elif filetype == "wind_processed":
            lendata = len(data)
            self.data0 = array.clone(self.uint_template, lendata, zero=False)
            self.data1 = array.clone(self.uint_template, lendata, zero=False)
            self.data2 = array.clone(self.float_template, lendata, zero=False)
        elif filetype == "generic_float":
            self.data0 = array.clone(self.float_template, self.num_cells, zero=False)

        # Read the data
        cdef unsigned ind = 0
        cdef unsigned i, j
        cdef float fl
        for line in data:
            if not line.strip():
                continue
            line = line.strip()
            if line.startswith("#"):
                continue
            try:
                line = line.split(",")
                if filetype == "active_inactive":
                    if ind >= self.num_cells:
                        raise ValueError("There must be " + str(self.ny) + " data lines in " + filename)
                    if len(line) != self.nx:
                        raise ValueError("There must be " + str(self.nx) + " entries per line in " + filename)
                    for i in range(self.nx):
                        j = int(line[i])
                        if not (j == 0 or j == 1):
                            raise ValueError("The data entries in " + filename + " must be either 0 or 1")
                        self.data0.data.as_uints[ind] = j
                        ind += 1
                elif filetype == "wind_raw":
                    if ind >= self.num_cells:
                        raise ValueError("There must be " + str(self.ny) + " data lines in " + filename)
                    if len(line) != 2 * self.nx:
                        raise ValueError("There must be " + str(2 * self.nx) + " entries per line in " + filename)
                    for i in range(self.nx):
                        self.data0.data.as_floats[ind] = float(line[2 * i])
                        self.data1.data.as_floats[ind] = float(line[2 * i + 1])
                        ind += 1
                elif filetype == "wind_processed":
                    if len(line) != 3:
                        raise ValueError("There must be 3 entries per line in " + filename)
                    self.data0.data.as_uints[ind] = int(line[0])
                    self.data1.data.as_uints[ind] = int(line[1])
                    self.data2.data.as_floats[ind] = float(line[2])
                    ind += 1
                if filetype == "generic_float":
                    if ind >= self.num_cells:
                        raise ValueError("There must be " + str(self.ny) + " data lines in " + filename)
                    if len(line) != self.nx:
                        raise ValueError("There must be " + str(self.nx) + " entries per line in " + filename)
                    for i in range(self.nx):
                        fl = float(line[i])
                        self.data0.data.as_floats[ind] = fl
                        ind += 1
            except:
                raise 

        if filetype == "active_inactive" or filetype == "wind_raw" or filetype == "generic_float":
            if ind != self.num_cells:
                raise ValueError("There must be " + str(self.ny) + " data lines in " + filename)
        elif filetype == "wind_processed":
            array.resize(self.data0, ind)
            array.resize(self.data1, ind)
            array.resize(self.data2, ind)

    cdef checkHeader(self, str header_data, str filename, list required_headers):
        headers_found = []
        for line in header_data.split("\n"):
            if not line.strip():
                continue
            line = line.strip()
            # all filetypes must have this header:
            if line.startswith("#"):
                if line in required_headers:
                    headers_found.append(line)
                continue
            else:
                # must be starting to read numerical data
                break
        if len(headers_found) != len(required_headers):
            raise ValueError("Header lines in " + filename + " must include " + "\n".join(required_headers))


    cpdef array.array getData0(self):
        return self.data0
    cpdef array.array getData1(self):
        return self.data1
    cpdef array.array getData2(self):
        return self.data2

    cpdef void restrictToActive(self, array.array global_index):
        cdef array.array c0
        cdef array.array c1
        cdef unsigned num_active = len(global_index)
        cdef unsigned i
        if self.filetype.startswith("active_inactive"):
            c0 = array.clone(self.uint_template, num_active, zero = False)
            for i in range(num_active):
                c0.data.as_uints[i] = self.data0.data.as_uints[global_index.data.as_uints[i]]
            self.data0 = c0
        elif self.filetype.startswith("wind_raw"):
            c0 = array.clone(self.float_template, num_active, zero = False)
            c1 = array.clone(self.float_template, num_active, zero = False)
            for i in range(num_active):
                c0.data.as_floats[i] = self.data0.data.as_floats[global_index.data.as_uints[i]]
                c1.data.as_floats[i] = self.data1.data.as_floats[global_index.data.as_uints[i]]
            self.data0 = c0
            self.data1 = c1
        elif self.filetype.startswith("generic_float"):
            c0 = array.clone(self.float_template, num_active, zero = False)
            for i in range(num_active):
                c0.data.as_floats[i] = self.data0.data.as_floats[global_index.data.as_uints[i]]
            self.data0 = c0
        

    cpdef outputCSV(self, str filename, array.array active, str inactive_value, str additional_header_lines):
        f = open(filename, 'w')
        f.write("#File written at: " + time.asctime() + "\n")
        f.write(self.required_header + "\n")
        f.write(additional_header_lines)
        cdef unsigned x_ind
        cdef unsigned y_ind
        cdef unsigned ind, active_ind
        for y_ind in range(self.ny):
            for x_ind in range(self.nx - 1):
                ind = y_ind * self.nx + x_ind
                active_ind = active[ind]
                if active_ind == self.num_cells:
                    # inactive cell
                    f.write(inactive_value + ",")
                else:
                    f.write(str(self.data0[active_ind]) + ",")
            x_ind = self.nx - 1
            ind = y_ind * self.nx + x_ind
            active_ind = active[ind]
            if active_ind == self.num_cells:
                # inactive cell
                f.write(inactive_value + "\n")
            else:
                f.write(str(self.data0[active_ind]) + "\n")
        f.close()

