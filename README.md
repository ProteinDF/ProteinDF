# ProteinDF

[![Build Status](https://www.travis-ci.org/ProteinDF/ProteinDF.svg?branch=master)](https://www.travis-ci.org/ProteinDF/ProteinDF)

## Manual

<https://proteindf.github.io/ProteinDF/>

## How to install ProteinDF

For build instructions refer to the INSTALL file in this directory. For further instructions on using ProteinDF,
refer to documents provided with ProteinDF.

ProteinDF works and has been tested on Linux system.
The standard configure scripts given will work for building it on any of the systems.

### Prerequisites

In order to build ProteinDF, you must have CMake, LAPACK, ScaLAPACK, MPI
and any related dev packages installed.

### build ProteinDF

Making building directory:

```bash
mkdir build
cd build
cmake ..
```

If you want to specify install path, you can run as following:

```bash
cmake -DCMAKE_INSTALL_PREFIX=[PATH] ..
```

### Building ProteinDF

To build ProteinDF, simply run:

```bash
make
```

#### OPTIONS for cmake

- CXXFLAGS

```bash
-DCMAKE_CXX_FLAGS="${CXX_FLAGS}"
```

The `-O3 -DNDEBUG` flag is usually recommended in order to improve performance.

- BLAS

```bash
-DBLAS_LIBRARIES="xxx"
```

- LAPACK

```bash
-DLAPACK_LIBRARIES="xxx"
```

- MPI

If you would like to specify MPI compiler, then

```bash
-DMPI_CXX_COMPILER="XXX"
```

otherwise,

```bash
-DMPI_CXX_COMPILER="xxx" \
-DMPI_CXX_INCLUDE_PATH="xxx" \
-DMPI_CXX_LIBRARIES="xxx"
```

- ScaLAPACK

```bash
-DSCALAPACK_LIBRARIES="xxx"
```

### Installing ProteinDF

To install ProteinDF, run:

```bash
make install
```

If you check the making progress, run:

```bash
make DO_VERBOSE=1
```

## License

Copyright (C) 2002-2024 ProteinDF developers.

This program is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation; either version 3 of the License, or (at your option) any later version.

This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more details.

You should have received a copy of the GNU General Public License along with this program.  If not, see <http://www.gnu.org/licenses/>

## Bugs

If you find any bugs, please let me know.
And if you have a suggestion for improvement, please let me know.
