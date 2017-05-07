ProteinDF
=========

[![Build Status](https://www.travis-ci.org/ProteinDF/ProteinDF.svg?branch=master)](https://www.travis-ci.org/ProteinDF/ProteinDF)

How to install ProteinDF
------------------------

For build instructions refer to the INSTALL file in this directory. For further instructions on using ProteinDF,
refer to documents provided with ProteinDF.


ProteinDF works and has been tested on Linux system.
The standard configure scripts given will work for building it 
on any of the systems.

### Prerequisites

In order to build ProteinDF, you must have CMake, LAPACK, ScaLAPACK, MPI 
and any related dev packages installed.


### build ProteinDF

Making building directory:

```
mkdir build
cd build
cmake ..
```

If you want to specify install path, you can run as following:

```
cmake -DCMAKE_INSTALL_PREFIX=[PATH] ..
```


### Building ProteinDF

To build ProteinDF, simply run:

```
make
```


### Installing ProteinDF

To install ProteinDF, run:

```
make install
```

License
-------

Copyright (C) 2002-2017 ProteinDF developers.

This program is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation; either version 3 of the License, or (at your option) any later version.

This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more details.

You should have received a copy of the GNU General Public License along with this program.  If not, see <http://www.gnu.org/licenses/>


Bugs
----

If you find any bugs, please let me know.
And if you have a suggestion for improvement, please let me know.

