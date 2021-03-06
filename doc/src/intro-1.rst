============
Introduction
============

About HIMG
----------

HIMG is an experimental image compression code based on 
adaptive higher-order finite element methods (hp-FEM). 
It is built upon the `Hermes2D library <http://hpfem.org/hermes2d>`_
and distributed under the `GNU General Public License 
<http://www.gnu.org/licenses/old-licenses/gpl-2.0.txt>`_. 

NOTICE
------

The Hermes library underwent substantial changes recently
(Hermes1D/2D/3D repos unified into a single repository,
hermes_common changed, assembling and solution procedures 
separated, several new matrix solver packages enabled, 
API changed). HIMG is now waiting to be updated accordingly. 
We will get to it as soon as we can.

Required Libraries
------------------

 * Hermes2D: Compatible with the version from May 25, 2010.
 * libpng: Compatible with a version 1.2.37.
 * UMFpack: This library is required by Hermes2D. See installation 
   instructions of Hermes2D `here <http://hpfem.org/hermes2d/doc/src/intro-2.html>`_.
 * cmake: This library is required by Hermes2D. See installation 
   instructions of Hermes2D.

Installation Instructions (Linux)
--------------------------------------

 1. Download, compile and install the library Hermes2D. You may need to  
    set specific variables in ``CMake.vars`` file.
 2. Install the libpng library.
 3. Switch to the folder containing HIMG and create the file ``CMake.vars``. 
    This file controls behavior of CMAKE. Supported variables are:

	- ``HERMES2D_ROOT``: A path to a folder containing an installed library Hermes2D. Use if the library Hermes2D is installed in a specific folder. The system assumes a sub-folder ``include``, which contains header files, and a sub-folder ``lib``, which contains library files. Ignore, otherwise.
	- ``HERMES2D_SOURCE_ROOT``: A path to a folder containing source code (including header files) of the library Hermes2D. This variable can be ignored if the variable ``HERMES2D_ROOT`` is used.
	- ``HERMES2D_LIBRARY_ROOT``: A path to a folder containing library files of the library Hermes2D. If not installed, the path is equal to ``HERMES2D_SOURCE_ROOT``, i.e., use ``set(HERMES2D_LIBRARY_ROOT ${HERMES2D_SOURCE_ROOT})``. This variable can be ignored if the variable ``HERMES2D_ROOT`` is used.
	- ``HERMES_COMMON_ROOT``: A path to a folder containing the library Hermes_common. Usually, it is a sub-folder of the library Hermes2D.
	
 4. Run CMAKE: ``cmake .``
 5. Compile using MAKE: ``make``
 
Running HIMG (Linux)
~~~~~~~~~~~~~~~~~~~~

The resulting application is controlled through command line parameters. For a complete list of parameters, run the application without any parameters. In order to try the sample image:

 1. Open a console if it is not already open.
 2. Switch to the folder containing HIMG.
 3. Run: ``./src/himg -i img/squares/squares.ppm``
	
Installation Instructions (MSVC)
--------------------------------

 1. Download and compile the library Hermes2D. Installation of the library Hermes2D will fail because CMAKE ignores UAC.
 2. Download and copy files to corresponding dependencies folders or any other folder that can be accessed by MSVC.
 3. Create and fill the file ``CMake.vars``. This file controls behavior of CMAKE. Supported variables are:

    - ``DEP_ROOT``: MSVC only. Set to a path of the dependencies folder (see installation instructions of Hermes2D).
	- ``HERMES2D_ROOT``: A path to a folder containing an installed library Hermes2D. Use if the library Hermes2D is installed in a specific folder. The system assumes a sub-folder ``include``, which contains header files, and a sub-folder ``lib``, which contains library files.
	- ``HERMES2D_SOURCE_ROOT``: A path to a folder containing a source code (including header files) of the library Hermes2D. Not necessary if the library was installed and it is accessible either at a standard path or at a path defined by the variable ``HERMES2D_ROOT``. Nevertheless, it is more comfortable to set this variable to a src/ folder of the library Hermes2D rather than installing the library.
	- ``HERMES2D_LIBRARY_ROOT``: A path to a folder containing library files (i.e. .lib) of the library Hermes2D. Set it to the folder in which the files .LIB are generated, e.g., use ``set(HERMES2D_LIBRARY_ROOT "my_hermes2d_root/src/Debug")``.
	- ``HERMES_COMMON_ROOT``: A path to a folder containing the library Hermes_common. Usually, it is a sub-folder of the library Hermes2D.
	- ``UMFPACK_ROOT``: Use ``set(UMFPACK_ROOT ${DEP_ROOT})``.
	- ``AMD_ROOT``: Use ``set(AMD_ROOT ${DEP_ROOT})``.
	- ``LIBPNG_ROOT``: Use ``set(LIBPNG_ROOT ${DEP_ROOT})``. If the library libpng was installed elsewhere, use an appropriate path. The script assumes that the path contains a sub-folder 'include' and a sub-folder 'lib'.
	
 4. Open a console, switch to a folder which contains HIMG and run ``cmake . -G "Visual Studio 9 2008"``.
 5. Open a created SLN file in MSVC and compile.
 
Running HIMG (MSVC)
~~~~~~~~~~~~~~~~~~~

The resulting application is controlled through command line parameters. For a complete list of parameters, run the application without any parameters. In order to try the sample image use the MSVC IDE. Set the command line parameters to ``-i ../img/squares/squares.ppm`` and run the application HIMG.


