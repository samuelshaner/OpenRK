================================================================================
OpenCSG                         VERSION 0.1.1                           8/4/2014
================================================================================
         _______  _______  _______  _        _______  _______  _______ 
        (  ___  )(  ____ )(  ____ \( (    /|(  ____ \(  ____ \(  ____ \
        | (   ) || (    )|| (    \/|  \  ( || (    \/| (    \/| (    \/
        | |   | || (____)|| (__    |   \ | || |      | (_____ | |      
        | |   | ||  _____)|  __)   | (\ \) || |      (_____  )| | ____ 
        | |   | || (      | (      | | \   || |            ) || | \_  )
        | (___) || )      | (____/\| )  \  || (____/\/\____) || (___) |
        (_______)|/       (_______/|/    )_)(_______/\_______)(_______)

================================================================================


================================================================================
                         AUTHORS & CONTACT INFORMATION                         
================================================================================

Contributors--------------------------------------------------------------------

        Will Boyd         -- wboyd@mit.edu
        Davis Tran        -- dvtran@mit.edu
        Samuel Shaner     -- shaner@mit.edu
        Qicang Shen       -- shenqicangst@gmail.com

Faculty Advisors----------------------------------------------------------------

        Benoit Forget     -- bforget@mit.edu
        Kord Smith        -- kord@mit.edu


================================================================================
                                  INTRODUCTION
================================================================================

OpenCSG is a Python package for Constructive Solid Geometry (CSG) models. The
goal for OpenCSG is to provide an easy-to-use, physics-agnostic library to
build geometry models of nuclear reactor cores. Compatibility modules for 
commonly used nuclear reactor physics codes, such as OpenMC, OpenMOC, Serpent,
MCNP are being concurrently developed for rapid and easy exportation of an
OpenCSG model directly into the relevant input file format for each code of
interest. 

OpenCSG can provide users a number of key benefits, including the following:

  1) Easy validation & verification between codes of a single geometric model
  2) Pythonic interface facilitates easy parametric optimization studies
  3) CSG tree traversal routines to facilitate downstream data processing

The OpenCSG development project was started by graduate and undergraduate 
students in the MIT Nuclear Science & Engineering Department working in the
Computational Reactor Physics Group.


================================================================================
                             DOWNLOAD & INSTALLATION
================================================================================

Make sure Git an the NumPy, h5py, and matplotlib python modules are installed 
on your machine. These packages can easily be installed on any Linux or Mac
machine using a package manager such as apt-get, yum or macports. For machines
running Debian-based Linux, such as Ubuntu, the following terminal command will
install all of the necessary pre-requisites for OpenCSG:

    > sudo get-apt install git python-numpy python-h5py python-matplotlib

Next, download OpenCSG from GitHub using the following command:
    
    > git clone https://github.com/mit-crpg/OpenCSG.git

This will create an "OpenCSG" directory with the code. Enter this directory and
install the OpenCSG Python package with the following commands:

    > cd OpenCSG
    > python setup.py install --user

To test that everything works as expected, open up an interactive Python console
and import OpenCSG:

    > python
    > import opencsg

If this worked without any errors, you now have OpenCSG installed as the 
"opencsg" Python package importable from any directory on your machine. 


================================================================================
                                  USING OPENCSG
================================================================================

At the time of this writing, there is no documentation on the API provided to 
users by OpenCSG. It is expected that such documentation will be provided in an 
early release version of OpenCSG. Currently, the recommended way to learn the
OpenCSG API is through the example input files in the "OpenCSG/sample-input"
directory. These files show how a user may create simple pin cell and 
rectangular lattice geometries and plot them in OpenCSG.

To run the "pin-cell.py" sample input file, simply enter the directory and run
it in Python:

    > cd sample-input/pin-cell
    > python pin-cell.py

The compatibility modules for nuclear reactor physics codes, such as OpenMC, 
are distributed separately from OpenCSG with the codes themselves. To learn how
to use these modules to interface with OpenCSG, please reference the 
documentation provided with your code of interest.


================================================================================
                               ONGOING DEVELOPMENT
================================================================================

At the time of this writing, some of the key areas of ongoing development 
include the following:

  1) Ray tracing for arbitrary geometries
  2) Robust and fast volume calculations for arbitrary Cells
  3) Intrinsics for hexagonal lattices as commonly used in fast reactor designs
  4) Optimization using vectorization with NumPy, as well as Cython and/or Numba

If you would like to contribute to the OpenCSG project, please contact the 
development team.


================================================================================
                                    LICENSE                                    
================================================================================

OpenCSG is currently going through review from the MIT Technology Licensing
Office for open source distribution under the MIT/X License:

Copyright © 2014-2015 Massachusetts Institute of Technology

Permission is hereby granted, free of charge, to any person obtaining a copy 
of this software and associated documentation files (the “Software”), to 
deal in the Software without restriction, including without limitation the 
rights to use, copy, modify, merge, publish, distribute, sublicense, and/or 
sell copies of the Software, and to permit persons to whom the Software is 
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in 
all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED “AS IS”, WITHOUT WARRANTY OF ANY KIND, EXPRESS OR 
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, 
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE 
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER 
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING 
FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS 
IN THE SOFTWARE.
