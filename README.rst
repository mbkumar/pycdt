=====
pycdt
=====

Corrections to DFT calculation for point defects in solids

Requirements
------------
pycdt requires pymatgen and optionally sxdefectalign packages.

Source Code
------------
If not available already, use the following steps.

#. Install `git <http://git-scm.com>`_, if not already packaged with your system.

#. Download the pycdt source code using the command::

    git clone https://github.com/mbkumar/pycdt.git

Installation
------------
#. Navigate to pycdt root directory::

    cd pycdt

#. Install the code, using the command::

    python setup.py install

The command tries to obtain the required packages and their dependencies and install them automatically.
Access to root may be needed if ``virtualenv`` is not used.

# The package can be installed at non-standard locations using the command::

    python setup.py install --prefix PYCDT_ROOTDIR

where PYCDT_ROOTDIR is your choice of directory. In UNIX/Linux environments,
add PYCDT_ROOTDIR to PATH and PYTHONPATH variables by the following commands::
    
    export PATH=$PATH:PYCDT_ROOTDIR
    export PYTHONPATH=$PYTHONPATH:PYCDT_ROOTDIR


Examples
--------

From the pycdt root folder, go to examples folder by typing::

    cd examples

