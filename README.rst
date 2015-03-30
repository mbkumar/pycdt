=====
pycdcd
=====

Corrections to DFT calculation for point defects in solids

Installation
------------
pycdcd requires pymatgen package.

Source Code
------------
If not available already, use the following steps.

#. Install `git <http://git-scm.com>`_, if not already packaged with your system.

#. Download the pydii source code using the command::

    git clone https://github.com/mbkumar/pycdcd.git

Installation
------------
#. Navigate to pycdcd root directory::

    cd pycdcd

#. Install the code, using the command::

    python setup.py install

The command tries to obtain the required packages and their dependencies and install them automatically.
Access to root may be needed if ``virtualenv`` is not used.

# The package can be installed at non-standard locations using the command::

    python setup.py install --prefix PYCDCD_ROOTDIR

where PYCDCD_ROOTDIR is your choice of directory. In UNIX/Linux environments,
add PYCDCD_ROOTDIR to PATH and PYTHONPATH variables by the following commands::
    
    export PATH=$PATH:PYCDCD_ROOTDIR
    export PYTHONPATH=$PYTHONPATH:PYCDCD_ROOTDIR


Examples
--------

From the pycdcd, go to examples folder by typing::

    cd examples

