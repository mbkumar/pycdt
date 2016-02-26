=====
PyCDT
=====

Python Charge Defects Tool is quite useful for processing the DFT calculations of charged point defects in semiconductors and insulators. It can generate the inputs for required DFT calculations and can processes the output of the DFT calculations to correct the errors arising due to the limitations of supercell methodology in evaluating the long range Coulomb interactions of defect charges. 

Requirements
------------
PyCDT requires pymatgen and optionally sxdefectalign packages.

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

    python setup.py install --prefix PYCDT_INSTALL_DIR

where PYCDT_INSTALL_DIR is your choice of directory. In UNIX/Linux environments,
add PYCDT_INSTALL_DIR to PATH and PYTHONPATH variables by the following commands::
    
    export PATH=$PATH:PYCDT_INSTALL_DIR
    export PYTHONPATH=$PYTHONPATH:PYCDT_INSTALL_DIR


Examples
--------

From the pycdt root folder, go to examples folder by typing::

    cd examples

