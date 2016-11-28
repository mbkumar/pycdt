=====
PyCDT
=====

Python Charge Defects Toolkit (PyCDT) is a python packaged aimed at making 
charged defects modeling simpler, high throughput ready, and also accessible 
to researchers who don't have the required background. At present it can do 
simple defect thermodynamic calculations and error corrections for periodic
boundary condition density functional calculations of charged defects in 
semiconductors and insulators. It can also generate the inputs for required 
DFT calculations and can processes the output of the DFT calculations.
The code is modular and any DFT code can be integrated into PyCDT for defect 
calculations. 

Requirements
------------
PyCDT requires pymatgen (and its dependencies) and optionally sxdefectalign packages.

Source Code
------------
If not available already, use the following steps.

#. Install `git <http://git-scm.com>`_, if not already packaged with your system.

#. Download the pycdt source code using the command::

    git clone https://bitbucket.org/mbkumar/pycdt.git

Installation
------------
#. Navigate to pycdt root directory::

    cd pycdt

#. Install the code, using the command::

    python setup.py install

    The command tries to obtain the required packages and their dependencies and install them automatically.
    Access to root may be needed if ``virtualenv`` is not used.

#. The package can be installed at non-standard locations using the command::

    python setup.py install --prefix PYCDT_INSTALL_DIR

    where PYCDT_INSTALL_DIR is your choice of directory. In UNIX/Linux environments,
    add PYCDT_INSTALL_DIR to PATH and PYTHONPATH variables by the following commands::
    
    export PATH=$PATH:PYCDT_INSTALL_DIR
    export PYTHONPATH=$PYTHONPATH:PYCDT_INSTALL_DIR

#. (Optional) Set the VASP pseudopotential directory in .pmgrc.yaml as follows::

     VASP_PSP_DIR: <location of vasp pseudopotential top directory>



Examples
--------

From the pycdt root folder, go to examples folder by typing::

    cd examples

