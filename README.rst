=====
PyCDT
=====

Python Charge Defects Toolkit (PyCDT) is a python package aimed at making 
charged defects modeling simpler, high throughput ready, and also accessible 
to researchers who don't have the required background. PyCDT can handle
thermodynamic calculations and error corrections in the context of periodic
boundary condition density functional calculations of charged defects in 
semiconductors and insulators. It can also generate the inputs for required 
DFT calculations and can process the output of the DFT calculations.
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
1. Navigate to pycdt root directory::

    cd pycdt

2. Install the code, using the command::

    python setup.py install

   The command tries to obtain the required packages and their dependencies and install them automatically.
   Access to root may be needed if ``virtualenv`` is not used.

3. The package can be installed at non-standard locations using the command::

    python setup.py install --prefix PYCDT_INSTALL_DIR

   where PYCDT_INSTALL_DIR is your choice of directory. In UNIX/Linux environments,
   add PYCDT_INSTALL_DIR to PATH and PYTHONPATH variables by the following commands::
    
    export PATH=$PATH:PYCDT_INSTALL_DIR
    export PYTHONPATH=$PYTHONPATH:PYCDT_INSTALL_DIR

4. (If not set) Set the VASP pseudopotential directory in $HOME/.pmgrc.yaml as follows::

     VASP_PSP_DIR: <Location of vasp pseudopotential top directory>

5. (If not set) Set the Materials Project API key in $HOME/.pmgrc.yaml as follows::

     MAPI_KEY: <Your mapi key obtained from www.materialsproject.org>



Examples
--------

From the pycdt root folder, go to examples folder by typing::

    cd examples

Questions?
---------
Post your questions on `PyCDT forum <https://groups.google.com/forum/#!forum/pycdt-forum>`_.

