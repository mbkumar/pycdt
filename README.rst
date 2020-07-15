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

------------
Requirements
------------
PyCDT requires pymatgen (and its dependencies) and optionally sxdefectalign packages.

------------
Source Code
------------
If not available already, use the following steps.

#. Install `git <http://git-scm.com>`_, if not already packaged with your system.

#. Download the pycdt source code using the command::

    git clone https://bitbucket.org/mbkumar/pycdt.git

------------
Installation
------------
0. (Optional step, but strongly suggested)
   Download and install anaconda and create a new virtual environment. 
   For example, to create a new anaconda based virtual environment with a name *pycdt_venv*, run the following steps.::

    conda create -n pycdt_venv
    conda install -n pycdt_venv matplotlib numpy scipy setuptools
   Activate this virtual environment, by running the command::

    conda activate pycdt_venv
   If everything goes well, you should see *pycdt_venv* at the command line prompt.
   Don't forget to activate your virtual environment whenever you are trying to 
   update or use pycdt.

1. Navigate to pycdt root directory::

    cd pycdt

2. Install the code, using the command::

    python setup.py install

   The command tries to obtain the required packages and their dependencies and install them automatically.
   (Access to root may be needed if ``virtualenv`` is not used.) You can run the ``pycdt`` now. 

3. (Optional for power users) The package can be installed at non-standard locations using the command::

    python setup.py install --prefix PYCDT_INSTALL_DIR

   where PYCDT_INSTALL_DIR is your choice of directory. In UNIX/Linux environments,
   add PYCDT_INSTALL_DIR to PATH and PYTHONPATH variables by the following commands::
    
    export PATH=$PATH:PYCDT_INSTALL_DIR
    export PYTHONPATH=$PYTHONPATH:PYCDT_INSTALL_DIR

4. Configure your VASP pseudopotentials for use with pymatgen.::

    pmg config -p <EXTRACTED_VASP_POTCAR> <MY_PSP>

   Here  *<EXTRACTED_VASP_POTCAR>* is the folder where your pseudopotentials are present and 
   *<MY_PSP>* is where the directory where the layout of pseudopotentials is organized  by pymatgen.
   For more information refer to `pymatgen installation instructions <https://pymatgen.org/installation.html>`_.

5. Set the VASP pseudopotential directory in $HOME/.pmgrc.yaml using the command::

    pmg config --add PMG_VASP_PSP_DIR <MY_PSP>

   where <MY_PSP> is the directory used in the previous step.


6. Set the Materials Project API key in $HOME/.pmgrc.yaml using the command::

     pmg config -add PMG_MAPI_KEY <Your mapi key obtained from www.materialsproject.org>


Documentation
------------
For help on how to use PyCDT, read the `article <https://doi.org/10.1016/j.cpc.2018.01.004>`_ published in Computer Physics Communication.

Examples
--------
From the pycdt root folder, go to examples folder by typing::

    cd examples


Questions?
----------
Post your questions on `PyCDT forum <https://groups.google.com/forum/#!forum/pycdt-forum>`_.

Citation
--------
To support PyCDT development, if you use pycdt in your research, please cite the following Comp Phys Comm article.

- Danny Broberg, Bharat Medasani, Nils E.R. Zimmermann, Guodong Yu, Andrew Canning, Maciej Haranczyk, Mark Asta, Geoffroy Hautier,
  PyCDT: A Python toolkit for modeling point defects in semiconductors and insulators,
  Computer Physics Communications, Volume 226, 2018, Pages 165-179.

