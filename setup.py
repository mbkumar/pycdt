import os
import glob

from setuptools import setup, find_packages

SETUP_PTH = os.path.dirname(os.path.abspath(__file__))

def readme():
    with open('README.rst') as f:
        return f.read()

setup(
        name="pycdcd",
        packages=find_packages(),
        version="0.0.1",
        install_requires=["pymatgen>=3.0.6", "numpy>=1.8"],
        author="Bharat Medasani",
        author_email="mbkumar@gmail.com",
        url="http://github.com/mbkumar/pycdcd",
        description="PyCDCD is a python tool to obtain corrections to the "
                  "DFT calculations for point defects in solids",
        long_description=readme(),
        classifiers=[
            "Programming Language :: Python :: 2",
            "Programming Language :: Python :: 2.7",
            "Programming Language :: Python :: 3.4",
            "Development Status :: 1 - Alpha",
            "Intended Audience :: Science/Research",
            "License :: OSI Approved :: MIT License",
            "Operating System :: OS Independent",
            "Topic :: Scientific/Engineering :: Information Analysis",
            "Topic :: Scientific/Engineering :: Physics",
            "Topic :: Scientific/Engineering :: Chemistry",
            "Topic :: Software Development :: Libraries :: Python Modules"
            ],
        license="MIT",
        scripts=glob.glob(os.path.join(SETUP_PTH, "scripts", "*")),
        test_suite='nose.collector',
        tests_require=['nose'],
        #entry_points={
        #    'console_scripts':[
        #        'gen_def_structure=pydii.scripts.gen_def_structure:'
        #                       'im_vac_antisite_def_struct_gen',
        #        'gen_def_energy=pydii.scripts.gen_def_energy:'
        #                       'im_vac_antisite_def_energy_parse',
        #        'gen_def_profile=pydii.scripts.gen_def_profile:'
        #                       'im_vac_antisite_def_profile',
        #        'gen_sol_pref_structure=pydii.scripts.gen_def_structure:'
        #                       'im_sol_sub_def_struct_gen',
        #        'gen_sol_def_energy=pydii.scripts.gen_def_energy:'
        #                       'im_sol_sub_def_energy_parse',
        #        'gen_sol_site_pref=pydii.scripts.gen_def_profile:'
        #                       'im_sol_sub_def_profile']}
)

