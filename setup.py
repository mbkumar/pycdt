import os
import glob

from setuptools import setup, find_packages

SETUP_PTH = os.path.dirname(os.path.abspath(__file__))

def readme():
    with open('README.rst') as f:
        return f.read()

setup(
        name="pycdt",
        packages=find_packages(),
        version="2.0.4",
        install_requires=["numpy>=1.18.1", "pymatgen>=2020.4.29", "matplotlib>=3.1""monty>=3.0.2"],
        package_data={"pycdt.utils": ["*.yaml"]},
        author="Danny Broberg, Bharat Medasani, Nils Zimmerman",
        author_email="mbkumar@gmail.com",
        maintainer="Bharat Medasani",
        maintainer_email="mbkumar@gmail.com",
        url="http://bitbucket.org/mbkumar/pycdt",
        description="PyCDT is a python package to facilitate "
                  "DFT calculations for point defects in solids",
        long_description=readme(),
        classifiers=[
            "Programming Language :: Python :: 3",
            "Programming Language :: Python :: 3.8",
            "Development Status :: 2 - Release",
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
        test_suite="nose.collector",
        tests_require=["nose"]
)

