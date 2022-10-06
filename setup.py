#!/usr/bin/env python

import setuptools
from numpy.distutils.core import Extension, setup

setup(
    name="ChemVapDep",
    version="1.0",
    author="Daniele Padula",
    author_email="dpadula85@yahoo.it",
    description='''
        A python package to run Chemical Vapout Deposition MD simulations with
        GROMACS.
        ''',
    url="https://github.com/dpadula85/ChemVapDep",
    packages=setuptools.find_packages(),
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: GNU GPL License",
        "Operating System :: OS Independent",
    ],
    entry_points={ 
        'console_scripts' : [
            'deposit=ChemVapDep.deposit:main',
            ]
        },
    zip_safe=False
)
