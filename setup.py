# -*- coding: utf-8 -*-
from setuptools import setup, find_packages

try:
    long_description = open("README.rst").read()
except IOError:
    long_description = ""

setup(
    name="gro_parser",
    version="0.1.0",
    description="Python package to parse and work with gro files",
    license="MIT",
    author="Cecile Hilpert",
    packages=find_packages(),
    install_requires=[],
    long_description=long_description,
    classifiers=[
        "Programming Language :: Python",
        "Programming Language :: Python :: 3.10",
    ]
)
