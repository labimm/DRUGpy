#!/usr/bin/env python

import pathlib

from setuptools import find_packages, setup

HERE = pathlib.Path(__file__).parent
README = (HERE / "README.md").read_text()


setup(
    name="DRUGpy",
    version="1.0.0",
    description="Some PyMOL utilities",
    long_description=README,
    long_description_content_type="text/markdown",
    url="https://github.com/labimm/DRUGpy",
    author="Pedro Sousa Lacerda",
    author_email="pslacerda@gmail.com",
    license="MIT",
    classifiers=[
        "Development Status :: 4 - Beta",
        "License :: OSI Approved :: MIT License",
        "Environment :: Plugins",
        "Intended Audience :: Science/Research",
        "Natural Language :: English",
        "Programming Language :: Python :: 3",
        "Programming Language :: Python :: 3.7",
    ],
    packages=find_packages(),
    install_requires=[
        "lxml",
        "pandas",
        "scipy",
        "requests",
        "cached_property",
        "matplotlib",
        "seaborn",
        "jinja2",
    ],
)
