#!/usr/bin/env python3
from setuptools import setup, find_packages

setup(
    name="goprobioticprofiler",
    version="1.0.0",
    description="GO enrichment and GO-Slim profiler for probiotic bacterial genomes",
    author="PROBEAT Team",
    packages=find_packages(),
    python_requires=">=3.11",
    install_requires=[
        "pandas>=2.0",
        "goatools>=1.4.0",
        "numpy>=1.24",
        "scipy>=1.10",
        "statsmodels>=0.14",
        "pyyaml>=6.0",
    ],
    entry_points={
        "console_scripts": [
            "goprobioticprofiler=goprobioticprofiler.cli:main",
        ],
    },
    classifiers=[
        "Programming Language :: Python :: 3.11",
        "License :: OSI Approved :: MIT License",
        "Topic :: Scientific/Engineering :: Bio-Informatics",
    ],
)
