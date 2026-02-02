#!/usr/bin/env python3
"""
Setup script for RhoNNO MLX integration package.
"""

from setuptools import setup, find_packages

with open("README.md", "r", encoding="utf-8") as fh:
    long_description = fh.read()

setup(
    name="rhonno_mlx",
    version="1.0.0",
    author="RhoNNO Team",
    author_email="",
    description="MLX-based neural network implementations for RhoNNO",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/rhonno/rhonno",
    packages=find_packages(),
    classifiers=[
        "Development Status :: 3 - Alpha",
        "Intended Audience :: Science/Research",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
        "Programming Language :: Python :: 3",
        "Programming Language :: Python :: 3.8",
        "Programming Language :: Python :: 3.9",
        "Programming Language :: Python :: 3.10",
        "Programming Language :: Python :: 3.11",
        "Topic :: Scientific/Engineering :: Artificial Intelligence",
    ],
    python_requires=">=3.8",
    install_requires=[
        "mlx>=0.5.0",
        "numpy>=1.21.0",
        "torch>=1.12.0",  # For compatibility/conversion
    ],
    extras_require={
        "dev": ["pytest>=6.0", "black", "flake8"],
    },
)