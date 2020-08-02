#Copyright 2013 Argonne UChicago LLC
#This file is part of MachBall

import setuptools

with open("README.md", "r") as fh:
    long_description = fh.read()

setuptools.setup(
    name="machball",
    version="0.2.0",
    author="Angel Yanguas-Gil",
    author_email="ayg@anl.gov",
    description="self-limited (ALD, ALE) reactive transport in nanostructures",
    long_description=long_description,
    url="https://github.com/aldsim/machball",
    packages=setuptools.find_packages(),
    install_requires=["numpy"],
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: BSD License",
        "Operating System :: OS Independent",
        "Intended Audience :: Science/Research",
        "Intended Audience :: Manufacturing"
        ],
    python_requires='>=3.6')
