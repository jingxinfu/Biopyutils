#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# Author            : Jingxin Fu <jingxin_fu@outlook.com>
# Date              : 09/02/2020
# Last Modified Date: 09/02/2020
# Last Modified By  : Jingxin Fu <jingxin_fu@outlook.com>
import setuptools
from Biopyutils import __version__

with open("README.md", "r") as fh:
    long_description = fh.read()
NAME='Biopyutils'
setuptools.setup(
    name=NAME,
    version=__version__,
    author="Jingxin Fu",
    author_email="jingxinfu.tj@gmail.com",
    description="A bionformatic tookits",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://jingxinfu.github.io/"+NAME,
    packages=setuptools.find_packages(),
    scripts=['bin/'+NAME],
    package_data={NAME: ["data/*"],},
    include_package_data=True,
    install_requires=['pandas','numpy'],
    python_requires='>=2.7, <4',
    keywords= ['Gene ID Transfer', 'Bioinformatics','Genomics','Computational Biologist'],
    classifiers=[
        "Programming Language :: Python",
        "License :: OSI Approved :: GNU General Public License v3 or later (GPLv3+)",
        "Operating System :: Unix",
        "Operating System :: MacOS",
    ]
)
