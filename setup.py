#!/usr/bin/env python

import setuptools

with open("README.md", "r") as f:
    long_description = f.read()

setuptools.setup(name='solsticepy'
	,version='0.3.7'
	,author='Ye Wang'
	,author_email='ye.wang@anu.edu.au'
	,description="Scripts and wrappers for Solstice ray-tracing software"
	,long_description=long_description
	,long_description_content_type="text/markdown"
	,url='https://github.com/anustg/solstice-scripts'
	,packages=setuptools.find_packages()
	,license="GPL v3.0 or later, see LICENSE file"
	,classifiers=[
		"Development Status :: 4 - Beta"
		,"Environment :: Console"
		,"Intended Audience :: Science/Research"
		,"License :: OSI Approved :: GNU General Public License v3 or later (GPLv3+)"
		,"Natural Language :: English"
		,"Operating System :: Microsoft :: Windows :: Windows 10"
		,"Operating System :: POSIX :: Linux"
		,"Programming Language :: Python :: 2"
		,"Programming Language :: Python :: 3"
		,"Topic :: Scientific/Engineering :: Physics"
	]
	,install_requires=['scipy','numpy','uncertainties','matplotlib','colorama']
	,python_requires='>=2.7'
)

