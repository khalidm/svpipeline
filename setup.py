#!/usr/bin/env python

from setuptools import setup

setup(
    name='svpipeline',
    version='0.1',
    author='Khalid Mahmood',
    author_email='khalid.mahmood@unimelb.edu.au',
    packages=['src'],
    entry_points={
        'console_scripts': ['svpipeline = src.main:main']
    },
    url='https://github.com/khalidm/svpipeline',
    license='LICENSE.txt',
    description='svpipeline is a bioinformatics pipeline to call structural \
     variants.',
    long_description=open('README.md').read(),
    install_requires=[
        "ruffus == 2.6.3",
        "drmaa == 0.7.6",
        "PyYAML == 3.11"
    ],
)
