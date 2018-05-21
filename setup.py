##packaging command python setup.py sdist

import sys
import os.path
from setuptools import setup, Extension, find_packages


setup(
    name='SplitFusion',
    version='0.0',
    packages=find_packages(),
    license='Creative Commons Attribution-Noncommercial-Share Alike license',
    long_description=open('README.md').read(),
    package_dir = {'SplitFusion': 'SplitFusion/'},
    package_data={'ReferenceGenomes': ['Reference_feat/*.txt']},
    py_modules = ['SplitFusion.Multiple_alignment', 'SplitFusion.Fusion_annotation', 'SplitFusion.Group_split'],
    install_requires=['subprocess32==3.2.7' ,
                        'HTSeq==0.9.1',
                        'statsmodels==0.8.0'],
    include_package_data=True

)


