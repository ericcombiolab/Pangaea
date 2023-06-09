"""
@Author: Yu Huang
@Email: yuhuang-cst@foxmail.com
"""

from setuptools import setup, find_packages, Extension
import numpy as np
import sys

if sys.version_info.major != 3:
    raise RuntimeError('RPH-KMeans requires Python 3')

# Set this to True to enable building extensions using Cython.
# Set it to False to build extensions from the C file (that
# was previously created using Cython).
USE_CYTHON = False

install_requires = open('requirements.txt').read().splitlines()
cmdclass = {}
include_dirs = ['rph_kmeans', np.get_include()]

if USE_CYTHON:
    from Cython.Distutils import build_ext
    install_requires.append('cython >= 0.26.1')
    cmdclass.update({'build_ext':build_ext})
    sources = ['rph_kmeans/_point_reducer_cy.pyx', 'rph_kmeans/_point_reducer_cy_lib.cpp']
    cmdclass['build_ext'] = build_ext
else:
    sources = ['rph_kmeans/_point_reducer_cy.by_cython.cpp', 'rph_kmeans/_point_reducer_cy_lib.cpp']

setup(
    name='rph_kmeans',
    description='KMeans algorithm with initial centers produced by random projection point reduction process.',
    version='1.0.0',
    author='Yu Huang',
    author_email='yuhuang-cst@foxmail.com',
    packages=find_packages(),
    zip_safe=False,
    url='https://github.com/tinglabs/rph_kmeans',
    license='LICENSE',
    long_description=open('README.md').read(),
    install_requires=install_requires,
    cmdclass = cmdclass,
    ext_modules=[
        Extension(
            'rph_kmeans/_point_reducer_cy',
            sources=sources,
            include_dirs=include_dirs,
            depends=["rph_kmeans/_point_reducer_cy_lib.h"],
            language='c++',
            extra_compile_args=['-std=c++11', '-O2', '-Wall'],
            extra_link_args=['-std=c++11'],
        ),
    ],
)

