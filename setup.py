#!/usr/bin/env python
from setuptools import setup

setup(name='yt_illustris_halos',
      packages=['yt_illustris_halos'],
      version='0.1.0',
      description='Download halos from the Illustris simulations into a format suitable for yt.',
      author='John ZuHone',
      author_email='jzuhone@gmail.com',
      url='http://github.com/jzuhone/yt_illustris_halos',
      download_url='https://github.com/jzuhone/yt_illustris_halos/tarball/0.1.0',
      install_requires=["six"],
      classifiers=[
          'Intended Audience :: Science/Research',
          'Operating System :: OS Independent',
          'Programming Language :: Python :: 2.7',
          'Programming Language :: Python :: 3.4',
          'Programming Language :: Python :: 3.5',
          'Topic :: Scientific/Engineering :: Visualization',
      ],
      )