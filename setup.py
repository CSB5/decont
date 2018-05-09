# https://packaging.python.org/tutorials/distributing-packages/
import sys
from distutils.core import setup

python_requires='>=3',

PACKAGE_NAME="decont"
PACKAGE_VERSION="0.5"
PACKAGE_BUGREPORT="wilma@gis.a-star.edu.sg"

setup(name = PACKAGE_NAME,
      version = PACKAGE_VERSION,
      description="Decontaminate FastQ files by mapping with BWA-MEM against a given source",
      author="Andreas Wilm",
      author_email=PACKAGE_BUGREPORT,
      #requires = [], samtools and bwa-mem?
      scripts = [
          "decont.py",
      ],
      # http://pypi.python.org/pypi?%3Aaction=list_classifiers
      classifiers=['Environment :: Console',
                   'Intended Audience :: Science/Research',
                   'Natural Language :: English',
                   'Operating System :: Unix',
                   'Programming Language :: Python :: 3',
                   'Topic :: Scientific/Engineering :: Bio-Informatics',
                   ],
      keywords='bioinformatics'
      )
