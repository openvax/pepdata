# Copyright (c) 2014. Mount Sinai School of Medicine
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#     http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.


from __future__ import print_function, division, absolute_import
import os

from setuptools import setup, find_packages

readme_dir = os.path.dirname(__file__)
readme_filename = os.path.join(readme_dir, 'README.md')

try:
    with open(readme_filename, 'r') as f:
        readme = f.read()
except:
    print("Failed to load README file")
    readme = ""

try:
    import pypandoc
    readme = pypandoc.convert(readme, to='rst', format='md')
except:
    print("Conversion of long_description from markdown to reStructuredText failed, skipping...")

if __name__ == '__main__':
    setup(
        name='pepdata',
        version="0.6.3",
        description="Python interface to IEDB and other immune epitope data",
        author="Alex Rubinsteyn",
        author_email="alex {dot} rubinsteyn {at} mssm {dot} edu",
        url="https://github.com/hammerlab/pepdata",
        license="http://www.apache.org/licenses/LICENSE-2.0.html",
        classifiers=[
            'Development Status :: 3 - Alpha',
            'Environment :: Console',
            'Operating System :: OS Independent',
            'Intended Audience :: Science/Research',
            'License :: OSI Approved :: Apache Software License',
            'Programming Language :: Python',
            'Topic :: Scientific/Engineering :: Bio-Informatics',
             ],
        install_requires=[
            'numpy>=1.7',
            'scipy>=0.9',
            'pandas>=0.13.1',
            'scikit-learn>=0.14.1',
            'progressbar33',
            'biopython>=1.65',
            'datacache>=0.4.4',
        ],
        long_description=readme,
        packages=find_packages(exclude="test"),
        package_data={'pepdata': ['data/*csv', 'matrices/*']},
        include_package_data=True
    )
