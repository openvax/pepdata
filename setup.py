import os

readme_filename = os.path.join(os.path.dirname(__file__), 'README.md')
with open(readme_filename, 'r') as f:
  readme = f.read()

try:
  import pypandoc
  readme = pypandoc.convert(readme, to='rst', format='md')
except:
  #print "Conversion of long_description from markdown to reStructuredText failed, skipping..."
  pass


from setuptools import setup

if __name__ == '__main__':
    setup(
        name='epitopes',
        version="0.1",
        description="Python interface to IEDB and other immune epitope data",
        author="Alex Rubinsteyn",
        author_email="alex {dot} rubinsteyn {at} mssm {dot} edu",
        url="https://github.com/hammerlab/epitope",
        license="http://www.opensource.org/licenses/BSD-3-Clause",
        classifiers=[
            'Development Status :: 3 - Alpha',
            'Environment :: Console',
            'Operating System :: OS Independent',
            'Intended Audience :: Science/Research',
            'License :: OSI Approved :: Apache License',
            'Programming Language :: Python',
            'Topic :: Scientific/Engineering :: Bio-Informatics',
             ],
        requires=[
            'numpy(>=1.7)',
            'pandas(>=0.13.1)',
            'scikit.learn(>=0.14.1)',
            'appdirs',
        ],
        long_description=readme,
        packages=['epitopes'],
        package_data = { 'epitopes' : ['data/*csv'] },
        include_package_data = True
    )
