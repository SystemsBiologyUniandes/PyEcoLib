import setuptools
import os
import glob

with open("README.md", "r") as fh:
    long_description = fh.read()

def add_datafiles(data_files, dest_dir, pattern):
    """Add directory structures to data files according to a pattern"""
    src_dir, elements = pattern
    def do_directory(root_dest_path, root_src_path, elements):
        files = []
        for e in elements:
            if isinstance(e, list):
                src_dir, elems = e
                dest_path = '/'.join([root_dest_path, src_dir])
                src_path = os.path.join(root_src_path, src_dir)
                do_directory(dest_path, src_path, elems)
            else:
                files.extend(glob.glob(os.path.join(root_src_path, e)))
        if files:
            data_files.append((root_dest_path, files))
    do_directory(dest_dir, src_dir, elements)

setuptools.setup(
    name="PyEcoLib",
    version="1.1.0",
    author="Cesar Vargas, César Nieto, Sergio Blanco",
    author_email="cavargar@gmail.com, canietoa@gmail.com, sergio.camilo.blanco@gmail.com",
    description="PyEcoLib (Python Library for E. coli size dynamics estimation) is library to estimate bacterial cell size stochastic dynamics including time-continuous growth process and division events.",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/SystemsBiologyUniandes/PyEcoLib",
    packages=setuptools.find_packages(),
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
    ],
    python_requires='>=3.6',
)

pyecolib_data_files = []
data_files = [('PyEcoLib', pyecolib_data_files)]
add_datafiles(data_files, 'PyEcoLib/PyEcoLib',
              ['*'])
add_datafiles(data_files, 'PyEcoLib/examples',
              ['examples'['*']])
