import setuptools

with open("README.md", "r") as fh:
    long_description = fh.read()

setuptools.setup(
    name="PyEcoLib",
    version="1.0.1",
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