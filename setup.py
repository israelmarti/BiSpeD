from setuptools import setup, find_packages

with open("README.md", "r", encoding="utf-8") as fh:
    long_description = fh.read()

setup(
    name="bisped",
    version="1.7",
    author="Cristian I. Martínez",
    description="Binary Spectral Disentangling (BiSpeD) for astronomical spectroscopy",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/israelmarti/BiSpeD",
    include_package_data=True,
    packages=find_packages(),
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
    ],
    python_requires=">=3.8",
    install_requires=[
        "astropy>=7.1.0",
        "matplotlib>=3.10.6",
        "numba>=0.62.1",
        "numpy>=2.3.2",
        "progress>=1.6.1",
        "PyAstronomy>=0.24.0",
        "scipy>=1.16.1",
        "specutils>=2.2.0",
        "psutil>=5.9.0",
    ],
)
