from setuptools import setup, find_packages

with open("README.md", "r", encoding="utf-8") as f:
    long_description = f.read()

setup(
    name="dimer_search",  # Your package name
    version="0.1.0",  # Initial version
    description="A package for creating, manipulating and computing molecular dimers.",
    long_description=long_description,
    long_description_content_type="text/markdown",
    author="Alejandro Santana-Bonilla",  # Replace with your name
    author_email="",  # Replace with your email
    url="https://github.com/alejandrosantanabonilla/dimer_search",  # Replace with your repo URL
    packages=find_packages(),  # Automatically find packages
    install_requires=[  # List your dependencies
        "numpy<2",
        "rdkit-pypi",
        "openbabel-wheel",
        "scipy",
        "parmed",
        "ase",
        "tblite",
        "wheel",
    ],
    classifiers=[  # Add classifiers for your package
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: GNU general License",  # Choose your license
        "Operating System :: OS Independent",
    ],
    python_requires=">=3.7",  # Specify minimum Python version
)
