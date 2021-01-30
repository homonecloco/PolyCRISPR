import setuptools
from distutils.core import setup, Extension


with open("README.md", "r", encoding="utf-8") as fh:
    long_description = fh.read()

setuptools.setup(
    name="PolyCRISPR", # Replace with your own username
    version="0.0.1",
    author="Ricardo H. Ramirez-Gonzalez",
    author_email="ricardo.ramirez-gonzalez@jic.ac.uk",
    description="A package to count amplicons and detect the origin",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/homonecloco/PolyCRISPR",
    packages=setuptools.find_packages(),
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
    ],
    python_requires='>=3.6',
    test_suite='setup.unit_tests',
    #ext_modules=[module1]
    entry_points={  # Optional
        'console_scripts': [
            'PolyCRISPR=PolyCRISPR:main'
        ],
    }
)
