from setuptools import setup, find_packages

with open("README.md", "r", encoding="utf-8") as fh:
    long_description = fh.read()

setup(
    name="geokitten",
    version="0.1.0",
    author="sape94",
    author_email="sergioapelayoe@gmail.com",
    description="Tools for geospatial data processing and visualization",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/sape94/GeoKitten",
    packages=find_packages(),
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
        "Topic :: Scientific/Engineering :: GIS",
    ],
    python_requires=">=3.9",
    install_requires=[
        "geopandas>=0.9.0",
        "numpy>=1.23.4",
        "geopandas>=0.14.2",
        "pandas>=2.2.2",
        "shapely>=2.0.2",
        "pytest>=8.3.5",
        "fiona>=1.8.22",
        "folium>=0.14.0",
        "branca>=0.6.0",
        "matplotlib>=3.7.1",
    ],
    extras_require={
        "dev": [
            "pytest>=8.3",
        ],
    },
    include_package_data=True,
)
