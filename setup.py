from setuptools import setup, find_packages

if __name__ == "__main__":

    with open("README.md", "r") as readme:
        long_desc = readme.read()

    setup(
        name="kaplan",
        version="1.3a.0",
        description="Conformer searching package.",
        long_description=long_desc,
        long_description_content_type="text/markdown",
        url="https://github.com/PeaWagon/Kaplan",
        author="Jennifer H. Garner",
        author_email="garnej2@mcmaster.ca",
        package_dir={"kaplan": "kaplan"},
        # user will need to install these manually since
        # psi4, vetee, and GOpt are not available on pypi
        # install_requires=[
        #    "numpy",
        #    "scipy",
        #    "openbabel",
        #    "psi4",
        #    "GOpt",
        #    "rmsd",
        #    "vetee",
        #    "matplotlib",
        #    "pubchempy",
        # ],
        python_requires=">=3.7",
        packages=find_packages(),
        classifiers=[
            "Programming Language :: Python :: 3.7",
            "License :: OSI Approved :: MIT License",
            "Operating System :: OS Independent",
            "Development Status :: 3 - Alpha",
            "Intended Audience :: Science/Research",
            "Topic :: Scientific/Engineering :: Chemistry",
        ],
        keywords="conformer optimisation geometry chemistry \
                  evolutionary-algorithm extinction \
                  ring-optimisation",

    )
