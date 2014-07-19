from setuptools import setup, find_packages

setup(
        # basic package data
        name = "Gimme tools",
        version = "0.97",

        # package structure
        packages=find_packages('src'),
        package_dir={'':'src'},

        install_requires = [
                            'networkx == 1.7',
                            'pygr == 0.8.2',
                            'bx-python == 0.7.1',
                            'bsddb3 == 6.0.1',
                            ]
        )
