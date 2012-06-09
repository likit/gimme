from setuptools import setup, find_packages

setup(
        # basic package data
        name = "Gimme tools",
        version = "0.9",

        # package structure
        packages=find_packages('src'),
        package_dir={'':'src'},
        )
