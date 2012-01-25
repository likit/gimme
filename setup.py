from setuptools import setup, find_packages

setup(
        name = 'Gimme',
        version = '0.7.4',
        author = 'Likit Preeyanon',
        author_email = 'preeyano@msu.edu',
        license = 'GPL',

        packages=find_packages('src'),
        package_dir={'':'scr'},

        install_requires = [
                            'networkx',
                            ],
        )
