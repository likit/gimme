from setuptools import setup, find_packages

setup(
        name = 'Gimme',
        version = '0.8',
        author = 'Likit Preeyanon',
        author_email = 'preeyano@msu.edu',
        license = 'GPL',

        packages=find_packages('src'),
        package_dir={'':'src'},

        entry_points = {
            'console_scripts': [
                        'gimme = gimme.gimme'
                        ]
                    },

        install_requires = [
                            'networkx',
                            ],
        )
