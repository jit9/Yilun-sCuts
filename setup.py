from setuptools import setup

setup(
    name='cutslib',
    version='0.1.0',
    packages=['cutslib'],
    include_package_data=True,
    scripts=['bin/cutspipe',
             'bin/run_cuts',
             'bin/submit_cuts',
             'bin/submit_cuts_tiger',
             'bin/update_crit',
             'bin/update_cuts',
             'bin/map_cuts.py',
             'bin/map_planet.py',
             'bin/plot_cuts',],
    install_requires=['dotmap',
                      'configparser',
                      'pycook',
                      'seaborn',
                      'ipdb',
                      'jinja2',
                      'future', # from here onwards will be moby2 deps
                      'numpy',
                      'matplotlib',
                      'scipy',
                      'astropy',
                      'ephem',
                      'pytz',
                      'h5py',
                      'profilehooks',
                      'fitsio',],
    entry_points={'console_scripts': [
        'cuts=cutslib.cli:main'
    ]},
)
