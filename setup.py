from setuptools import setup

setup(
    name='cutslib',
    version='0.1.0',
    packages=['cutslib'],
    scripts=['bin/promote_version',
             'bin/cutspipe',
             'bin/run_cuts',
             'bin/submit_cuts'],
    install_requires=['dotmap',
                      'configparser',
                      'pycook']
)
