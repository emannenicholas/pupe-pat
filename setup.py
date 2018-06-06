"""
pupepat - Pupil Plate Ellipticity (PUPE)- Photometric Analysis Tool (PAT)

Author
    Curtis McCully (cmccully@lco.global)
    Brook Taylor
    Rob Siverd

License
    GPL v3.0

February 2018
"""
from setuptools import setup

setup(name='pupepat',
      author=['Curtis McCully', 'Brook Taylor'],
      author_email=['cmccully@lco.global', ''],
      version='0.2.3',
      packages=['pupepat'],
      package_dir={'pupepat': 'pupepat'},
      setup_requires=['pytest-runner'],
      install_requires=['numpy', 'astropy', 'lcogt_logging', 'sep', 'matplotlib', 'watchdog', 'scipy', 'corner',
                        'emcee', 'pypdf2'],
      tests_require=['pytest'],
      entry_points={'console_scripts': ['run_pupepat_listener=pupepat.main:run_watcher',
                                        'run_pupepat=pupepat.main:analyze_directory']})
