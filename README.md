# target_selection

![Versions](https://img.shields.io/badge/python->3.7-blue)
[![Documentation Status](https://readthedocs.org/projects/sdss-target-selection/badge/?version=latest)](https://sdss-target-selection.readthedocs.io/en/latest/?badge=latest)
<!-- [![Travis (.org)](https://img.shields.io/travis/sdss/target_selection)](https://travis-ci.org/sdss/target_selection)
[![codecov](https://codecov.io/gh/sdss/target_selection/branch/main/graph/badge.svg)](https://codecov.io/gh/sdss/target_selection) -->

Code to perform target selection for BHM/MWM using catalogdb.

## Installation

To install `target_selection` do

```console
$ pip install numpy
$ pip install sdss-target-selection
```

Due to a misconfiguration in `pymange`, `numpy` needs to be installed before you can install `target_selection`.

## Development

New code must follow the [SDSS Coding Standards](https://sdss-python-template.readthedocs.io/en/latest/standards.html). Linting checks run as a GitHub workflow after each commit. The workflow also check that the package and dependencies can be installed. Pull Requests that don't pass the checks cannot be merged.

To test the code offline, you must install [flake8](https://flake8.pycqa.org/en/latest/) and [isort](https://pycqa.github.io/isort/) and run the following commands from the root of the package

```console
$ flake8 . --count --show-source --statistics
$ isort -c python/ bin/
```

If the command exit without error, everything should be fine.
