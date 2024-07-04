# target_selection

![Versions](https://img.shields.io/badge/python->=3.10-blue)
[![Documentation Status](https://readthedocs.org/projects/sdss-target-selection/badge/?version=latest)](https://sdss-target-selection.readthedocs.io/en/latest/?badge=latest)

Code to perform target selection for BHM/MWM using catalogdb.

## Installation

To install `target_selection` do

```console
pip install sdss-target-selection
```

or from the GitHub repository

```console
git clone git@github.com:sdss/target_selection
cd target_selection
pip install .
```

## Development

This code adheres to the [SDSS Coding Standards](https://sdss-python-template.readthedocs.io/en/latest/standards.html).

We use [poetry](https://python-poetry.org) as the PEP517 backend and for dependency specification and resolution. To install `poetry` follow [these instructions](https://python-poetry.org/docs/#installation). The you can install the project for development with

```console
cd target_selection
poetry install
```

Note that as long as you don't need to install or update dependencies you can still install `target_selection` in editable mode with

```console
pip install -e .
```

(this will *not* install the development packages under the `dev` dependencies group, those need to be manually pip-installed in this case). Please, **do not add new dependencies without updating the `poetry.lock` file**. One of the workflows that run on commit will check that the lockfile is still valid.

We use [ruff](https://docs.astral.sh/ruff/) for both linting, import sorting, and code formatting. The configuration is stored in `pyproject.toml` and it's mainly the default `ruff` configuration (and similar to `flake8`) but with a line length of 99 characters for historical reasons and because it simplifies writing long Peewee/SQLAlchemy query statements. The formatting is similar to [black](https://github.com/psf/black), and thus quite opinionated.

A workflow checks for linting and formatting errors on each commit, and pull requests are blocked until the workflow succeeds. The easiest way to fix these problems is by installing `ruff` and letting it format the code, and then checking if any linting errors remain

```console
$ pip install ruff
$ ruff format ./python/
1 file reformatted, 46 files left unchanged
$ ruff check ./python/
All checks passed!
```

Updating the package version can be done directly in the `pyproject.toml` file and doesn't require having `poetry` installed or otherwise updating the lockfile.
