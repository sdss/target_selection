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

A workflow checks for linting and formatting errors on each commit, and pull requests are blocked until the workflow succeeds. The easiest way to fix these problems is by installing `ruff` and letting it format the code, and then checking if any linting errors remain. The commands `ruff format` and `ruff check` are independent and both must be run.

```console
$ pip install ruff
$ ruff format ./python/
1 file reformatted, 46 files left unchanged
$ ruff check ./python/
All checks passed!
```

Updating the package version can be done directly in the `pyproject.toml` file and doesn't require having `poetry` installed or otherwise updating the lockfile.

### Visual Studio Code configuration

If using Visual Studio Code, it is recommended to install the [ruff](https://marketplace.visualstudio.com/items?itemName=charliermarsh.ruff) and [prettier](https://marketplace.visualstudio.com/items?itemName=esbenp.prettier-vscode) extensions. Then you can create a workspace file inside the cloned repo, under `.vscode/settings.json` with the following configuration

```json
{
    "[python]": {
        "editor.formatOnSave": true,
        "editor.codeActionsOnSave": {
            "source.fixAll": "explicit",
            "source.organizeImports.ruff": "explicit"
        },
        "editor.wordWrap": "off",
        "editor.tabSize": 4,
        "editor.defaultFormatter": "charliermarsh.ruff"
    },
    "[markdown]": {
        "editor.wordWrapColumn": 88
    },
    "[restructuredtext]": {
        "editor.wordWrapColumn": 88
    },
    "[json]": {
        "editor.quickSuggestions": {
            "strings": true
        },
        "editor.suggest.insertMode": "replace",
        "editor.formatOnSave": true,
        "editor.defaultFormatter": "esbenp.prettier-vscode",
        "editor.tabSize": 2
    },
    "[yaml]": {
        "editor.insertSpaces": true,
        "editor.formatOnSave": true,
        "editor.defaultFormatter": "esbenp.prettier-vscode",
        "editor.tabSize": 2,
        "editor.autoIndent": "advanced",
    },
    "prettier.tabWidth": 2,
    "editor.rulers": [99],
    "editor.wordWrapColumn": 99,
    "python.analysis.typeCheckingMode": "off",
    "ruff.nativeServer": true
}
```

which will apply the formatting and linting automatically on save.

### Working with `sdssdb`

Often developing `target_selection` requires concurrent changes to [sdssdb](https://github.com/sdss/sdssdb). This presents a bit of a challenge to keep everything in sync and provide a convenient developing environment.

The recommended way to develop with `target_selection` and `sdssdb` is to install `sdssdb` in editable mode in the `target_selection` environment. To do this run

```bash
pip install -e <PATH-TO-SDSSDB-ROOT-DIR>
```

After this if you run

```bash
pip list | grep sdssdb
```

you should see a path next to the version, indicating that `sdssdb` is being imported from that path. Any local changes in that path will be applied when `target_selection` imports `sdssdb`.

### Updating `sdssdb` and other dependencies before tagging

When tagging `target_selection` please make sure that you've tagged `sdssdb` (if needed) and that the lockfile reflects the change. First you'll need to install [poetry](https://python-poetry.org/docs/#installation) (this should only be required once). Follow the installation instructions there, but a recommended way to install is:

1. Install [pipx](https://pipx.pypa.io/stable/installation/). Follow the instructions for your system, but generally `pip install pipx` works fine.
2. Use `pipx` to install `poetry` with

    ```bash
    pipx install poetry
    ```

3. Ensure that the `poetry` executable can be found. If you execute `poetry` in your shell and that files, make sure that you have `$HOME/.local/bin` in your `$PATH` or follow other `pipx` instructions.

Once `poetry` is installed follow these instructions:

1. If there are changes in `sdssdb` that affect `target_selection`, tag `sdssdb`. Please do update the `sdssdb` [change log](https://github.com/sdss/sdssdb/blob/main/CHANGELOG.rst) adding the header indicating the date in which the release was made. Note that you don't necessarily need to do this every time that `target_selection` is tagged, only when the new `target_selection` tag depends on untagged changes in `sdssdb`.

2. Update the `version` in `pyproject.toml` to the release version.

3. Update the dependency for `sdssdb` in `pyproject.toml` for example if `sdssdb` 0.14.1 has just been tagged, change

    ```toml
    sdssdb = "0.14.0"
    ```

    to

    ```toml
    sdssdb = "0.14.1"
    ```

4. Recreate the lockfile with

    ```bash
    poetry lock
    ```

    It may take a bit of time for `poetry` to detect a recently released version of `sdss`. You may want to run `poetry lock --no-cache` which will take a bit longer but should grab the new package as long as PyPI has updated its internal cache.

    Note that `poetry lock` *will not* update the version of `sdssdb` in your virtual environment, only the lockfile. If you do want to update the dependency and the lockfile use `poetry update sdssdb` (or `poetry update` to update all the dependencies).

5. Commit the new changes (including the `poetry.lock` file) and tag the new version of `target_selection`.
