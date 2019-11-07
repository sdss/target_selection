# -*- coding: utf-8 -*-
#
# Licensed under a 3-clause BSD license.
#
# @Author: Brian Cherinka
# @Date: 2017-09-27
# @Filename: tasks.py
# @License: BSD 3-Clause
# @Copyright: Brian Cherinka


from __future__ import absolute_import, division, print_function, unicode_literals

import os

from invoke import Collection, task


# This file contains tasks that can be easily run from the shell terminal using the Invoke
# python package. If you do not have invoke, install it with pip install
# To list the tasks available, type invoke --list from the top-level repo directory

@task
def clean_docs(ctx):
    """Cleans up the docs"""
    print('Cleaning the docs')
    ctx.run('rm -rf docs/sphinx/_build')


@task
def build_docs(ctx, clean=False):
    """Builds the Sphinx docs"""

    if clean:
        print('Cleaning the docs')
        ctx.run('rm -rf docs/sphinx/_build')

    print('Building the docs')
    os.chdir('docs/sphinx')
    ctx.run('make html', pty=True)


@task
def show_docs(ctx):
    """Shows the Sphinx docs"""
    print('Showing the docs')
    os.chdir('docs/sphinx/_build/html')
    ctx.run('open ./index.html')


@task
def clean(ctx):
    """Cleans up the crap before a Pip build"""

    print('Cleaning')
    ctx.run('rm -rf htmlcov')
    ctx.run('rm -rf build')
    ctx.run('rm -rf dist')
    ctx.run('rm -rf *.egg-info')
    ctx.run('rm -rf python/*.egg-info')


@task(clean, default=True)
def deploy(ctx):
    """Deploy the project to PyPI"""
    print('Deploying to PyPI!')
    ctx.run('python setup.py sdist bdist_wheel --universal')
    ctx.run('twine upload dist/*')


@task(clean)
def deploy_test(ctx):
    """Deploy the project to the test version of  PyPI"""
    print('Deploying to Test PyPI!')
    ctx.run('python setup.py sdist bdist_wheel --universal')
    ctx.run('twine upload --repository-url https://test.pypi.org/legacy/ dist/*')


@task(name='install-deps')
def install_deps(ctx, extras=None):
    """Install only dependencies from setup.cfg."""

    import setuptools

    if not os.path.exists('setup.cfg'):
        raise RuntimeError('setup.cfg cannot be found. If your project uses '
                           'requirement files use pip install -r instead.')

    if extras:
        extras = extras.split(',')
    else:
        extras = []

    config = setuptools.config.read_configuration('setup.cfg')

    if not config['options']:
        return

    options = config['options']

    setup_requires = options.get('setup_requires', [])
    install_requires = options.get('install_requires', [])

    requires = setup_requires + install_requires
    requires_str = (' '.join('"' + item + '"' for item in requires))
    if len(requires) > 0:
        ctx.run(f'pip install --upgrade {requires_str}', pty=True)

    for extra in extras:
        print(f'Installing extras={extra}')
        if 'extras_require' not in options:
            raise RuntimeError('extras_require is not defined')
        extra_deps = options['extras_require'].get(extra, [])
        if len(extra_deps) > 0:
            extra_deps_str = (' '.join('"' + item + '"' for item in extra_deps))
            ctx.run(f'pip install --upgrade {extra_deps_str}', pty=True)


os.chdir(os.path.dirname(__file__))

# create a collection of tasks
ns = Collection(clean, install_deps)

# create a sub-collection for the doc tasks
docs = Collection('docs')
docs.add_task(build_docs, 'build')
docs.add_task(clean_docs, 'clean')
docs.add_task(show_docs, 'show')
ns.add_collection(docs)

deploy_task = Collection('deploy')
deploy_task.add_task(deploy, 'pypi')
deploy_task.add_task(deploy_test, 'test')
ns.add_collection(deploy_task)
