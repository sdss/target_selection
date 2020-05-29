.. highlight:: sh

.. _cli:

Command-line interface
======================

The command ``target_selection`` provides access to some of the features of the library.

Examples of use
---------------

When ``target_selection`` is imported a database connection to ``sdss5db`` is attempted using the profile that sdssdb identifies as matching the domain. To use a specific profile ::

   target_selection --profile operations

or one can provide the database connection parameters. Most commands require a read-write connection to the database so we set the user to ``sdss`` ::

   target_selection --user sdss --host operations.sdss.org

We can also enable the verbose mode by passing ``-v``. To run a cross-match run we use the ``xmatch`` subcommand with the cross-match plan version in the configuration file ::

   target_selection --user sdss xmatch "0.1.2"

The ``run`` subcommand will execute target selection for a given plan and load the results into ``targetdb`` ::

   target_selection --verbose --user sdss run --overwrite "0.4.1"

This will execute target selection plan ``0.4.1`` from the configuration file, overwriting the intermediate table if it exists.

We can remove the results of a target selection run from ``targetdb`` ::

   target_selection --user sdss clear "0.4.1"

Note that this won't remove entries from the ``target`` or ``magnitude`` table because a target entry can be related to multiple target selection plans.

Available commands
------------------

This is an auto-generated list the subcommands and parameters available for the ``target_selection`` CLI. It is equivalent of running ``target_selection <subcommand> --help`` for each subcommand.

.. click:: target_selection.__main__:target_selection
   :prog: target_selection
   :show-nested:
