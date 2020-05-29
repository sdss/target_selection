
.. _database-server:

The database server
===================

For SDSS-V target selection we use a Postgresql database server running on a dedicated machine, ``operations.sdss.org``, at University of Utah. Details on the server and its configuration are given below.

The parent catalogues and the results of :ref:`cross-matching <cross-matching>` and :ref:`target selection <target-selection>` live in the ``sdss5db`` database, whose schema can be found :ref:`here <sdssdb:available-databases>`. The files used to create and populate the parent catalogues are part of the `sdssdb <https://github.com/sdss/sdssdb>`_ product. Each parent catalogue must be a well documented publication; the downloaded files are in general maintained at Utah at ``/uufs/chpc.utah.edu/common/home/sdss50/sdsswork/target/catalogs``. Additional information about each catalogue can be found in the `wiki <https://wiki.sdss.org/x/Y4DzAQ>`__.

Database preparation
--------------------

The ``operations`` server has 32 CPU cores and approximately 250 GB of RAM with 80 TB of conventional hard drive space. The disks are configured as RAID-6 (?) with about 3TB of SSD disk space used as cache. The OS and software runs on a small SSD hard drive partition.

We run a Postgresql 12 server in which we create the ``sdss5db`` database. We define two main roles in the server, ``sdss`` with full privileges on the database, and ``sdss_user`` with read only permissions to ``catalogdb`` and ``targetdb``. We create an additional schema ``sandbox``, writeable by ``sdss_user`` that can be used for testing by any user.

Each parent catalogue is loaded manually. SQL files and scripts to load each table exists in the directory ``schema/sdss5db/catalogdb`` of sdssdb_. Overall, the loading process is as follows:

- Create the new table without data, primary key, or indexes to speed up data insertion.
- Load the data from CSV files using `COPY <https://www.postgresql.org/docs/12/sql-copy.html>`__. For catalogues provided as FITS tables, sdssdb includes a series of tools for efficient ingestion. For large catalogues with multiple input files, the process can be parallelised but since this is an IO-heavy process, it does not make sense to have many more than 6 or 7 concurrent processes (the number of independent disks in the RAID).
- Create indexes and add the primary key. Add foreign keys, remembering to index the referring column.

For each table with coordinates we create a `Q3C <https://github.com/segasai/q3c>`__ index on them for efficient spatial querying, and cluster the table on that index (note that this may decrease the efficiency of other indexes). Alternatively it's sometimes faster to recreate the table sorted by the index, ``CREATE TABLE AS SELECT * FROM tab ORDER BY q3c_ang2ipix(ra,dec);``, with the added advantage that no locks are created. See `this issue <https://github.com/segasai/q3c/issues/24#issuecomment-610716846>`__ for details.

We also add the `pg_healpix <https://github.com/segasai/pg_healpix>`__ and `pgunit <https://github.com/petere/pguint>`__ extensions. Although it's not installed at Utah, `pg_repack <https://github.com/reorg/pg_repack>`__ provides a way of keeping tables clustered on an index without locks.

Server configuration
--------------------

The ``postgresql.conf`` file with the configuration for the database server is kept under version control in the ``operations`` table of the `config product <https://github.com/sdss/config>`__. The default Postgresql configuration is very conservative and not well suited for a database with table with billions of rows. In general is a good idea to start by using a website such as `PGTune <https://pgtune.leopard.in.ua/>`__ to determine the best default configuration depending on the server parameters. Here are the final configurations values:

.. code-block:: text

    search_path = '"$user",public,catalogdb,targetdb'
    max_connections = 20
    shared_buffers = 64GB
    effective_cache_size = 3000GB
    maintenance_work_mem = 2GB
    checkpoint_completion_target = 0.9
    wal_buffers = 16MB
    effective_io_concurrency = 9
    work_mem = 100MB
    min_wal_size = 4GB
    max_wal_size = 16GB
    max_worker_processes = 32
    max_parallel_workers_per_gather = 16
    max_parallel_workers = 32
    max_parallel_maintenance_workers = 16

    random_page_cost = 0.2
    seq_page_cost = 0.1
    cpu_index_tuple_cost 0.0001
    cpu_operator_cost 0.0025
    default_statistics_target = 500
    temp_buffers = 500MB

    autovacuum_max_workers = 3
    autovacuum_vacuum_threshold = 50
    autovacuum_vacuum_scale_factor = 0.001
    autovacuum_analyze_scale_factor = 0.002

    shared_preload_libraries = 'auto_explain'
    auto_explain.log_analyze = true
    auto_explain.log_min_duration = '100s'
    auto_explain.log_buffers = 'on'
    auto_explain.log_format = 'text'

    log_destination = 'stderr'
    logging_collector = 'on'
    log_statement = 'all'
    log_directory = '/scratch/pg_logs'
    log_filename = 'postgresql-%Y-%m-%d_%H%M%S.log'
    log_file_mode = '0666'

All the parameters are explained in the Postgresql documentation but we offer comments on a few key ones with suggestions for the values to use:

- *shared_buffers* is the amount of memory used by all the Postgresql server processes. Postgresql loads data from tables and indexes here to operate on them. A reasonable value is between 25% and one third of the total RAM. Postgresql will eventually use all this memory and won't return it. A larger value doesn't seem to improve efficiency. A good explanation of the different types of memory used by Postgresql can be found `here <https://severalnines.com/database-blog/architecture-and-tuning-memory-postgresql-databases>`__

- *work_mem* is the memory used by individual processes for hash operations such as sorts or joins. Note that this is the maximum memory allowed to *each* of such operations so if a query has three hash joins and a sort it will use four times the amount of work_mem. Because of that it's better to set a conservative value (but larger than the default) and change it locally inside specific transactions (using ``SET LOCAL work_mem = 'X'``). This can help efficiency substantially but must be used with care. In some queries involving large tables setting it to 10GB or even larger helps but one must keep an eye on the memory usage lest the server runs out of RAM and crashes.

- *effective_cache_size* is only used by the query planner to determine the approximate size of the disk cache and optimise queries. For a system without fast disk caching this should be the rest of the memory not used by shared_buffers (between two thirds and 75%). In our case we set it to about 3TB. The value in itself is not important and you won't run out of memory because of it, but too low a value will lean the planner towards sequential scans and a very large value will use more indexes (which is only good if fast access to to the index is possible).

- *temp_buffers* is the memory allowed for temporary tables. Any temporary table that requires more than this amount of memory is written to disk. As with work_mem it's best to have a conservative default value and modify it locally within a transaction.

- *maintenance_work_mem* is the memory use by *each* of the maintenance processes such as ``VACUUM`` or ``CREATE INDEX``. 2GB is a reasonable value that can be increased locally within a transaction. In general it's recommended to increase this value by a lot during the initial database loading.

- *effective_io_concurrency* indicates how many concurrent disk I/O operations are allowed. This is a complicated value to set in a system with a RAID and SSD cache but in general it seems that setting it to the number of disks in the RAID is reasonable value.

- *random_page_cost* and *seq_page_cost* indicate the relative cost of performing a sequential read of a table versus a random (index) access. In spinning disks random accesses are up to four times more costly than sequential ones; in SSD disks they are almost equivalent. This values are used by the query planner to calculate the cost associated to operations and determine whether to use sequential scans or indexes. We decrease their relative cost while lowering their absolute value with respect to CPU operations. This results in the planner using indexes for most cases except for the scanning of very large tables for which most or all rows need to be returned. More details are given `here <https://www.postgresql.org/docs/current/runtime-config-query.html>`__.

- *cpu_index_tuple_cost* and *cpu_operator_cost* are the costs of processing each index entry during an index scan, and each operator or function, respectively. They don't seem to impact the query planner very heavily but we reduce them to about a tenth of their original value to account for faster, modern CPUs.

- *default_statistics_target* is the fraction of the table that is read during ``ANALYZE`` to create statistics about table and index sizes. The default value is 100 and we increase it to 500 which seems to be a good compromise between reasonable fast runs of ``ANALYZE`` and accurate statistics.

- The autovacuum parameters are changed to make sure that up to three autovacuum workers are spun when tables are modified. ``AUTOVACUUM`` does not get triggered until a certain fraction of the table has changed. The default values usually fail to trigger a vacuum in large tables so we increase their sensitivity by decreasing the value of ``autovacuum_vacuum_scale_factor`` and ``autovacuum_analyze_scale_factor``.

- We configure ``auto_explain`` to log to file the ``EXPLAIN`` of each query that takes more than 100 seconds. This, along with tools such as `PEV2 <https://dalibo.github.io/pev2/#/>`__ are very useful to determine why slow queries are so and what the query planner is doing. Here is a `blog post <https://www.depesz.com/2013/04/16/explaining-the-unexplainable/>`__ explaining how to read an ``EXPLAIN ANALYZE``.

- Finally we enable file logging to ``/scratch/pg_logs``.

- For production we do not modify the *fsync* or *synchronous_commit* parameters since we don't see a very significant improvement and they entail some risk. During the initial database loading it's probably a good idea to at least set ``synchronous_commit=off``. More details are available `here <https://www.postgresql.org/docs/12/runtime-config-wal.html#RUNTIME-CONFIG-WAL-ARCHIVING>`__.

Connecting and using the database
---------------------------------

These instructions assume that you have access to the Utah system. The operations machine is only accessible by a small group of users but the server is available from any other virtual machine. Once you have ssh'd to your favourite Utah machine create or edit your `~/.pgpass <https://www.postgresql.org/docs/12/libpq-pgpass.html>`__ with the line ::

    *:*:sdss5db:sdss:XXX
    *:*:sdss5db:sdss_user:XXX

where ``XXX`` is the password that needs to be requested from an administrator. Set the permissions to the file by running ``chmod 0600 ~/.pgpass``. Once you've done that you should be able to connect to the server as ::

    psql -U sdss_user -h operations.sdss.org sdss5db

Alternatively you can create an ssh tunnel to any machine at Utah and forward the Postgresql port in operations, e.g. ::

    ssh -l {your_utah_username} -L {a_port_of_your_choosing}:operations.sdss.org:5432 manga.sdss.org

(you can use any machine, not only the manga VM). Then add a ``.pgpass`` file as the above in your system and do ::

    psql -U sdss_user -h localhost -p {a_port_of_your_choosing} sdss5db

Apart from test, raw SQL queries, the easiest way to work with the database is by using sdssdb_. From a machine at Utah it's best to use the ``operations`` miniconda Python install by doing ``module load miniconda/3.7.7_operations`` which includes the latest version of sdssdb. Once that is loaded you can import the database connection and models for catalogdb or targetdb by doing ::

    >>> from sdssdb.peewee.sdss5db import database
    >>> database.set_profile('operations')
    True
    >>> database.connected
    True
    >>> from sdssdb.peewee.sdss5db import catalogdb, targetdb

Refer to the sdssdb_ documentation for details on how to use the profiles, models, and other connection options.

If you are connecting via an ssh tunnel the ``operations`` profile won't work, instead do ::

    >>> from sdssdb.peewee.sdss5db import database
    >>> database.connect_from_parameters(user='sdss_user', host='localhost', port={a_port_of_your_choosing})
    True

There are some more details and tips on using the database server in the `wiki <https://wiki.sdss.org/x/oIBsAw>`__.

Tips for running queries efficiently
------------------------------------

(These tips are written in raw SQL but they are equally applicable if you're using sdssdb/ORM). While testing queries, especially long-running ones, it's important to make sure a limit is applied in some way. The easiest way is to add a ``LIMIT`` to the query to return only the first N results (make sure to order your query if you want the results to be reproducible). For example:

.. code-block:: postgresql

    SELECT * FROM catalog c
        INNER JOIN catalog_to_tic_v8 ctic USING (catalogid)
        INNER JOIN tic_v8 tic ON tic.id = ctic.target_id
        INNER JOIN gaia_dr2_source gaia ON gaia.source_id = tic.gaia_int
    WHERE gaia.parallax < 0.5
    ORDER BY gaia.parallax DESC
    LIMIT 100;

will return the catalog information for the 100 Gaia targets with the largest parallaxes as long as those are < 0.5. This query runs in ~20s while the query without the ``LIMIT`` could take more than one hour. You can also use aggregate functions such as ``COUNT(*)`` to get statistics from your queries.

Alternatively, it's possible to limit your query by doing a radial query:

.. code-block:: postgresql

    SELECT * FROM catalog c
        INNER JOIN catalog_to_tic_v8 ctic USING (catalogid)
        INNER JOIN tic_v8 tic ON tic.id = ctic.target_id
        INNER JOIN gaia_dr2_source gaia ON gaia.source_id = tic.gaia_int
    WHERE q3c_radial_query(c.ra, c.dec, 100, 20, 1);

This query will return all the catalog rows that are cross-matched with Gaia DR2 and that fall within a radius of 1 degree around (100, 20) deg.

For very large queries it's best to avoid using a naked SELECT statement that output to the terminal. For example ``SELECT * FROM unwise`` will return 2 billion rows and 300 columns. What's more, the output will probably be larger than the RAM size available and you'll crash the database server. And even if the query works you won't be able to process it in any useful way from the screen. Instead, save the results to a new table:

.. code-block:: postgresql

    CREATE TABLE sandbox.temp_results AS SELECT * FROM unwise;

All users can use the ``sandbox`` schema for this purpose (it's writeable even by the ``sdss_user`` role). Remember to drop your table once you're done with it. You can use temporary tables but note those will disappear automatically once you close the connection.
