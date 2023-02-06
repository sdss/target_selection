
import os
import time
from itertools import combinations

import numpy as np
import peewee
import yaml
from peewee import fn
from playhouse.postgres_ext import ArrayField

from sdssdb.peewee.sdss5db.catalogdb import database
from sdssdb.peewee.sdss5db.targetdb import Target

import target_selection


database.connect(dbname='sdss5db')  # Use 'sdss' to get read/write permissions


class TempMatch(peewee.Model):
    """Model for the first output table with all the matches."""
    lead = peewee.TextField(null=False)
    lowest_catalogid = peewee.BigIntegerField(index=True, null=False)
    highest_catalogid = peewee.BigIntegerField(index=True, null=False)

    class Meta:
        database = database
        schema = 'sandbox'


class UniqueMatch(peewee.Model):
    """Model for the final table with unique pairs of catalogids"""
    lowest_catalogid = peewee.BigIntegerField(index=True, null=False)
    highest_catalogid = peewee.BigIntegerField(index=True, null=False)
    match_tables = ArrayField(peewee.TextField, null=False, index=False)
    n_tables = peewee.SmallIntegerField(null=False)

    class Meta:
        database = database
        schema = 'sandbox'


class MetaXMatch:
    """Correlates two crossmatch runs.

    This class correlates to catalogdb crossmatch runs by taking all the catalogids from the
    older run restricting to the ones used in target selection and giving for each the
    corresponding catalogid(s) from the second run, independent if the catalogid from the
    second run has been included in target selection or not. To do this it searches all
    the catalogids linked to the targetid associated to the target catalogid from the
    first run in the relational table.

    To instantiate the object we need the name of the configarion file and the database
    with which we will work.


    The configuration file contains the xmatch plans to be matched, the name of the log file,
    split_insert_number to indicate how many entries are ingested simultaneously to the first
    output table, the name of another configuration file indicating details about
    each crossmatch run, and 4 optional parameters. There parameters are sample_region to test
    the code in a single region, database options, show first to display the first N results
    in the log for each table, and split_query to indicate the tables in which the main query
    is split into sub queries restricted to target_id ranges. The individual crossmatch run
    config file on the other hand indicates the version_id and tables of each xmatch run,
    mainly to calculate the intersecting tables between the 2 runs, to look for the matches
    in those tables.

    Parameters
    ----------
    database : peewee:PostgresqlDatabase
        A `peewee:PostgresqlDatabase` to the database the tables to
        cross-match.
    config_filename : str
        Configuration file containing at least parameters xmatch_plans, individual_xmatch_config
        log_file, and split_inrest_number.
    """

    def __init__(self, config_filename, database):
        self.database = database
        config = yaml.load(open(config_filename, 'r'), Loader=yaml.SafeLoader)
        self.first_plan = config['xmatch_plans'][0]
        self.second_plan = config['xmatch_plans'][1]
        log_filename = config['log_file']
        self.log = target_selection.log
        log_path_and_name = os.path.realpath('./' + log_filename)
        self.log.start_file_logger(log_path_and_name, rotating=False, mode='a')
        self.log.sh.setLevel(0)

        self.log.info(f'Matching Plans: {self.first_plan} '
                      f'with {self.second_plan}')
        ind_xmatch_info_filename = config['individual_xmatch_config']
        ind_xmatch_config = yaml.load(open(ind_xmatch_info_filename, 'r'),
                                      Loader=yaml.SafeLoader)
        first_plan_config = ind_xmatch_config[self.first_plan]
        second_plan_config = ind_xmatch_config[self.second_plan]
        self.first_version = first_plan_config['version_id']
        self.second_version = second_plan_config['version_id']
        self.log.info(f'Corresponding to versions {self.first_version}, and {self.second_version}')

        optional_parameters = ['sample_region', 'database_options', 'show_first', 'split_query']

        # All the optional parameters present in the configuration file are stored directly
        # As attributes of the MetaXMatch object and for the ones that are not the attribute
        # get a None value.

        for param in optional_parameters:
            if param in config.keys():
                setattr(self, param, config[param])
            else:
                setattr(self, param, None)

        # If sample_region is included in the configuration file, then the word full is ommited
        # from the output table and instead we include the ra and dec centers and the radius used

        if self.sample_region:
            racen, deccen, radius = self.sample_region
            st_dec = str(deccen)
            if '-' not in st_dec:
                st_dec = 'pos' + st_dec
            else:
                st_dec = st_dec.replace('-', 'neg')
                output_name = (f'catalog_ver{self.first_version}_to_ver{self.second_version}_all'
                               f'_ra{racen}_dec{st_dec}_rad{int(radius)}')
            log_message = f'###  Using ra={racen:5.1f}  dec={deccen:5.1f}  rad={radius:04.1f}  ###'
        else:
            output_name = f'catalog_ver{self.first_version}_to_ver{self.second_version}_full_all'
            log_message = '###            Using Full SKY             ###'

        self.output_name = output_name
        self.log.info(' ')
        self.log.info('#' * 45)
        self.log.info(log_message)
        self.log.info('#' * 45)

        first_tables = set(first_plan_config['order'])
        second_tables = set(second_plan_config['order'])
        intersecting_tables = list(first_tables & second_tables)
        intersecting_tables.sort()
        self.intersecting_tables = intersecting_tables
        self.log.info(f'The intersecting tables to match are: {intersecting_tables}')
        self.split_insert_nunmber = config['split_insert_nunmber']

        if self.split_query:
            self.split_query_info = {el[0]: (el[1], el[2]) for el in self.split_query}

    def match_in_table(self, name, split=False, min_targetid=0, max_targetid=0):

        """This method takes as argument the name of the table to be matched to
        do the query and find the catalogid pairs. This method can also be used
        with argument split=True in which case argument min_targetid, and max_targetid are used
        to do a sub query restricted to a target_id range. For each input table ingested
        in both crossmatch runs the main query starts by creating a cte querying the target
        table and storing the target_id in the relational table associated with the catalogids
        from the target table that belong to the older crossmatch run version_id.
        Then the other query takes the target_ids from the cte and checks which catalogid’s from
        both versions are associated with that target_id. Then for each target_id from the query
        it looks for all the possible pairs of catalogids in which one catalogid comes from the
        first crossmatch run and the other comes from the second crossmatch run.
        Finally these catalogid pairs for each input table are ingested in a table with suffix all

        Parameters
        ----------
        name : str
            Name of the table crossmatched in both runs in which the catalogid pais will be looked
        split : bool
            If True the query is split into multiple subqueries restricted to ranges of targetid
            To use this for a specific table the configuration file parameter split query must
            contain as an element a list where the first element is the name of the table,
            the second one is a number larger than the maximum targetid of the table,
            and the third one is the target id bin width to be used to split the queries
        min_targetid : int
            minimum target_id to be used for this specific query
        max_targetid : int
            maximum target_id to be used for this specific query

            """

        # First we set the database options for the query
        if self.database_options:
            for param in self.database_options:
                value = self.database_options[param]
                setting_message = f'SET {param}={value!r};'
                self.database.execute_sql(setting_message)
                self.log.info(setting_message)

        temp_fields = [TempMatch.lead, TempMatch.lowest_catalogid, TempMatch.highest_catalogid]

        db_models = self.database.models

        tables = [db_models['catalogdb.' + name] for name in self.intersecting_tables]
        rel_tables = [db_models['catalogdb.catalog_to_' + name]
                      for name in self.intersecting_tables]

        index = self.intersecting_tables.index(name)
        table = tables[index]
        rel_table = rel_tables[index]
        tablename = table._meta.table_name
        cte_targetids = (Target.
                         select(rel_table.target_id, Target.ra, Target.dec)
                         .join(rel_table, on=(Target.catalogid == rel_table.catalogid))
                         .where(rel_table.version == self.first_version)
                         .distinct(rel_table.target_id)
                         .cte('cte_targets'))

        query = (rel_table
                 .select(rel_table.target_id,
                         peewee.Value(tablename).alias('table'),
                         fn.array_agg(rel_table.catalogid)
                           .over(partition_by=rel_table.target_id).alias('catalogids'),
                         fn.array_agg(rel_table.version_id)
                           .over(partition_by=rel_table.target_id).alias('version_ids'))
                 .join(cte_targetids, on=(rel_table.target_id == cte_targetids.c.target_id))
                 .where((rel_table.version == self.first_version) |
                        (rel_table.version == self.second_version))
                 .distinct(rel_table.target_id)
                 .with_cte(cte_targetids))

        # If sample_region is included in the configuration file, then the entire process
        # restricted to a (single) sample region where the tuple in the configuration file
        # indicates the ra center, the dec center, and the radius to use

        if self.sample_region:
            racen, deccen, radius = self.sample_region
            query = query.where(fn.q3c_radial_query(cte_targetids.c.ra,
                                                    cte_targetids.c.dec,
                                                    racen, deccen, radius))

        if split:
            query = query.where(rel_table.target_id >= min_targetid,
                                rel_table.target_id < max_targetid)

        if not split:
            self.log.info(' ')
            self.log.info(' ')

        # In this step we loop over the results and for each targetid we check all the pair
        # of catalogids associated to it where one comes from the first crossmatch run and
        # the other one from the second crossmatch run

        res = query.dicts()
        results_list = []
        npair = 0
        nres = 0
        for entry in res:
            nres += 1
            catids = entry['catalogids']
            verids = entry['version_ids']
            catid_combs = list(combinations(catids, 2))  # Check all the pairs of catalogid's
            verid_combs = list(combinations(verids, 2))  # Check all the pairs of version_id's
            ind_valid = [ind for ind in range(len(catid_combs))
                         if verid_combs[ind][0] != verid_combs[ind][1]]
            for ind_comb in ind_valid:  # For each pair comming from different version_id's
                npair += 1
                sorted_catids = np.sort(catid_combs[ind_comb])
                curr_tuple = (entry['table'], sorted_catids[0], sorted_catids[1])
                results_list.append(curr_tuple)
                if self.show_first and npair <= self.show_first and not split:
                    tuple_log = str(npair) + '  ' + '  '.join([str(el) for el in curr_tuple])
                    self.log.info(tuple_log)

        # In this final part we go over the results and split it in chunks of N entries
        # each time ingesting those entries to the output table

        parts = int((len(results_list) - 1) / self.split_insert_nunmber) + 1
        for npart in range(parts):
            with self.database.atomic():
                if len(results_list) == 0:
                    continue
                if not split:
                    self.log.info(f'Ingesting part {npart+1} out of {parts}')
                first_el = npart * self.split_insert_nunmber
                last_el = (npart + 1) * self.split_insert_nunmber
                curr_results = results_list[first_el:last_el]
                TempMatch.insert_many(curr_results, fields=temp_fields).execute()

        return nres, npair

    def run(self, starting_table=None):
        """Loops over the input tables to check the catalogid pairs associated to each target id.

        This method sets the models creates the tables and then goes over all the tables
        to look in each one the catalogid matches from both runs. It can be used with
        the starting_table name option to skip all the tables before that one in the process.

        Parameters
        ----------
        starting_table : str
            Name of the table in which the process would be restarted. By default tables are
            ordered alphabetically.
        """

        t0 = time.time()
        tables_to_process = self.intersecting_tables
        if starting_table:
            ind = tables_to_process.index(starting_table)
            tables_to_process = tables_to_process[ind:]

        TempMatch._meta.table_name = self.output_name
        self.database.bind([TempMatch])
        if TempMatch.table_exists():
            self.database.drop_tables([TempMatch])
            self.log.info(f'Dropped table {TempMatch._meta.table_name}')

        self.database.create_tables([TempMatch])
        self.log.info(f'Created table {TempMatch._meta.table_name}')

        for curr_table in tables_to_process:
            t1 = time.time()
            logname = 'Running table ' + curr_table
            self.log.info('-' * 52)
            self.log.info(f'|{logname:^50}|')
            self.log.info('-' * 52)

            if self.split_query and curr_table in self.split_query_info.keys():
                ntargetids, npairs = 0, 0
                self.log.info(f'Doing the query of {curr_table} splitting by targetid')
                max_tar, bin_tar = self.split_query_info[curr_table]
                self.log.info(f'Using binning={bin_tar} and max value ={max_tar}')
                min_targetid, max_targetid = 0, bin_tar
                while (min_targetid < max_tar):
                    ti_iter = time.time()
                    curr_targetids, curr_pairs = self.match_in_table(curr_table, split=True,
                                                                     min_targetid=min_targetid,
                                                                     max_targetid=max_targetid)
                    ntargetids += curr_targetids
                    npairs += curr_pairs
                    min_targetid += bin_tar
                    max_targetid += bin_tar
                    tf_iter = time.time()
                    self.log.info(f'Ingested {curr_pairs} to output table using targetid range '
                                  f'{min_targetid} to {max_targetid} '
                                  f'in {(tf_iter-ti_iter):.2f} seconds')
            else:
                ntargetids, npairs = self.match_in_table(curr_table)
            t2 = time.time()
            self.log.info(f'Found {ntargetids} targetids with multiple catalogids associated')
            self.log.info(f'{npairs} valid catalogid pairs were ingested to output table ')
            self.log.info(f'Table {curr_table} was processed in {(t2-t1):.2f} seconds')
        t3 = time.time()
        self.log.info(f'The entire match took {(t3-t0):.2f} seconds')


def create_unique_table(first_xmatch_version, second_xmatch_version):
    """Function to create the final output table with unlist of unique pairs of catalogid.

    This function takes the already created table with all the catalogid pairs for each input
    table and selects the list of unique pairs indicating for each pair the table in which it
    was found and the number of tables in which it was found

    Parameters
    ----------
    first_xmatch_version : int
        version_id number of the first crossmatch run
    second_xmatch_version : int
        version_id number of the second crossmatch run

    """
    ti = time.time()
    input_tablename = (f'catalog_ver{first_xmatch_version}_to_'
                       f'ver{second_xmatch_version}_full_all')
    output_tablename = (f'catalog_ver{first_xmatch_version}_to_'
                        f'ver{second_xmatch_version}_full_unique')
    UniqueMatch._meta.table_name = output_tablename
    database.bind([UniqueMatch])
    log = target_selection.log
    log_path_and_name = os.path.realpath(f'./catalog_ver{first_xmatch_version}'
                                         f'_to_ver{second_xmatch_version}'
                                         '_unique_table_creation.log')
    log.start_file_logger(log_path_and_name, rotating=False, mode='a')

    if UniqueMatch.table_exists():
        database.drop_tables([UniqueMatch])
        log.info(f'Dropped table {output_tablename}')
    database.create_tables([UniqueMatch])
    log.info(f'Created table {output_tablename}')
    TempMatch._meta.table_name = input_tablename
    query = (TempMatch
             .select(TempMatch.lowest_catalogid,
                     TempMatch.highest_catalogid,
                     fn.array_agg(TempMatch.lead),
                     fn.count(TempMatch.lead))
             .join(Target, on=(TempMatch.lowest_catalogid == Target.catalogid))
             .group_by(TempMatch.lowest_catalogid,
                       TempMatch.highest_catalogid))

    fields = [UniqueMatch.lowest_catalogid, UniqueMatch.highest_catalogid,
              UniqueMatch.match_tables, UniqueMatch.n_tables]

    insert_query = UniqueMatch.insert_from(query, fields).returning()
    n_unique = insert_query.execute()
    tf = time.time()
    log.info(f'Created unique pairs table with {n_unique} '
             f'entries in {(tf-ti):.2f} seconds')


"""
EXAMPLE OF USAGE
metax = MetaXMatch('cat_to_cat.yml', database)
metax.run()
create_unique_table(25,31)
"""
